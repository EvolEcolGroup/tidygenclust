// The quadratic program of the Newton step of one column (see
// radmix_qp_sweep.h for the problem), solved exactly by an active-set
// method based on the sweep operator (Jennrich and Sampson 1978; Lange,
// "Numerical Analysis for Statisticians"). The solver follows K. Lange's
// sweep-operator code as implemented in OpenADMIXTURE.jl (MIT licence, see
// inst/COPYRIGHTS).
//
// Why a small dense solver: K, the number of ancestral populations, is
// small (typically below 20), while each block update solves one problem per
// individual or per locus. Each problem costs at most O(K^3), fits in cache,
// and is solved by one thread with preallocated scratch memory (QPWork) and
// no library calls; a general-purpose solver would spend longer setting up
// than computing.
//
// The method works on a "tableau", a symmetric matrix that bundles the
// Hessian, the constraint and the gradient. Sweeping a diagonal entry brings
// the corresponding variable into the Newton system, and the inverse sweep
// takes it out again, each at O(K^2) cost:
//   - after sweeping every parameter (and the sum constraint), the last
//     column of the tableau holds the Newton step of the problem with the
//     equality constraint but without the bounds;
//   - the solver walks along that step until a parameter reaches a bound,
//     un-sweeps (fixes) that parameter, recomputes the step of the others,
//     and so on; a fixed parameter whose gradient points back into the box
//     is swept in again. This is an active-set strategy, in which "swept"
//     means free and "unswept" means held on a bound.

#include "radmix_qp_sweep.h"
#include <algorithm>
#include <cmath>

namespace radmix {

// (the scratch vectors are sized for the tableau of a column of q, with
// K + 2 rows; that of a column of p, with K + 1 rows, fits in them)
QPWork::QPWork(int K) : K(K), tab((K + 2) * (K + 2)), diag(K + 2),
  swept(K + 2), tmp(K + 2), delta(K) {}

// Sweep (or inverse sweep) the symmetric n x n matrix a on diagonal k.
//
// The sweep operator is Gauss-Jordan elimination on one pivot, arranged to
// keep the matrix symmetric. With pivot a_kk, a sweep replaces
//   a_kk -> -1 / a_kk
//   a_ik -> a_ik / a_kk,  a_ki -> a_ki / a_kk          (i != k)
//   a_ij -> a_ij - a_ik a_kj / a_kk                    (i, j != k)
// and the inverse sweep is the same except that the off-diagonal entries of
// row and column k change sign (a_ik -> -a_ik / a_kk). An inverse sweep
// undoes a sweep, and sweeps on different pivots commute.
//
// If the matrix is partitioned as [[A, B], [B', D]] and every pivot of the
// A block is swept, the result is
//   [[-A^-1, A^-1 B], [B' A^-1, D - B' A^-1 B]]:
// the A block is inverted (with a minus sign), B is solved against A, and D
// becomes the Schur complement. The solver uses this to get the Newton step
// of any subset of free variables. tmp (n values) is scratch.
static void sweep(double* a, int n, int k, double* tmp, bool inverse) {
  const double p = 1.0 / a[k + k * n];
  // save column k, then zero row and column k, so that the rank-one update
  // below produces the new row and column k directly
  for (int i = 0; i < n; i++) {
    tmp[i] = a[i + k * n];
    a[i + k * n] = 0;
    a[k + i * n] = 0;
  }
  // With tmp[k] = -1 (sweep), a_ij -= p tmp[i] tmp[j] gives:
  //   i, j != k : a_ij - a_ik a_jk / a_kk      (Schur complement)
  //   i or j = k: 0 + a_ik / a_kk               (scaled row and column k)
  //   i = j = k : 0 - 1 / a_kk                  (negated reciprocal)
  // With tmp[k] = +1 (inverse sweep) the middle case becomes -a_ik / a_kk.
  tmp[k] = inverse ? 1.0 : -1.0;
  for (int i = 0; i < n; i++) {
    const double pv = p * tmp[i];
    double* col = a + i * n;
    for (int j = 0; j < n; j++) col[j] -= pv * tmp[j];
  }
}

// Build the tableau of one problem; returns its dimension n.
//
// For a column of p (box constraints only), n = K + 1:
//
//   [ H      grad ]      K rows
//   [ grad'    0  ]      1 row
//
// Sweeping the K diagonal entries of H turns the last column into H^-1 grad,
// the unconstrained Newton step.
//
// For a column of q (box plus the sum constraint), n = K + 2:
//
//   [ H + mu 11'   1   grad ]   K rows    (parameters)
//   [ 1'           0   r    ]   1 row     (constraint: 1' delta = r)
//   [ grad'        r   0    ]   1 row     (right-hand side)
//
// where r = 1 - sum(par) is the change in the total that delta must make
// (0 when par already sums to one). Write A = H + mu 11'. Sweeping the K
// parameters puts A^-1 grad in the last column and -1' A^-1 1 (negative when
// A is positive definite) in the diagonal entry of the constraint; sweeping
// the constraint as well turns the last column into
//   delta = A^-1 (grad + lambda 1),
//   lambda = (r - 1' A^-1 grad) / (1' A^-1 1),
// the solution of the equality-constrained problem by Lagrange multipliers
// (lambda is the multiplier of the sum constraint, which makes 1' delta = r).
//
// Why mu: sweeping the parameters inverts the K x K block, which must
// therefore be positive definite, whereas the constrained problem only
// needs H to be positive definite on the subspace 1' delta = 0. H is a
// negative Hessian of the log-likelihood, a sum of terms d2 v v' with
// d2 >= 0 (radmix_kernels.cpp), so it is positive semi-definite, and then
// H + mu 11' is positive definite for any mu > 0 as soon as H is positive
// definite on the subspace: for x with 1' x != 0, x' (H + mu 11') x >=
// mu (1' x)^2 > 0. (H itself can be singular off the subspace, e.g. when the
// frequencies of two populations are proportional at every locus.) The
// shift does not change the solution: every feasible delta has 1' delta = r,
// so it only adds the constant mu r^2 / 2 to the objective. We take
// mu = trace(H) / K^2, so that mu 11', whose only non-zero eigenvalue is
// mu K, adds the mean eigenvalue of H along the direction of the ones, which
// keeps A as well conditioned as H.
static int create_tableau(const double* H, const double* grad,
                          const double* par, int K, bool simplex,
                          QPWork& work) {
  const int n = simplex ? K + 2 : K + 1;
  double* tab = work.tab.data();
  // the index of the right-hand side (the gradient)
  const int last = n - 1;
  // (this also sets the diagonal entry of the constraint and the
  // bottom-right corner to 0)
  for (int i = 0; i < n * n; i++) tab[i] = 0;
  // H into the top-left K x K block (leading dimension n)
  for (int j = 0; j < K; j++)
    for (int i = 0; i < K; i++) tab[i + j * n] = H[i + j * K];
  // the gradient into the last row and column
  for (int i = 0; i < K; i++) {
    tab[last + i * n] = grad[i];
    tab[i + last * n] = grad[i];
  }
  if (!simplex) return n;
  // mu = trace(H) / K^2 (floored at 0, in case H is not a valid negative
  // Hessian, e.g. in the tests)
  double trace = 0.0;
  for (int i = 0; i < K; i++) trace += H[i + i * K];
  const double mu = std::max(0.0, trace) / ((double) K * K);
  for (int j = 0; j < K; j++)
    for (int i = 0; i < K; i++) tab[i + j * n] += mu;
  // the constraint row and column (index K): the coefficients of 1' delta
  for (int i = 0; i < K; i++) {
    tab[K + i * n] = 1.0;
    tab[i + K * n] = 1.0;
  }
  // its right-hand side, r = 1 - sum(par)
  double sum_par = 0;
  for (int i = 0; i < K; i++) sum_par += par[i];
  tab[K + last * n] = tab[last + K * n] = 1.0 - sum_par;
  return n;
}

// Solve one problem (see radmix_qp_sweep.h) by the active-set sweep method.
//
// Invariant of the loop, for the current set F of free (swept) parameters:
//   - for i in F, tab[i, last] is the remaining Newton step of parameter i,
//     how far it still has to move to reach the minimiser of the quadratic
//     with the fixed parameters held where they are;
//   - for i not in F, tab[i, last] is the gradient of the quadratic model
//     with respect to parameter i at that minimiser (positive: increasing
//     par[i] would still improve the log-likelihood). A step along the
//     Newton direction does not change these entries, which is why the step
//     below only rescales the rows of the free parameters.
// For q, where the sum constraint is swept too, "gradient" means that of the
// Lagrangian, grad - A delta + lambda 1: 0 for a free parameter, and for a
// fixed one the rate at which the model improves when par[i] increases and
// the free parameters absorb the change in the sum. (mu changes every entry
// of A delta by the same amount, mu r, which lambda absorbs, so it does not
// change these signs either.)
int qp_solve_column(const double* H, const double* grad, const double* par,
                    int K, bool simplex, double* delta, QPWork& work) {
  // the bounds of par + delta (the projections then tighten them to
  // [1e-5, 1 - 1e-5])
  const double pmin = 0.0, pmax = 1.0;
  // small: a parameter within 1e-5 of a bound counts as on the bound (the
  // floor of the projections, so that a value left at the floor by the
  // previous update is treated as being on the bound). tol: the relative
  // size below which a pivot means that the matrix is not numerically
  // positive definite.
  const double small = 1e-5, tol = 1e-8;
  // p parameters and c equality constraints (0 or 1)
  const int p = K, c = simplex ? 1 : 0;
  const int n = create_tableau(H, grad, par, K, simplex, work);
  const int last = n - 1;
  double* tab = work.tab.data();
  double* d = work.diag.data();
  double* tmp = work.tmp.data();
  char* swept = work.swept.data();
  // delta accumulates the steps taken
  for (int i = 0; i < K; i++) delta[i] = 0;
  work.passes = 0;
  // Sweep the parameters, then the constraint; give up (delta = 0) if the
  // matrix is not numerically positive definite. Sweeping the pivots one
  // after the other is a Cholesky-like elimination: just before pivot i is
  // swept, tab[i, i] is the i-th pivot of the factorisation, positive for
  // every i exactly when the matrix is positive definite. Requiring it to
  // stay above tol times its original value also rejects matrices that are
  // positive definite only in exact arithmetic (a pivot that has shrunk by
  // eight orders of magnitude would make the step mostly rounding error).
  for (int i = 0; i < n; i++) d[i] = tab[i + i * n];
  for (int i = 0; i < p; i++) {
    if (d[i] <= 0 || tab[i + i * n] < d[i] * tol) return 0;
    sweep(tab, n, i, tmp, false);
  }
  // all the parameters start free
  for (int i = 0; i < p; i++) swept[i] = 1;
  // After the parameter sweeps, the diagonal entry of the constraint is
  // -1' A^-1 1, which must be negative; sweeping it imposes 1' delta = r.
  // (For p, c = 0 and this loop does nothing.)
  for (int i = p; i < p + c; i++) {
    if (tab[i + i * n] >= 0) return 0;
    sweep(tab, n, i, tmp, false);
  }
  // The active-set loop. Each pass takes a (possibly partial) step, then
  // either fixes a free parameter that has reached a bound, frees a fixed
  // parameter whose gradient points into the box, or stops at a KKT point.
  // `freed` is the parameter freed in the previous pass (-1 if none), used
  // to detect cycling at a vertex (see below).
  int freed = -1;
  for (int iteration = 1; iteration <= 1000; iteration++) {
    work.passes = iteration;
    const int freed_last_pass = freed;
    freed = -1;
    // The longest step a in [0, 1] along the remaining Newton step
    // u = tab[F, last] for which every free parameter stays within [0, 1]
    // (a parameter moving up is limited by the upper bound, one moving down
    // by the lower bound; a = 1 is the full step). Components with
    // |u_i| <= 1e-10 cannot reach a bound and are ignored. (a can come out
    // slightly negative, a small step back, when rounding has left a free
    // parameter just beyond its bound, e.g. when two parameters reached
    // their bounds in the same step and only one has been fixed so far.)
    double a = 1.0;
    for (int i = 0; i < p; i++) {
      if (!swept[i]) continue;
      const double ui = tab[i + last * n];
      const double ai = (ui > 0 ? pmax : pmin) - par[i] - delta[i];
      if (std::fabs(ui) > 1e-10) a = std::min(a, ai / ui);
    }
    // Take the step a u for the free parameters; only the fraction (1 - a)
    // of the Newton step remains. The gradient of the quadratic is linear in
    // delta, so rescaling the remaining step keeps the tableau consistent
    // without sweeping again. The fixed parameters and, for q, the entry of
    // the constraint (-lambda) describe the minimiser itself, which a step
    // towards it does not change.
    for (int i = 0; i < p; i++) {
      if (!swept[i]) continue;
      const double ui = tab[i + last * n];
      delta[i] += a * ui;
      tab[i + last * n] = (1.0 - a) * ui;
      tab[last + i * n] = tab[i + last * n];
    }
    // Fix a free parameter that has reached a bound (within `small`): the
    // inverse sweep removes it from the Newton system, and the last column
    // becomes the Newton step of the remaining free parameters (still with
    // the sum constraint for q, which stays swept). One parameter per pass.
    // The |tab[i, i]| > 1e-10 test avoids dividing by a vanishing pivot; for
    // q, it also keeps the last free parameter free, as its value is then
    // fixed by the others and its pivot is 0 up to rounding.
    bool cycle = false;
    for (int i = 0; i < p; i++) {
      const double est = par[i] + delta[i];
      const bool critical = pmin >= est - small || pmax <= est + small;
      if (swept[i] && std::fabs(tab[i + i * n]) > 1e-10 && critical) {
        // Cycling at a vertex of the simplex: the parameter freed in the
        // previous pass could not move (a step of exactly 0), and would be
        // fixed again, then freed again, and so on. delta is the solution:
        // stop (see the note at the end).
        if (i == freed_last_pass && a == 0.0) return 0;
        sweep(tab, n, i, tmp, true);
        swept[i] = 0;
        cycle = true;
        break;
      }
    }
    if (cycle) continue;
    // Free a fixed parameter that violates the KKT conditions. At the
    // optimum, each fixed parameter must be held by its bound: its gradient
    // must point out of the box (u_i <= 0 at the lower bound, u_i >= 0 at
    // the upper bound). If it points into the box, the objective can still
    // be improved by moving that parameter, which is swept back into the
    // Newton system. A parameter whose pivot has vanished (no curvature
    // left once the free parameters are accounted for) cannot usefully
    // re-enter the system, and sweeping it would divide by zero: it stays
    // on its bound.
    for (int i = 0; i < p; i++) {
      const double ui = tab[i + last * n];
      const double est = par[i] + delta[i];
      const bool violation = (ui > 0 && pmin >= est - small) ||
        (ui < 0 && pmax <= est + small);
      if (!swept[i] && std::fabs(tab[i + i * n]) > 1e-10 && violation) {
        sweep(tab, n, i, tmp, false);
        swept[i] = 1;
        freed = i;
        cycle = true;
        break;
      }
    }
    if (cycle) continue;
    // no free parameter on a bound, and no fixed parameter that wants to
    // leave its bound: the KKT conditions hold, and delta is the solution
    return iteration;
  }
  // The iteration limit. A solution at a vertex of the simplex, where every
  // parameter is on a bound, makes the active set cycle: once the sum
  // constraint is swept, the last free parameter is determined by the
  // others, its pivot is 0, and the guard above keeps it free even on a
  // bound; a fixed parameter whose gradient points into the box is then
  // freed, takes a step of 0 (the last free parameter cannot move past its
  // bound), is fixed again, and so on, with delta unchanged. This is
  // detected above and stops the loop at once (it would otherwise take 1000
  // passes, some 200 times the cost of a typical column). If the steps are
  // not exactly 0, the loop can still reach the limit, with delta usually at
  // the solution already, which is why the kernels use delta regardless of
  // the return value.
  return 0;
}

}  // namespace radmix

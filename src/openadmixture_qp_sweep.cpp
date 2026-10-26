// Port of OpenADMIXTURE.jl src/qp.jl (create_tableau!, quadratic_program!,
// sweep!), originally written by K. Lange and modified by S. Ko. Only the
// two cases used by update_q!/update_p! are ported: box constraints alone
// (P) and box plus the single sum constraint (Q); the p_double variant is
// not needed. Bounds are fixed at pmin = 0, pmax = 1 as in the callers.
//
// Overview
// --------
// Each Newton step of ADMIXTURE updates Q one individual at a time and P one
// locus at a time, holding everything else fixed. For a column with current
// value par, log-likelihood gradient grad and negative Hessian H, the step
// delta minimises the local quadratic model
//   0.5 * delta' H delta - grad' delta
// subject to 0 <= par + delta <= 1 (and sum(delta) = 1 - sum(par) for Q, so
// that the new ancestry proportions still sum to one). See
// openadmixture_qp_sweep.h for the callers in src/openadmixture_kernels.cpp.
//
// Why a small dense active-set solver: K (the number of ancestral
// populations) is small, typically below 20, while there are I + J
// independent problems per outer iteration (one per individual, one per
// locus). Each problem therefore costs O(K^3) at most, fits in cache, and
// can be solved by one thread with preallocated scratch (QPWork) and no
// library calls. A general-purpose QP solver would spend far longer on
// set-up than on the arithmetic.
//
// The method (Jennrich and Sampson 1978; Lange, "Numerical Analysis for
// Statisticians") works on a "tableau": a symmetric matrix that bundles the
// Hessian, the constraint and the gradient. The sweep operator is applied
// to diagonal entries of the tableau to bring variables into, or take them
// out of, the Newton system:
//   * after sweeping every parameter (and the sum constraint), the last
//     column of the tableau holds the Newton step for the equality
//     constrained problem, ignoring the box bounds;
//   * the solver then walks along that step until a parameter hits a bound,
//     un-sweeps (fixes) that parameter, recomputes the step for the
//     remaining free parameters, and so on. Fixed parameters whose
//     gradient points back into the feasible box are swept in again. This
//     is a classic active-set strategy: "swept" = free, "unswept" = held on
//     a bound.
#include "openadmixture_qp_sweep.h"
#include <algorithm>
#include <cmath>

// internal code of the OpenADMIXTURE implementation (see the note on the
// namespace in openadmixture_geno_access.h)
namespace openadmixture {


// V: an orthonormal basis of R^K whose first column is ones(K) / sqrt(K),
// playing the role of svd(ones(1, K), full = true).V in Julia. We use the
// Householder reflector I - 2 u u' / (u' u), u = e1 - ones / sqrt(K), which
// maps e1 to ones / sqrt(K). It is symmetric, as is the reflector LAPACK
// returns, so V and V' are interchangeable in the mu formula below.
//
// The scratch vectors are sized for the simplex (Q) tableau, which has
// K + 2 rows; the P tableau (K + 1 rows) fits inside the same buffers.
QPWork::QPWork(int K) : K(K), tab((K + 2) * (K + 2)), diag(K + 2),
  swept(K + 2), V(K * K), tmp(K + 2) {
  const double s = 1.0 / std::sqrt((double) K);
  // u = e1 - ones / sqrt(K)
  std::vector<double> u(K, -s);
  u[0] += 1.0;
  double uu = 0;
  for (int k = 0; k < K; k++) uu += u[k] * u[k];
  // uu == 0 only for K == 1, where V = 1
  const double scale = uu > 0 ? 2.0 / uu : 0.0;
  // V = I - (2 / u'u) u u', filled column by column (column-major)
  for (int j = 0; j < K; j++)
    for (int i = 0; i < K; i++)
      V[i + j * K] = (i == j ? 1.0 : 0.0) - scale * u[i] * u[j];
}

// Sweep (or inverse sweep) the symmetric n x n tableau a on diagonal k.
//
// The sweep operator is Gauss-Jordan elimination on one pivot, arranged so
// that a symmetric matrix stays symmetric. Writing a_kk for the pivot,
// a forward sweep on k replaces
//   a_kk -> -1 / a_kk
//   a_ik -> a_ik / a_kk,  a_ki -> a_ki / a_kk          (i != k)
//   a_ij -> a_ij - a_ik a_kj / a_kk                    (i, j != k)
// and the inverse sweep is identical except that the off-diagonal entries
// of row/column k change sign: a_ik -> -a_ik / a_kk. Sweeping and then
// inverse-sweeping the same k restores the original matrix, and sweeps on
// different pivots commute.
//
// What this achieves: if the tableau is partitioned as [[A, B], [B', D]]
// and every pivot in the A block is swept, the result is
//   [[-A^{-1}, A^{-1} B], [B' A^{-1}, D - B' A^{-1} B]],
// i.e. the A block is inverted (with a minus sign), the off-diagonal block
// is "solved" against A, and the rest becomes the Schur complement. The
// solver uses this to obtain Newton steps for any chosen subset of free
// variables: sweeping a variable adds it to the linear system, inverse
// sweeping removes it again, each at O(n^2) cost.
//
// tmp (length n) is scratch for the pivot column.
static void sweep(double* a, int n, int k, double* tmp, bool inverse) {
  const double p = 1.0 / a[k + k * n];
  // save column k, then zero row and column k so that the rank-one update
  // below produces the new row/column k directly
  for (int i = 0; i < n; i++) {
    tmp[i] = a[i + k * n];
    a[i + k * n] = 0;
    a[k + i * n] = 0;
  }
  // With tmp[k] = -1 (forward) the update a_ij -= p * tmp[i] * tmp[j]
  // gives, case by case:
  //   i, j != k : a_ij - a_ik a_jk / a_kk      (Schur complement)
  //   i or j = k: 0 + a_ik / a_kk               (scaled pivot row/column)
  //   i = j = k : 0 - 1 / a_kk                  (negated reciprocal)
  // With tmp[k] = +1 (inverse) the middle case becomes -a_ik / a_kk, which
  // is exactly the sign change that undoes a forward sweep.
  tmp[k] = inverse ? 1.0 : -1.0;
  for (int i = 0; i < n; i++) {
    const double pv = p * tmp[i];
    double* col = a + i * n;
    for (int j = 0; j < n; j++) col[j] -= pv * tmp[j];
  }
}

// Build the tableau as create_tableau! does. Julia is called with
// r = -grad (the caller negates Xtz, the log-likelihood gradient), and
// stores -r, so the last row/column of the tableau holds grad itself.
//
// Layout for P (box constraints only), n = K + 1:
//
//   [ H     grad ]      K rows
//   [ grad'   0  ]      1 row
//
// Sweeping the K diagonal entries of H turns this into
//
//   [ -H^{-1}            H^{-1} grad        ]
//   [ grad' H^{-1}   -grad' H^{-1} grad     ]
//
// so the last column now holds H^{-1} grad, the unconstrained Newton step
// (the minimiser of 0.5 delta' H delta - grad' delta).
//
// Layout for Q (box plus the sum constraint), n = K + 2:
//
//   [ H + mu 11'   1   grad ]   K rows    (parameters)
//   [ 1'           0   r    ]   1 row     (constraint: 1' delta = r)
//   [ grad'        r   0    ]   1 row     (right-hand side)
//
// where r = 1 - sum(par) is the change in the total that delta must make
// (0 when par already sums to one). Write A = H + mu 11'. Sweeping the K
// parameter pivots gives the unconstrained step A^{-1} grad in the last
// column and, in the constraint's diagonal entry, -1' A^{-1} 1, which is
// negative whenever A is positive definite. Sweeping that constraint pivot
// as well then turns the last column into
//   delta = A^{-1} (grad + lambda 1),
//   lambda = (r - 1' A^{-1} grad) / (1' A^{-1} 1),
// which is the solution of the equality-constrained problem by Lagrange
// multipliers: lambda is the multiplier of the sum constraint, chosen so
// that 1' delta = r exactly.
//
// Returns the tableau dimension n.
static int create_tableau(const double* H, const double* grad,
                          const double* par, int K, bool simplex,
                          QPWork& work) {
  const int n = simplex ? K + 2 : K + 1;
  double* tab = work.tab.data();
  // index of the right-hand-side row/column (the gradient)
  const int last = n - 1;
  // start from zero: this also sets the constraint diagonal and the
  // bottom-right corner to 0
  for (int i = 0; i < n * n; i++) tab[i] = 0;
  // copy H into the top-left K x K block (note leading dimension n, not K)
  for (int j = 0; j < K; j++)
    for (int i = 0; i < K; i++) tab[i + j * n] = H[i + j * K];
  // gradient into the last row and last column (kept symmetric)
  for (int i = 0; i < K; i++) {
    tab[last + i * n] = grad[i];
    tab[i + last * n] = grad[i];
  }
  if (!simplex) return n;
  // mu from the Gerschgorin circle theorem so that H + mu * 11' is
  // positive definite: w = V * H * V[:, 1], mu = 2 * (|w|_1 - 2 w_1) / K.
  // Deviation: Julia does not floor mu at 0, so for a well-conditioned H
  // mu can be negative and make the tableau indefinite (the original
  // Lange code, kept commented out in qp.jl, floors it at 0 and so do we).
  // In practice Julia's mu is 0 anyway: mul!(tmp_k, v, tmp_k) aliases its
  // input and output, and BLAS gemv with beta = 0 zeroes y (= x) first.
  // On the feasible set mu * (1' delta)^2 is constant, so mu only affects
  // whether the initial sweeps succeed, never the solution.
  //
  // Why the penalty is needed: the sweeps on the parameter block invert
  // the full K x K matrix, so it must be positive definite on all of R^K.
  // But the constrained problem only needs H to be positive definite on the
  // subspace 1' delta = 0; along the ones direction H may be singular or
  // nearly so. Adding mu 11' stiffens only that direction. Since every
  // feasible delta has 1' delta = r, the objective changes by the constant
  // 0.5 * mu * r^2 there, so the minimiser is the same.
  //
  // Where the formula comes from: in the basis V (first vector
  // ones / sqrt(K)), V' (mu 11') V = mu K e1 e1', so the penalty adds mu K
  // to the (1, 1) entry of V' H V and nothing else. w is the first column
  // of V' H V. Gerschgorin's theorem places an eigenvalue within
  // sum_{j > 1} |w_j| of the diagonal entry w_1 + mu K, so choosing
  // mu K >= |w|_1 - 2 w_1 (= sum_{j > 1} |w_j| - w_1 when w_1 >= 0)
  // keeps that disc out of the negative half-line; the factor 2 is a
  // safety margin inherited from Lange's code. When H is already
  // comfortably diagonally dominant in that direction the formula goes
  // negative, which would subtract curvature, hence the floor at 0.
  const double* V = work.V.data();
  // Hv = H * V[:, 1] (the first column of V is ones / sqrt(K))
  double* Hv = work.tmp.data();
  for (int i = 0; i < K; i++) {
    double s = 0;
    for (int j = 0; j < K; j++) s += H[i + j * K] * V[j];
    Hv[i] = s;
  }
  // w = V * Hv, accumulated on the fly: we only need |w|_1 and w_1
  double norm1 = 0, w1 = 0;
  for (int i = 0; i < K; i++) {
    double s = 0;
    for (int j = 0; j < K; j++) s += V[i + j * K] * Hv[j];
    norm1 += std::fabs(s);
    if (i == 0) w1 = s;
  }
  const double mu = std::max(0.0, 2.0 * (norm1 - 2.0 * w1) / K);
  // add mu * 11' to the parameter block
  for (int j = 0; j < K; j++)
    for (int i = 0; i < K; i++) tab[i + j * n] += mu;
  // constraint row/column (index K): the coefficients of 1' delta
  for (int i = 0; i < K; i++) {
    tab[K + i * n] = 1.0;
    tab[i + K * n] = 1.0;
  }
  // right-hand side of the constraint, r = 1 - sum(par), placed where the
  // constraint row meets the last column (and symmetrically)
  double sum_par = 0;
  for (int i = 0; i < K; i++) sum_par += par[i];
  tab[K + last * n] = tab[last + K * n] = 1.0 - sum_par;
  return n;
}

// Solve one column QP (problem statement in openadmixture_qp_sweep.h) by the
// active-set sweep method. Inputs: H (K x K negative Hessian, column-major),
// grad and par (length K), simplex (true for Q). Output: delta (length K).
// Returns the number of active-set iterations, or 0 if the initial sweeps fail
// (H not numerically positive definite, delta = 0) or if the loop does not
// settle within 1000 iterations (delta as far as it got).
//
// Invariant maintained by the loop, for the current free set F (swept
// parameters):
//   * for i in F, tab[i, last] is the remaining Newton step for parameter
//     i, i.e. how far it still has to move to reach the minimiser of the
//     quadratic with the fixed parameters held where they are;
//   * for i not in F, tab[i, last] is the gradient of the quadratic model
//     with respect to parameter i at that minimiser (positive means
//     increasing par[i] would still improve the log-likelihood). These
//     entries are not changed by a step along the Newton direction, which
//     is why the step loop below only rescales the swept rows.
// For Q, where the sum constraint is swept too, "gradient" means the
// gradient of the Lagrangian: with A = H + mu 11' and lambda as in
// create_tableau, tab[i, last] = (grad - A delta + lambda 1)_i at the
// minimiser. It is 0 for every free parameter; for a fixed one it is the
// rate at which the model improves when par[i] increases and the free
// parameters absorb the change in the sum. mu changes every entry of
// A delta by the same amount (mu r, as 1' delta = r), which lambda absorbs,
// so mu does not affect these signs either.
int qp_solve_column(const double* H, const double* grad, const double* par,
                    int K, bool simplex, double* delta, QPWork& work) {
  // box bounds for par + delta; the projections later tighten these to
  // [1e-5, 1 - 1e-5]
  const double pmin = 0.0, pmax = 1.0;
  // small: a parameter within 1e-5 of a bound counts as on the bound (this
  // matches the 1e-5 floor applied by project_q / project_p, so a column
  // value sitting at the output floor is treated as being at the bound).
  // tol: relative pivot size below which H is declared not positive
  // definite.
  const double small = 1e-5, tol = 1e-8;
  // p: number of parameters, c: number of equality constraints (0 or 1)
  const int p = K, c = simplex ? 1 : 0;
  const int n = create_tableau(H, grad, par, K, simplex, work);
  const int last = n - 1;
  double* tab = work.tab.data();
  double* d = work.diag.data();
  double* tmp = work.tmp.data();
  char* swept = work.swept.data();
  // delta accumulates the steps taken so far
  for (int i = 0; i < K; i++) delta[i] = 0;
  // sweep on the parameters, then on the constraint; bail out (delta = 0)
  // if the Hessian is not numerically positive definite
  //
  // Sweeping pivots one after another is a Cholesky-like elimination: just
  // before pivot i is swept, tab[i, i] is the i-th pivot of the
  // factorisation, which is positive for every i exactly when the matrix is
  // positive definite. Requiring it to stay above tol times its original
  // value also rejects matrices that are positive definite only in exact
  // arithmetic (a pivot that has collapsed by eight orders of magnitude
  // means the step would be dominated by rounding error).
  for (int i = 0; i < n; i++) d[i] = tab[i + i * n];
  for (int i = 0; i < p; i++) {
    if (d[i] <= 0 || tab[i + i * n] < d[i] * tol) return 0;
    sweep(tab, n, i, tmp, false);
  }
  // all parameters start free
  for (int i = 0; i < p; i++) swept[i] = 1;
  // After the parameter sweeps the constraint diagonal holds
  // -1' A^{-1} 1, which must be strictly negative; a non-negative value
  // signals a numerically broken A. Sweeping it imposes 1' delta = r (see
  // create_tableau). For P, c = 0 and this loop does nothing.
  for (int i = p; i < p + c; i++) {
    if (tab[i + i * n] >= 0) return 0;
    sweep(tab, n, i, tmp, false);
  }
  // Active-set loop. Each pass either takes a (possibly partial) step,
  // fixes one parameter that has reached a bound, frees one fixed parameter
  // whose gradient points into the box, or stops at a KKT point. `freed` is
  // the parameter freed by the KKT test in the previous pass (-1 if none),
  // used to detect cycling at a vertex (see the note at the end).
  int freed = -1;
  work.passes = 0;
  for (int iteration = 1; iteration <= 1000; iteration++) {
    work.passes = iteration;
    const int freed_last_pass = freed;
    freed = -1;
    // maximal fractional step along the swept directions
    //
    // u = tab[F, last] is the remaining Newton step. Find the largest
    // a in [0, 1] such that par + delta + a * u stays inside [0, 1] for
    // every free parameter: a parameter moving up (u_i > 0) is limited by
    // the upper bound, one moving down by the lower bound. a = 1 means the
    // full step is feasible. Components with |u_i| <= 1e-10 are ignored,
    // as they cannot reach a bound and would make ai / ui ill-defined.
    // (a can come out marginally negative, a small step backwards, when
    // rounding has left a free parameter just beyond its bound, e.g. when
    // two parameters reached their bounds in the same step and only one of
    // them has been fixed so far.)
    double a = 1.0;
    for (int i = 0; i < p; i++) {
      if (!swept[i]) continue;
      const double ui = tab[i + last * n];
      const double ai = (ui > 0 ? pmax : pmin) - par[i] - delta[i];
      if (std::fabs(ui) > 1e-10) a = std::min(a, ai / ui);
    }
    // Take the step a * u for the free parameters and record that only
    // the fraction (1 - a) of the Newton step is still outstanding. Since
    // the quadratic's gradient is linear in delta, scaling the remaining
    // step by (1 - a) keeps the tableau consistent without re-sweeping.
    // Fixed parameters do not move (their delta stays at its bound value).
    // For Q, the constraint row (tab[K, last], which holds -lambda) is not
    // rescaled either: like the entries of the fixed parameters, it
    // describes the minimiser itself, which a step towards it leaves
    // unchanged.
    for (int i = 0; i < p; i++) {
      if (!swept[i]) continue;
      const double ui = tab[i + last * n];
      delta[i] += a * ui;
      tab[i + last * n] = (1.0 - a) * ui;
      tab[last + i * n] = tab[i + last * n];
    }
    // inverse sweep a swept parameter that has hit a bound
    //
    // A free parameter now within `small` of 0 or 1 is fixed there:
    // un-sweeping removes it from the Newton system, and the last column
    // is updated to the Newton step for the remaining free parameters
    // (still honouring the sum constraint for Q, because the constraint
    // pivot stays swept). Only one parameter is released per pass, after
    // which the step is recomputed. The |tab[i, i]| > 1e-10 test avoids
    // dividing by a vanishing pivot. (For Q this also protects the last
    // free parameter: with the constraint swept its value is fixed by the
    // others and its pivot is 0 up to rounding, so it is not released; see
    // the note on cycling at the end of this function.)
    bool cycle = false;
    for (int i = 0; i < p; i++) {
      const double est = par[i] + delta[i];
      const bool critical = pmin >= est - small || pmax <= est + small;
      if (swept[i] && std::fabs(tab[i + i * n]) > 1e-10 && critical) {
        // Cycling at a vertex: the parameter freed in the previous pass
        // could not move (a step of exactly 0) and would be fixed again,
        // after which the KKT test would free it again, and so on. delta
        // is the solution, which was otherwise returned (with 0) after
        // the 1000 passes, unchanged since this one.
        if (i == freed_last_pass && a == 0.0) return 0;
        sweep(tab, n, i, tmp, true);
        swept[i] = 0;
        cycle = true;
        break;
      }
    }
    if (cycle) continue;
    // sweep an unswept parameter that violates the KKT conditions
    //
    // At the optimum each fixed parameter must be held by its bound: the
    // gradient must push it outwards (u_i <= 0 at the lower bound, u_i >= 0
    // at the upper bound). If instead the gradient points into the box
    // (u_i > 0 at the lower bound, u_i < 0 at the upper bound), the
    // objective can still be improved by freeing that parameter, so it is
    // swept back into the Newton system and the loop continues.
    //
    // Deviation from OpenADMIXTURE.jl: Julia sweeps the violating parameter
    // without checking its pivot tab[i, i], whereas the inverse sweep above
    // is guarded by |tab[i, i]| > 1e-10. A (near) zero pivot here would
    // divide by zero inside sweep() and fill the tableau, and therefore
    // delta, with Inf/NaN. We apply the same 1e-10 guard: a parameter whose
    // pivot has vanished cannot usefully re-enter the Newton system (its
    // curvature, once the other free parameters are accounted for, is zero),
    // so it is left on its bound. When the pivot is healthy, which is
    // always the case in our tests, the behaviour is identical to Julia's.
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
    // No free parameter on a bound and no fixed parameter that wants to
    // leave its bound: the KKT conditions hold and delta is the solution.
    return iteration;
  }
  // Iteration limit reached. At a vertex solution, a parameter can be
  // freed by the KKT test, take no step, be fixed again by the critical
  // test, and so on: this is detected above, and stops the loop at once
  // with the same delta (which used to take the 1000 passes, about 4% of
  // the columns of Q near convergence, making update_q 6-21% slower). If
  // the step is not exactly 0, the loop can still reach the limit; delta
  // has usually already converged to the right step by then, so callers
  // must not treat 0 as fatal; the kernels use delta regardless and rely
  // on the projections to restore feasibility.
  //
  // Why it cycles (traced on the vertex example of
  // test_openadmixture_qp.R): once the sum constraint is swept, the last
  // free parameter is determined by the others, its pivot tab[i, i] is 0,
  // and the |tab[i, i]| > 1e-10 guard keeps it free even when it sits on a
  // bound. At a vertex every parameter is on a bound. A fixed parameter
  // whose gradient (relative to that last free one) points into the box is
  // swept in, the step length a is 0 (up to rounding) because the last free
  // parameter cannot move past its bound, and the critical test fixes the
  // parameter again. delta no longer changes, so only time is lost: about
  // 1000 passes, some 200 times the cost of a typical column (K = 3).
  return 0;
}

}  // namespace openadmixture

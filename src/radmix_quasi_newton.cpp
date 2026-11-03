// Quasi-Newton acceleration of the block updates (Zhou, Alexander and Lange
// 2011).
//
// One sweep of block updates (all of Q given P, then all of P given the new
// Q) is a map x -> F(x) of the parameters x = (q, p), and the maximum
// likelihood estimate is one of its fixed points. The fitting loop
// (radmix_fit.cpp) applies F twice per iteration, giving x0 (the current
// estimate), x1 = F(x0) and x2 = F(x1), and keeps the last n_qn secant pairs
// u = x1 - x0 and v = x2 - x1 as the columns of U and V. Near the fixed
// point, F(x) ~ F(x0) + M (x - x0), where M approximates the Jacobian of F;
// the secant conditions M U = V are satisfied by M = V (U'U)^-1 U', the
// smallest such matrix. Solving x = F(x0) + M (x - x0) with this M (by the
// Woodbury identity, without ever forming M) gives the quasi-Newton point
//   x_qn = x1 - V (U'(U - V))^-1 U'(x0 - x1).
// Only an n_qn x n_qn system is solved (3 x 3 by default), so the step costs
// two passes over the parameters, much less than one pass over the
// genotypes. The point ignores the constraints, so it is projected back onto
// the feasible set. Whether it is kept is decided by the fitting loop.
//
// Results do not depend on the number of threads: the long inner products
// are split into chunks of fixed size, whose partial sums are added up in
// chunk order.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include "radmix_projections.h"
#include "radmix_quasi_newton.h"

namespace radmix {

// number of parameters per chunk of the long sums (and of the parallel
// loops over the parameters): large enough to amortise the scheduling, and
// fixed, so that the partial sums do not depend on the number of threads
static const size_t QN_CHUNK = 65536;

// Store a secant pair (see radmix_quasi_newton.h): u = x1 - x0 and
// v = x2 - x1 into column col of U and V, the q block in rows 0 to nq - 1
// and the p block in rows nq to nq + np - 1.
void SecantPairs::store(size_t col, const double* q, const double* p,
                        const double* q1, const double* p1, const double* q2,
                        const double* p2, size_t nq, size_t np) {
  double* u = &U[col * n_par];
  double* v = &V[col * n_par];
  for (size_t t = 0; t < nq; t++) {
    u[t] = q1[t] - q[t];
    v[t] = q2[t] - q1[t];
  }
  for (size_t t = 0; t < np; t++) {
    u[nq + t] = p1[t] - p[t];
    v[nq + t] = p2[t] - p1[t];
  }
}

// The quasi-Newton point (see radmix_quasi_newton.h and the top of this
// file), in three steps: the small system, its solution, and the point
// itself with its projection.
bool qn_step(const double* U, const double* V, size_t n_par, size_t m,
             const double* q, const double* p, const double* q1,
             const double* p1, int K, size_t n_indiv, size_t n_loci,
             double* q_out, double* p_out, int n_cores) {
  // the size of the q block: element t of a secant vector is in the q
  // block if t < nq, and is element t - nq of the p block otherwise
  const size_t nq = (size_t) K * n_indiv;
  // 1. The small system lhs coef = rhs, with lhs = U'(U - V) (m x m) and
  // rhs = U'(x0 - x1), in one pass over the parameters. Each chunk stores
  // its m * m + m partial sums (lhs column-major, then rhs) in its own slot
  // of `partial`.
  const size_t n_chunks = (n_par + QN_CHUNK - 1) / QN_CHUNK;
  const size_t n_sums = m * m + m;
  std::vector<double> partial(n_chunks * n_sums, 0.0);
#pragma omp parallel for num_threads(n_cores) schedule(static)
  for (ptrdiff_t c = 0; c < (ptrdiff_t) n_chunks; c++) {
    const size_t t0 = c * QN_CHUNK, t1 = std::min(n_par, t0 + QN_CHUNK);
    double* sums = &partial[c * n_sums];
    for (size_t t = t0; t < t1; t++) {
      // element t of x0 - x1, from the q or the p block
      const double d = t < nq ? q[t] - q1[t] : p[t - nq] - p1[t - nq];
      for (size_t a = 0; a < m; a++) {
        const double ua = U[a * n_par + t];
        for (size_t b = 0; b < m; b++)
          sums[a + m * b] += ua * (U[b * n_par + t] - V[b * n_par + t]);
        sums[m * m + a] += ua * d;
      }
    }
  }
  // add up the partial sums in chunk order (outside the parallel region)
  arma::mat lhs(m, m, arma::fill::zeros);
  arma::vec rhs(m, arma::fill::zeros);
  for (size_t c = 0; c < n_chunks; c++) {
    const double* sums = &partial[c * n_sums];
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) lhs(a, b) += sums[a + m * b];
      rhs(a) += sums[m * m + a];
    }
  }
  // 2. Solve it (lhs is not symmetric, so by LU). Secant pairs become nearly
  // collinear close to convergence, and a system that is singular to
  // machine precision (by the same test as R's solve()) would turn rounding
  // errors into a meaningless step: it is rejected, as is anything that is
  // not finite. (no_approx stops Armadillo from falling back on a
  // least-squares solution; the negated comparison also rejects a NaN
  // condition number.)
  if (!lhs.is_finite() || !rhs.is_finite()) return false;
  const double rc = arma::rcond(lhs);
  if (!(rc >= DBL_EPSILON)) return false;
  arma::vec coef;
  if (!arma::solve(coef, lhs, rhs, arma::solve_opts::no_approx)) return false;
  if (!coef.is_finite()) return false;
  // 3. x_qn = x1 - V coef, then projected onto the feasible set column by
  // column. Element t of the output depends only on element t of x1 and V,
  // so the output may overwrite x0, which is no longer needed. (The
  // coefficients are copied into a plain vector, read by every thread.)
  std::vector<double> cf(coef.begin(), coef.end());
#pragma omp parallel for num_threads(n_cores) schedule(static)
  for (ptrdiff_t c = 0; c < (ptrdiff_t) n_chunks; c++) {
    const size_t t0 = c * QN_CHUNK, t1 = std::min(n_par, t0 + QN_CHUNK);
    for (size_t t = t0; t < t1; t++) {
      double step = 0.0;
      for (size_t a = 0; a < m; a++) step += V[a * n_par + t] * cf[a];
      if (t < nq) q_out[t] = q1[t] - step;
      else p_out[t - nq] = p1[t - nq] - step;
    }
  }
  // the point ignores the constraints: project each column of q onto the
  // simplex and clamp each column of p, as after a Newton step
#pragma omp parallel for num_threads(n_cores) schedule(static)
  for (ptrdiff_t i = 0; i < (ptrdiff_t) n_indiv; i++)
    project_q(q_out + (size_t) i * K, K);
#pragma omp parallel for num_threads(n_cores) schedule(static)
  for (ptrdiff_t j = 0; j < (ptrdiff_t) n_loci; j++)
    project_p(p_out + (size_t) j * K, K);
  return true;
}

}  // namespace radmix

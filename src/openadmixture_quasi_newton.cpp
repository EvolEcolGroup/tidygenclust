// Quasi-Newton acceleration of Zhou, Alexander and Lange (2011), as in
// update_UV!() and update_QN!() of OpenADMIXTURE.jl (src/quasi_newton.jl).
//
// The fitting loop (src/openadmixture_fit.cpp) applies the block-update map F
// twice per iteration, giving x0 (the current estimate), x1 = F(x0) and x2 =
// F(x1), where x = (q, p) stacks all the parameters. The last n_qn secant pairs
// u = x1 - x0 and v = x2 - x1 are kept as the columns of U and V, and the
// quasi-Newton step jumps towards the fixed point of F:
//   x_qn = x1 - V (U'(U - V))^-1 U'(x0 - x1).
//
// Derivation. We look for the fixed point of F. Near the current estimate,
// F(x) ~ F(x0) + M (x - x0), where M approximates the Jacobian of F. The
// secant pairs require M U = V (each column of V is the change in F caused by
// the change in x in the corresponding column of U), and the smallest such
// matrix is M = V (U'U)^-1 U'. Solving x = F(x0) + M (x - x0) with this M
// (using the Woodbury identity to avoid ever forming the large matrix M)
// gives the formula above, with x1 = F(x0). Only an n_qn x n_qn system needs
// to be solved, so the step is cheap even with millions of parameters.
//
// Why this is in C++
// ------------------
// The step only involves vectors of length K (I + J) and an n_qn x n_qn
// linear system (3 x 3 by default), so it is cheap; it was first written in
// R. But in R every operation allocates: stacking q and p into x0, x1 and x2,
// the differences x1 - x0 and x2 - x1, U - V inside crossprod(), and so on.
// With many loci, each of these is a copy of P, and since R frees memory
// only when the garbage collector runs, the peak memory of a fit reached
// about 30 copies of the parameters. Here U and V are kept in a single
// object for the whole fit, and the step is computed directly from the Q and
// P blocks, without stacking them and without temporary vectors, writing the
// result into buffers owned by the fitting loop. This is purely an
// implementation change: the arithmetic is the same as in the Julia code (up
// to rounding, as sums are accumulated in a different order).
//
// Determinism: the long sums (inner products over K (I + J) values) are
// split into fixed chunks whose partial sums are added up in order, so the
// result does not depend on the number of threads.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include "openadmixture_projections.h"
#include "openadmixture_quasi_newton.h"

using namespace Rcpp;

// internal code of the OpenADMIXTURE implementation (see the note on the
// namespace in openadmixture_geno_access.h)
namespace openadmixture {


// number of values per chunk of the long sums (fixed, so that results do not
// depend on the number of threads)
static const size_t QN_CHUNK = 65536;

// Store the secant pair of one iteration (see openadmixture_quasi_newton.h).
// Column col of U receives u = x1 - x0 (how far the first sweep moved the
// parameters) and column col of V receives v = x2 - x1 (how far the second
// sweep moved them). The fitting loop cycles col through 0, ..., n_qn - 1, so
// the oldest pair is overwritten and U and V always hold the last n_qn pairs.
void qn_push(QNState& st, size_t col, const double* q, const double* p,
             const double* q1, const double* p1, const double* q2,
             const double* p2, size_t nq, size_t np) {
  // start of column col of U and V (columns are n_par values long)
  double* u = &st.U[col * st.n_par];
  double* v = &st.V[col * st.n_par];
  // x = (q, p) is never stacked in memory: the first nq positions of a
  // column hold the q block and the next np positions the p block, so we
  // fill them from the two blocks in turn
  for (size_t t = 0; t < nq; t++) {
    u[t] = q1[t] - q[t];
    v[t] = q2[t] - q1[t];
  }
  for (size_t t = 0; t < np; t++) {
    u[nq + t] = p1[t] - p[t];
    v[nq + t] = p2[t] - p1[t];
  }
}

// Compute the quasi-Newton point x_qn = x1 - V (U'(U - V))^-1 U'(x0 - x1)
// from the first n_used secant pairs, and project it onto the feasible set
// (see openadmixture_quasi_newton.h). The work has three stages:
//   1. one pass over the parameters to form the small m x m system
//      (m = n_used, usually 3);
//   2. solving that system for the m coefficients;
//   3. one pass over the parameters to form x_qn and project it.
// Stages 1 and 3 touch every parameter and run in parallel; stage 2 is tiny.
bool qn_step(const QNState& st, size_t n_used, const double* q,
             const double* p, const double* q1, const double* p1, int K,
             size_t n_indiv, size_t n_loci, double* q_out, double* p_out,
             int n_cores) {
  // nq: size of the q block (positions 0 to nq - 1 of x); the p block
  // follows. m: number of secant pairs used.
  const size_t nq = (size_t) K * n_indiv, n_par = st.n_par, m = n_used;
  const double *U = st.U.data(), *V = st.V.data();

  // The small system: lhs = U'(U - V) (m x m) and rhs = U'(x0 - x1) (m),
  // i.e. for each pair of columns (a, b), sum_t U[t, a] (U[t, b] - V[t, b]),
  // and for each column a, sum_t U[t, a] (x0[t] - x1[t]). These are computed
  // in one pass over the parameters, chunk by chunk, with the partial sums
  // of each chunk stored separately and added up in chunk order afterwards.
  // Each chunk needs m * m sums for lhs and m sums for rhs; they are stored
  // next to each other (lhs column-major, then rhs) in `partial`, one group
  // of n_sums values per chunk, so that threads never write to the same
  // place.
  const size_t n_chunks = (n_par + QN_CHUNK - 1) / QN_CHUNK;
  const size_t n_sums = m * m + m;
  std::vector<double> partial(n_chunks * n_sums, 0.0);
#pragma omp parallel for num_threads(n_cores) schedule(static)
  for (ptrdiff_t c = 0; c < (ptrdiff_t) n_chunks; c++) {
    const size_t t0 = c * QN_CHUNK, t1 = std::min(n_par, t0 + QN_CHUNK);
    double* sums = &partial[c * n_sums];
    for (size_t t = t0; t < t1; t++) {
      // x0 - x1 at position t, read from the q or the p block
      const double d = t < nq ? q[t] - q1[t] : p[t - nq] - p1[t - nq];
      // U[a * n_par + t] is element t of column a of U (column-major)
      for (size_t a = 0; a < m; a++) {
        const double ua = U[a * n_par + t];
        // lhs(a, b) += U[t, a] (U[t, b] - V[t, b])
        for (size_t b = 0; b < m; b++)
          sums[a + m * b] += ua * (U[b * n_par + t] - V[b * n_par + t]);
        // rhs(a) += U[t, a] (x0[t] - x1[t])
        sums[m * m + a] += ua * d;
      }
    }
  }
  // add up the chunks, always in chunk order (deterministic for any number
  // of threads), into the small system held as Armadillo objects
  arma::mat lhs(m, m, arma::fill::zeros);
  arma::vec rhs(m, arma::fill::zeros);
  for (size_t c = 0; c < n_chunks; c++) {
    const double* sums = &partial[c * n_sums];
    for (size_t a = 0; a < m; a++) {
      for (size_t b = 0; b < m; b++) lhs(a, b) += sums[a + m * b];
      rhs(a) += sums[m * m + a];
    }
  }

  // Solve the small system. As with solve() in R (which the first version of
  // this code used), a system that is singular to machine precision is
  // rejected: this happens when secant pairs are (nearly) collinear, e.g.
  // close to convergence, and the step is then skipped.
  // rcond() estimates the reciprocal condition number of lhs: values close
  // to 0 mean that the system is nearly singular, and solving it would
  // amplify rounding errors into a meaningless step. R's solve() uses the
  // same test, with the same threshold (machine epsilon). The negated
  // comparison also rejects NaN. no_approx stops Armadillo from falling back
  // on an approximate (least-squares) solution.
  if (!lhs.is_finite() || !rhs.is_finite()) return false;
  const double rc = arma::rcond(lhs);
  if (!(rc >= DBL_EPSILON)) return false;
  arma::vec coef;
  if (!arma::solve(coef, lhs, rhs, arma::solve_opts::no_approx)) return false;
  if (!coef.is_finite()) return false;

  // x_qn = x1 - V coef, written straight into the output buffers, then
  // projected column by column (columns of K values)
  // (the coefficients are copied into a plain vector so that the parallel
  // loop below only reads raw memory)
  std::vector<double> cf(coef.begin(), coef.end());
#pragma omp parallel for num_threads(n_cores) schedule(static)
  for (ptrdiff_t c = 0; c < (ptrdiff_t) n_chunks; c++) {
    const size_t t0 = c * QN_CHUNK, t1 = std::min(n_par, t0 + QN_CHUNK);
    for (size_t t = t0; t < t1; t++) {
      // element t of V coef
      double step = 0.0;
      for (size_t a = 0; a < m; a++) step += V[a * n_par + t] * cf[a];
      // position t of x belongs to the q block or to the p block
      if (t < nq) q_out[t] = q1[t] - step;
      else p_out[t - nq] = p1[t - nq] - step;
    }
  }
  // The extrapolated point ignores the constraints: bring each individual's
  // ancestry proportions back onto the simplex, and each allele frequency
  // into [1e-5, 1 - 1e-5] (openadmixture_projections.cpp). Columns are
  // independent, so these loops are parallel too.
#pragma omp parallel for num_threads(n_cores) schedule(static)
  for (ptrdiff_t i = 0; i < (ptrdiff_t) n_indiv; i++)
    project_q(q_out + (size_t) i * K, K);
#pragma omp parallel for num_threads(n_cores) schedule(static)
  for (ptrdiff_t j = 0; j < (ptrdiff_t) n_loci; j++)
    project_p(p_out + (size_t) j * K, K);
  return true;
}

}  // namespace openadmixture

using namespace openadmixture;

// ---------------------------------------------------------------------------
// R interface, used by the tests (the fitting loop calls the functions above
// directly)

static QNState& get_state(SEXP state) {
  XPtr<QNState> ptr(state);
  if (ptr.get() == NULL)
    stop("the quasi-Newton state is no longer valid");
  return *ptr;
}

// Create the secant pairs for n_par = K (I + J) parameters, keeping the last
// n_qn pairs.
// [[Rcpp::export]]
SEXP openadmixture_qn_new_cpp(double n_par, int n_qn) {
  if (n_par < 1 || n_qn < 1) stop("'n_par' and 'n_qn' must be positive");
  XPtr<QNState> ptr(new QNState((size_t) n_par, (size_t) n_qn), true);
  return ptr;
}

// Store a secant pair in column `col` (1-based); see qn_push().
// [[Rcpp::export]]
void openadmixture_qn_push_cpp(SEXP state, int col, const NumericMatrix& q,
                               const NumericMatrix& p, const NumericMatrix& q1,
                               const NumericMatrix& p1, const NumericMatrix& q2,
                               const NumericMatrix& p2) {
  QNState& st = get_state(state);
  const size_t nq = q.size(), np = p.size();
  if (nq + np != st.n_par || (size_t) q1.size() != nq ||
      (size_t) q2.size() != nq || (size_t) p1.size() != np ||
      (size_t) p2.size() != np)
    stop("the dimensions of q and p do not match the quasi-Newton state");
  if (col < 1 || (size_t) col > st.n_qn) stop("'col' is out of range");
  qn_push(st, col - 1, q.begin(), p.begin(), q1.begin(), p1.begin(),
          q2.begin(), p2.begin(), nq, np);
}

// The projected quasi-Newton step from (q1, p1) using the first n_used
// pairs; returns list(q, p), or NULL if the system is singular (see
// qn_step()).
// [[Rcpp::export]]
SEXP openadmixture_qn_step_cpp(SEXP state, int n_used, const NumericMatrix& q,
                               const NumericMatrix& p, const NumericMatrix& q1,
                               const NumericMatrix& p1, int n_cores) {
  QNState& st = get_state(state);
  const size_t nq = q.size(), np = p.size();
  if (nq + np != st.n_par || (size_t) q1.size() != nq ||
      (size_t) p1.size() != np)
    stop("the dimensions of q and p do not match the quasi-Newton state");
  if (n_used < 1 || (size_t) n_used > st.n_qn) stop("'n_used' is out of range");
  if (n_cores < 1) stop("'n_cores' must be at least 1");
  const int K = q.nrow();
  NumericMatrix q_new(K, q.ncol()), p_new(K, p.ncol());
  if (!qn_step(st, n_used, q.begin(), p.begin(), q1.begin(), p1.begin(), K,
               q.ncol(), p.ncol(), q_new.begin(), p_new.begin(), n_cores))
    return R_NilValue;
  return List::create(Named("q") = q_new, Named("p") = p_new);
}

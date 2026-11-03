// Starting values for gt_radmix() from the SVD of the genotypes
// (init = "svd"): alternating least squares (ALS) on the factors of a
// truncated SVD, following ADAMIXTURE (model/als.py) and, more loosely,
// fastmixture (factorALS() in fastmixture/utils.py).
//
// The factorisation
// -----------------
// Under the ADMIXTURE model, the expected frequency of the alternate allele
// in individual i at locus j is pi_ji = sum_k p_jk q_ik, i.e. Pi = P Q'
// (loci x individuals). The truncated SVD of the centred and scaled
// genotypes, computed in R (radmix_svd_factors() in R/radmix_svd.R), gives a
// low-rank approximation Pi ~ Za Ua', where Za (loci x R) and Ua
// (individuals x R) are the scaled singular vectors, plus a last column
// holding the allele frequencies (in Za) and ones (in Ua). ALS looks for P
// and Q such that P Q' ~ Za Ua', solving alternately for Q given P and for
// P given Q. Each step only involves products with the factors, never the
// genotypes, so an iteration costs O((I + J) K R), against O(I J K) for a
// pass over the genotypes.
//
// As in ADAMIXTURE:
// - each individual (and each locus) is solved by least squares bounded in
//   [0, 1] (ADAMIXTURE uses scipy's BVLS; we use solve_rows(), below); Q
//   is then clipped to [1e-5, 1 - 1e-5] and its columns normalised to sum
//   to one, and P is clipped to [1e-5, 1 - 1e-5];
// - the Q step weights each locus by the inverse of its binomial variance,
//   1 / (f (1 - f)) (bounded and normalised to mean one in R), so that the
//   loci are weighted roughly as in the likelihood;
// - the P step uses Q extrapolated along its last change (momentum 0.8);
// - iterations stop when the root mean square change in Q falls below tol.
//
// Layout and linear algebra
// -------------------------
// As in the kernels, q is K x I and p is K x J, column-major, so the K
// values of an individual (or a locus), which are solved for together, are
// contiguous. The products with the factors are tall and skinny matrix
// products, computed with Armadillo, hence with R's BLAS (which
// radmix_svd_start() limits to one thread: BLAS threads gain nothing on
// these products, and compete with the OpenMP threads); the small bounded
// least squares problems, one per individual or locus, are solved in
// parallel with OpenMP. Each
// column is solved on its own, so the results do not depend on the number
// of threads. Most of the time goes into these small problems (about 90%
// on EUR_subset), hence the direct solution in solve_rows().

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include "radmix_projections.h"

using namespace Rcpp;

// internal code of gt_radmix() (see the note on the namespace in
// radmix_geno_access.h)
namespace radmix {

// coordinate descent: maximum number of sweeps over the K coordinates, and
// tolerance on the largest change in a sweep
static const int CD_MAX_SWEEPS = 20;
static const double CD_TOL = 1e-7;
// the weight of the last change of Q in its extrapolation (ADAMIXTURE's
// momentum)
static const double ALS_MOMENTUM = 0.8;

// the bounds of the estimates, as in the rest of gt_radmix()
static inline double als_clip(double v) {
  return std::min(std::max(v, RADMIX_EPS), 1.0 - RADMIX_EPS);
}

// The x in [0, 1]^K that minimises 0.5 x' A x - b' x, with A (K x K,
// column-major) symmetric and positive semi-definite, by cyclic coordinate
// descent starting from the value in x: each coordinate in turn is set to
// its exact minimiser given the others, clamped to [0, 1]. A coordinate
// with A_kk = 0 (a flat direction) keeps its value.
static void box_cd(const double* A, const double* b, double* x, int K) {
  for (int s = 0; s < CD_MAX_SWEEPS; s++) {
    double change = 0.0;
    for (int k = 0; k < K; k++) {
      // column k of A, which is also row k
      const double* a = A + (size_t) K * k;
      if (!(a[k] > 0.0)) continue;
      // the gradient of the objective along x_k is (A x)_k - b_k
      double r = b[k];
      for (int l = 0; l < K; l++) r -= a[l] * x[l];
      const double xk = std::min(std::max(x[k] + r / a[k], 0.0), 1.0);
      change = std::max(change, std::fabs(xk - x[k]));
      x[k] = xk;
    }
    if (change < CD_TOL) break;
  }
}

// Solve the bounded least squares problems of all the columns of xt (K x n),
// minimising 0.5 x' A x - b' x in [0, 1]^K for each column x of xt and the
// matching column b of bt. All the problems share A, so their unconstrained
// solutions A^-1 b are obtained at once, with a single matrix product; a
// solution within [0, 1]^K is also the solution of the bounded problem (the
// objective is convex), and is used as such (this saves about half of the
// time of the ALS on EUR_subset, where 80% of the loci have interior
// solutions near convergence). The other columns are solved by coordinate
// descent, from their current value. If A is ill-conditioned, its inverse
// would be inaccurate, so all the columns use coordinate descent.
static void solve_rows(const arma::mat& A, const arma::mat& bt,
                       arma::mat& xt, int n_cores) {
  const int K = xt.n_rows;
  arma::mat a_inv;
  const bool direct = arma::rcond(A) > 1e-8 && arma::inv(a_inv, A);
  arma::mat x_direct;
  if (direct) x_direct = a_inv * bt;
#pragma omp parallel for num_threads(n_cores) schedule(static)
  for (arma::uword i = 0; i < xt.n_cols; i++) {
    double* x = xt.colptr(i);
    bool inside = direct;
    if (direct) {
      const double* xd = x_direct.colptr(i);
      for (int k = 0; k < K; k++) inside = inside && xd[k] >= 0.0 &&
        xd[k] <= 1.0;
      if (inside) std::copy(xd, xd + K, x);
    }
    if (!inside) box_cd(A.memptr(), bt.colptr(i), x, K);
  }
}

// Clip the K proportions of one individual and normalise them to sum to one
// (ADAMIXTURE's map_Q)
static void map_q_col(double* q, int K) {
  double sum = 0.0;
  for (int k = 0; k < K; k++) {
    q[k] = als_clip(q[k]);
    sum += q[k];
  }
  for (int k = 0; k < K; k++) q[k] /= sum;
}

// Q step: for each individual i, the q_i in [0, 1]^K that minimises
//   sum_j w_j (pi_ji - sum_k p_jk q_ik)^2,
// i.e. 0.5 q' A q - b_i' q with A = P' W P and b_i = P' W Za ua_i, starting
// from (and overwriting) qt; then mapped onto the simplex.
static void q_step(const arma::mat& pt, const arma::mat& za,
                   const arma::mat& ua, const arma::vec& w, arma::mat& qt,
                   int n_cores) {
  const int K = qt.n_rows;
  // P' W, K x J
  const arma::mat wpt = pt.each_row() % w.t();
  const arma::mat A = wpt * pt.t();
  // b_i for all the individuals, K x I
  const arma::mat bt = (wpt * za) * ua.t();
  solve_rows(A, bt, qt, n_cores);
  for (arma::uword i = 0; i < qt.n_cols; i++) map_q_col(qt.colptr(i), K);
}

// P step: for each locus j, the p_j in [0, 1]^K that minimises
//   sum_i (pi_ji - sum_k p_jk q_ik)^2,
// i.e. 0.5 p' A p - b_j' p with A = Q' Q and b_j = Q' Ua za_j, starting
// from (and overwriting) pt; then clipped.
static void p_step(const arma::mat& qt, const arma::mat& za,
                   const arma::mat& ua, arma::mat& pt, int n_cores) {
  const arma::mat A = qt * qt.t();
  // b_j for all the loci, K x J
  const arma::mat bt = (qt * ua) * za.t();
  solve_rows(A, bt, pt, n_cores);
  pt.transform(als_clip);
}

}  // namespace radmix

using namespace radmix;

// ALS on the factors of a truncated SVD.
//
// za: loci x R, ua: individuals x R, so that Pi ~ za ua' (see above).
// w: the weight of each locus in the Q step.
// init_p: the starting allele frequencies (loci x K, the user layout),
//   clipped to [1e-5, 1 - 1e-5]; Q starts from its Q step, from 1 / K.
// max_iter, tol: maximum number of iterations, and threshold on the root
//   mean square change in Q between iterations.
// Returns list(q = K x I, p = K x J, n_iter, converged), in the layout of
// the kernels (e.g. radmix_loglik_cpp()).
// [[Rcpp::export]]
List radmix_als_cpp(NumericMatrix za, NumericMatrix ua, NumericVector w,
                    NumericMatrix init_p, int max_iter, double tol,
                    int n_cores) {
  // the arguments are prepared in R (radmix_svd_start()); these checks only
  // guard against misuse from within the package
  const arma::uword J = za.nrow(), I = ua.nrow(), R = za.ncol(),
    K = init_p.ncol();
  if (ua.ncol() != (int) R || w.size() != (R_xlen_t) J ||
      init_p.nrow() != (int) J || K < 1 || I == 0 || J == 0 ||
      max_iter < 1 || n_cores < 1)
    stop("invalid arguments to radmix_als_cpp()");
  // views of the R matrices (no copies)
  const arma::mat Za(za.begin(), J, R, false, true);
  const arma::mat Ua(ua.begin(), I, R, false, true);
  const arma::vec W(w.begin(), J, false, true);
  // the estimates, in the layout of the kernels
  arma::mat pt = arma::mat(init_p.begin(), J, K, false, true).t();
  pt.transform(als_clip);
  arma::mat qt(K, I);
  qt.fill(1.0 / K);
  q_step(pt, Za, Ua, W, qt, n_cores);
  // the previous Q (for the extrapolation), the extrapolated Q used in the
  // P step, and the new Q
  arma::mat q_prev = qt, q_use(K, I), q_new(K, I);
  bool converged = false;
  int iter = 0;
  for (iter = 1; iter <= max_iter; iter++) {
    if (iter == 1) {
      q_use = qt;
    } else {
      q_use = qt + ALS_MOMENTUM * (qt - q_prev);
      for (arma::uword i = 0; i < I; i++) map_q_col(q_use.colptr(i), K);
    }
    p_step(q_use, Za, Ua, pt, n_cores);
    // the Q step starts from the current Q (not the extrapolated one)
    q_new = qt;
    q_step(pt, Za, Ua, W, q_new, n_cores);
    const double rmse =
      std::sqrt(arma::accu(arma::square(q_new - qt)) / (double) (K * I));
    // q_prev <- qt, qt <- q_new (the old q_prev becomes a spare buffer)
    q_prev.swap(qt);
    qt.swap(q_new);
    // (a NaN rmse is never below tol, so such a fit runs to max_iter)
    if (rmse < tol) {
      converged = true;
      break;
    }
    checkUserInterrupt();
  }
  // without convergence, the loop has left iter at max_iter + 1
  if (iter > max_iter) iter = max_iter;
  // map_q_col() keeps each column of Q on the simplex, but its
  // normalisation can take entries slightly below 1e-5 (when clipping has
  // raised the sum above one), so Q is projected onto the feasible set of
  // the fit, as are starting values given by the user (radmix_check_init())
  for (arma::uword i = 0; i < I; i++) project_q(qt.colptr(i), K);
  return List::create(
    Named("q") = wrap(qt), Named("p") = wrap(pt), Named("n_iter") = iter,
    Named("converged") = converged);
}

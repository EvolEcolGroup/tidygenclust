// R interface of the building blocks of gt_radmix(): the passes over the
// genotypes, the QP solver, the projections, the quasi-Newton step and the
// decoding of the genotypes. Besides radmix_fit_cpp() (radmix_fit.cpp) and
// radmix_als_cpp() (radmix_als.cpp), the R code of gt_radmix() only calls
// three functions of this file: radmix_loglik_cpp() (to rank the SVD
// starts), radmix_deviance_cpp() (cross-validation) and
// radmix_project_cpp() (starting values of Q). The others exist so that
// each building block can be tested on its own.
//
// The C++ functions trust their inputs (raw pointer arithmetic, no R calls
// inside parallel regions), so every wrapper checks the dimensions of its
// arguments, with Rcpp::stop(), before calling them, and allocates the R
// outputs beforehand. q and p are passed as numeric matrices, K x I and
// K x J (the internal layout); Rcpp passes a double matrix without copying
// it, and the C++ functions never write into their inputs.

#include <Rcpp.h>
#include "radmix_geno_access.h"
#include "radmix_kernels.h"
#include "radmix_projections.h"
#include "radmix_qp_sweep.h"
#include "radmix_quasi_newton.h"

using namespace Rcpp;
using namespace radmix;

// Check q (K x I) and p (K x J) against the genotypes (I individuals and J
// loci, both at least 1), and n_cores >= 1, with an R error otherwise;
// returns K. Called on the main thread, before any kernel.
static int check_dims(const GenoView& geno, const NumericMatrix& q,
                      const NumericMatrix& p, int n_cores) {
  const int K = q.nrow();
  if (K < 1) stop("'q' must have at least one row.");
  if (p.nrow() != K) stop("'q' and 'p' must have the same number of rows (K).");
  if ((size_t) q.ncol() != geno.nrow())
    stop("'q' must have one column per individual (%d), not %d.",
         (int) geno.nrow(), q.ncol());
  if ((size_t) p.ncol() != geno.ncol())
    stop("'p' must have one column per locus (%d), not %d.",
         (int) geno.ncol(), p.ncol());
  if (geno.nrow() == 0 || geno.ncol() == 0)
    stop("There must be at least one individual and one locus.");
  if (n_cores < 1) stop("'n_cores' must be at least 1.");
  return K;
}

// The mask of fixed individuals of em() and update_q() (radmix_kernels.h,
// "Fixed individuals") from an optional logical vector with one value per
// individual (TRUE = fixed): fills `mask` and returns its data, or NULL
// when fixed is NULL. Checked on the main thread, before any kernel.
static const unsigned char* fixed_mask(const GenoView& geno,
                                       const Nullable<LogicalVector>& fixed,
                                       std::vector<unsigned char>& mask) {
  if (fixed.isNull()) return nullptr;
  LogicalVector fv(fixed.get());
  if ((size_t) fv.size() != geno.nrow())
    stop("'fixed' must have one value per individual (%d), not %d.",
         (int) geno.nrow(), (int) fv.size());
  mask.assign(geno.nrow(), 0);
  for (size_t i = 0; i < mask.size(); i++) {
    if (fv[i] == NA_LOGICAL) stop("'fixed' must not contain NA.");
    mask[i] = fv[i] ? 1 : 0;
  }
  return mask.data();
}

// Log-likelihood (a scalar) of q (K x I) and p (K x J) for the genotypes of
// data (the output of radmix_prepare(), with an optional cross-validation
// mask); the same arguments are used by the other passes below.
// [[Rcpp::export]]
double radmix_loglik_cpp(List data, NumericMatrix q, NumericMatrix p,
                         int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  return loglik(geno, q.begin(), p.begin(), K, n_cores);
}

// Binomial deviance of the visible genotypes, for cross-validation (data
// then carries a held-out mask); returns c(deviance, number of genotypes).
// [[Rcpp::export]]
NumericVector radmix_deviance_cpp(List data, NumericMatrix q,
                                  NumericMatrix p, int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  double dev, count;
  deviance(geno, q.begin(), p.begin(), K, n_cores, &dev, &count);
  return NumericVector::create(dev, count);
}

// One EM step; returns list(q = K x I, p = K x J). The individuals marked
// TRUE in the optional `fixed` keep their column of q.
// [[Rcpp::export]]
List radmix_em_cpp(List data, NumericMatrix q, NumericMatrix p,
                   int n_cores,
                   Nullable<LogicalVector> fixed = R_NilValue) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  std::vector<unsigned char> mask;
  const unsigned char* fx = fixed_mask(geno, fixed, mask);
  NumericMatrix q_new(K, q.ncol()), p_new(K, p.ncol());
  em(geno, q.begin(), p.begin(), K, q_new.begin(), p_new.begin(), n_cores,
     fx);
  return List::create(_["q"] = q_new, _["p"] = p_new);
}

// Gradient and negative Hessian of the log-likelihood with respect to each
// column of q; returns list(grad = K x I, hess = K x K x I array).
// [[Rcpp::export]]
List radmix_grad_q_cpp(List data, NumericMatrix q, NumericMatrix p,
                       int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  const int I = q.ncol();
  NumericMatrix grad(K, I);
  // (the size in size_t, to avoid an int overflow for large K * K * I)
  NumericVector hess((size_t) K * K * I);
  hess.attr("dim") = IntegerVector::create(K, K, I);
  grad_q(geno, q.begin(), p.begin(), K, grad.begin(), hess.begin(), n_cores);
  return List::create(_["grad"] = grad, _["hess"] = hess);
}

// The same for the columns of p; returns list(grad = K x J,
// hess = K x K x J array).
// [[Rcpp::export]]
List radmix_grad_p_cpp(List data, NumericMatrix q, NumericMatrix p,
                       int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  const int J = p.ncol();
  NumericMatrix grad(K, J);
  NumericVector hess((size_t) K * K * J);
  hess.attr("dim") = IntegerVector::create(K, K, J);
  grad_p(geno, q.begin(), p.begin(), K, grad.begin(), hess.begin(), n_cores);
  return List::create(_["grad"] = grad, _["hess"] = hess);
}

// Block update of q (a constrained Newton step for each column); returns the
// new K x I matrix. The individuals marked TRUE in the optional `fixed`
// keep their column of q.
// [[Rcpp::export]]
NumericMatrix radmix_update_q_cpp(List data, NumericMatrix q,
                                  NumericMatrix p, int n_cores,
                                  Nullable<LogicalVector> fixed =
                                    R_NilValue) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  std::vector<unsigned char> mask;
  const unsigned char* fx = fixed_mask(geno, fixed, mask);
  NumericMatrix q_new(K, q.ncol());
  update_q(geno, q.begin(), p.begin(), K, q_new.begin(), n_cores, fx);
  return q_new;
}

// Block update of q and log-likelihood of (q, p) in one pass; returns
// list(q = the new K x I matrix, loglik).
// [[Rcpp::export]]
List radmix_update_q_loglik_cpp(List data, NumericMatrix q,
                                NumericMatrix p, int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  NumericMatrix q_new(K, q.ncol());
  const double ll = update_q_loglik(geno, q.begin(), p.begin(), K,
                                    q_new.begin(), n_cores);
  return List::create(_["q"] = q_new, _["loglik"] = ll);
}

// Block update of p; returns the new K x J matrix.
// [[Rcpp::export]]
NumericMatrix radmix_update_p_cpp(List data, NumericMatrix q,
                                  NumericMatrix p, int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  NumericMatrix p_new(K, p.ncol());
  update_p(geno, q.begin(), p.begin(), K, p_new.begin(), n_cores);
  return p_new;
}

// Number of segments of loci of the Q-side passes for I individuals and J
// loci (1 = none), so that the tests can exercise both code paths. (I and J
// are R numbers, converted to size_t, so they must be non-negative.)
// [[Rcpp::export]]
int radmix_q_segments_cpp(double I, double J) {
  return (int) q_segments((size_t) I, (size_t) J);
}

// Solve the QP of one column (see radmix_qp_sweep.h): H is the K x K
// negative Hessian, grad and par the gradient and the current value of the
// column, simplex TRUE for a column of q. Returns list(delta = the step,
// iter = the return value of qp_solve_column(): the number of iterations of
// the active-set loop, or 0 on failure, on cycling at a vertex of the
// simplex or at the iteration limit, passes = the number of passes of the
// loop, also when it stopped early).
// [[Rcpp::export]]
List radmix_qp_solve_cpp(const NumericMatrix& H, const NumericVector& grad,
                         const NumericVector& par, bool simplex) {
  const int K = grad.size();
  if (H.nrow() != K || H.ncol() != K || par.size() != K)
    stop("H must be K x K and grad, par of length K");
  QPWork work(K);
  NumericVector delta(K);
  int iter = qp_solve_column(H.begin(), grad.begin(), par.begin(), K,
                             simplex, delta.begin(), work);
  return List::create(Named("delta") = delta, Named("iter") = iter,
                      Named("passes") = work.passes);
}

// Project each column of x onto the feasible set: the simplex (sum = 1, all
// values at least 1e-5) if simplex is TRUE, as for q, or the box
// [1e-5, 1 - 1e-5] otherwise, as for p. Returns a new matrix (x is not
// modified). Used by radmix_check_init() for the starting values of Q, and
// by the tests.
// [[Rcpp::export]]
NumericMatrix radmix_project_cpp(const NumericMatrix& x, bool simplex) {
  if (x.nrow() < 1) stop("'x' must have at least one row.");
  NumericMatrix out = clone(x);
  const int K = out.nrow(), n = out.ncol();
  for (int i = 0; i < n; i++) {
    double* col = out.begin() + (size_t) i * K;
    if (simplex) project_q(col, K);
    else project_p(col, K);
  }
  return out;
}

// The quasi-Newton point from the secant pairs in the columns of u and v
// (each n_par x m, the q block first, then the p block), the current
// estimate (q, p) and the estimate after one sweep (q1, p1), all in the
// internal layout (K x I and K x J). Returns list(q, p), the projected
// point in the same layout, or NULL if the step fails (see qn_step()).
// [[Rcpp::export]]
SEXP radmix_qn_step_cpp(const NumericMatrix& u, const NumericMatrix& v,
                        const NumericMatrix& q, const NumericMatrix& p,
                        const NumericMatrix& q1, const NumericMatrix& p1,
                        int n_cores) {
  const int K = q.nrow();
  const size_t nq = q.size(), np = p.size(), m = u.ncol();
  if (K < 1 || p.nrow() != K || q1.nrow() != K || p1.nrow() != K)
    stop("q, p, q1 and p1 must have the same number of rows (K)");
  if ((size_t) q1.size() != nq || (size_t) p1.size() != np)
    stop("q1 and p1 must have the dimensions of q and p");
  if ((size_t) u.nrow() != nq + np || v.nrow() != u.nrow() ||
      v.ncol() != u.ncol() || m < 1)
    stop("u and v must have one row per parameter and the same columns");
  if (n_cores < 1) stop("'n_cores' must be at least 1");
  NumericMatrix q_new(K, q.ncol()), p_new(K, p.ncol());
  if (!qn_step(u.begin(), v.begin(), nq + np, m, q.begin(), p.begin(),
               q1.begin(), p1.begin(), K, q.ncol(), p.ncol(), q_new.begin(),
               p_new.begin(), n_cores))
    return R_NilValue;
  return List::create(Named("q") = q_new, Named("p") = p_new);
}

// The genotypes of the view as a dense matrix (individuals x loci, NA for
// missing or hidden), so that the tests can check the decoding and the
// cross-validation masks (through radmix_decode() in R). The fit never
// builds this matrix.
// [[Rcpp::export]]
NumericMatrix radmix_decode_cpp(const List& data) {
  const GenoView geno(data);
  const size_t n = geno.nrow(), m = geno.ncol();
  NumericMatrix out(n, m);
  double g;
  for (size_t j = 0; j < m; j++)
    for (size_t i = 0; i < n; i++)
      out(i, j) = geno.get(i, j, g) ? g : NA_REAL;
  return out;
}

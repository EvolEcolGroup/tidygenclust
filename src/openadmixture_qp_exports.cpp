// Rcpp wrappers around the QP solver and the projections. They are not on
// the hot path: the Newton updates call qp_solve_column, project_q and
// project_p directly from C++ (src/openadmixture_kernels.cpp), and so does the
// quasi-Newton step (src/openadmixture_quasi_newton.cpp). R/gt_openadmixture.R
// uses openadmixture_project_q_cpp on user-supplied starting values; the other
// wrappers exist so that the solver and the projections can be checked from
// testthat.
#include <Rcpp.h>
#include "openadmixture_qp_sweep.h"
#include "openadmixture_projections.h"
using namespace Rcpp;
using namespace openadmixture;

// Solve one column QP (see openadmixture_qp_sweep.h for the problem statement).
// H: K x K negative Hessian; grad, par: length-K gradient and current
// column; simplex: TRUE for a Q column (sum constraint), FALSE for P.
// Returns list(delta = step, iter = active-set iterations, 0 on failure,
// passes = number of passes of the active-set loop). iter is also 0 when
// the active set cycles at a vertex of the simplex, although delta is then
// correct (see the end of qp_solve_column), so tests should judge the
// solver by delta rather than by iter.
// [[Rcpp::export]]
List openadmixture_qp_solve_cpp(const NumericMatrix& H,
                                const NumericVector& grad,
                                const NumericVector& par, bool simplex) {
  const int K = grad.size();
  // the C++ solver does no checking of its own, so validate shapes here
  if (H.nrow() != K || H.ncol() != K || par.size() != K)
    stop("H must be K x K and grad, par of length K");
  // a fresh scratch object per call; fine for testing, the kernels reuse
  // one per thread instead
  QPWork work(K);
  NumericVector delta(K);
  // R matrices are column-major, as the solver expects, so the data
  // pointers can be passed straight through
  int iter = qp_solve_column(H.begin(), grad.begin(), par.begin(), K,
                             simplex, delta.begin(), work);
  return List::create(Named("delta") = delta, Named("iter") = iter,
                      Named("passes") = work.passes);
}

// Project each column of q onto {sum == 1, entries >= 1e-5}; returns a copy.
// q is K x I (one column per individual, the layout used internally).
// [[Rcpp::export]]
NumericMatrix openadmixture_project_q_cpp(const NumericMatrix& q) {
  // project_q() needs at least one component (it reads the last one)
  if (q.nrow() < 1) stop("'q' must have at least one row.");
  // clone so that the caller's matrix is not modified in place
  NumericMatrix out = clone(q);
  const int K = out.nrow(), n = out.ncol();
  double* x = out.begin();
  // columns are contiguous in column-major storage; the size_t cast keeps
  // the offset from overflowing int for large I
  for (int i = 0; i < n; i++) project_q(x + (size_t) i * K, K);
  return out;
}

// Clamp every entry of p to [1e-5, 1 - 1e-5]; returns a copy.
// p is K x J (one column per locus, the layout used internally).
//
// Clamping is entry-wise, so we could clamp the whole matrix in one call,
// and the first version did exactly that, passing out.size() as the length.
// But project_p() takes an int (it was written for a single column of K
// values, as in Julia), and out.size() is an R_xlen_t: with more than
// 2^31 - 1 entries (e.g. K = 10 and over 215 million loci) the length would
// overflow and the call would clamp the wrong number of entries. Clamping
// one column at a time, as project_q() does above, keeps every length equal
// to K, and the size_t cast keeps the offset from overflowing. (This
// wrapper is specific to this package, so there is no Julia counterpart.)
// [[Rcpp::export]]
NumericMatrix openadmixture_project_p_cpp(const NumericMatrix& p) {
  NumericMatrix out = clone(p);
  const int K = out.nrow(), n = out.ncol();
  double* x = out.begin();
  for (int j = 0; j < n; j++) project_p(x + (size_t) j * K, K);
  return out;
}

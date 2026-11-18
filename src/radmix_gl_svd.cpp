// The SVD start of gt_radmix() for genotype likelihoods (R/radmix_svd.R,
// radmix_gl_svd_factors()). The genotypes are not observed, so the SVD is
// that of their standardised posterior dosages under Hardy-Weinberg
// proportions with the allele frequency of each locus, as in the first step
// of PCAngsd (Meisner and Albrechtsen 2018):
//   X_ij = (E[g_ij | D_ij, f_j] - 2 f_j) / sqrt(2 f_j (1 - f_j)),
// where g_ij is the number of copies of allele2, f_j its frequency at locus
// j, estimated from the likelihoods by EM (emMAF, Kim et al. 2011), and
// E[g | D, f] = (w1 + 2 w2) / (w0 + w1 + w2) with w0 = (1 - f)^2,
// w1 = 2 e1 f (1 - f) and w2 = e2 f^2 (the expected count g of slice_gl()
// in radmix_kernels.cpp, with qp = f). An individual without data at a
// locus has X_ij = 0 (its dosage is the prior mean 2 f_j), as a missing
// genotype imputed with the mean.
//
// X is never stored: the randomised SVD in R only needs its products with
// blocks of vectors, X B and X' B, each computed by one pass over the
// likelihoods of the source (decoded slice by slice, as in the kernels).
// Every output value is summed by one thread in a fixed order, so the
// results do not depend on the number of threads. Like the kernels, the
// passes run in OpenMP threads and never call R.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "radmix_geno_access.h"

using namespace Rcpp;
using namespace radmix;

namespace {

// the slices of individuals read from the source at once (as the kernels)
const size_t SLICE = 64;

// E[g | D, x]: the expected number of copies of allele2 given the
// likelihoods (1, e1, e2) and a frequency x of allele2
inline double gl_dosage(double e1, double e2, double x) {
  const double y = 1.0 - x;
  const double w0 = y * y, w1 = 2.0 * e1 * x * y, w2 = e2 * x * x;
  return (w1 + 2.0 * w2) / (w0 + w1 + w2);
}

// the standard deviation of the genotypes of a locus under Hardy-Weinberg
// proportions, 1 for a monomorphic locus (whose standardised dosages are
// then 0)
inline double dosage_scale(double f) {
  const double s = std::sqrt(2.0 * f * (1.0 - f));
  return s > 0 ? s : 1.0;
}

// the likelihoods of all the individuals at locus j, slice by slice
inline void read_locus(const GenoView& geno, size_t j, double* e1,
                       double* e2, double* obs) {
  const size_t I = geno.nrow();
  for (size_t i0 = 0; i0 < I; i0 += SLICE) {
    const size_t ni = std::min(SLICE, I - i0);
    geno.gl_slice(i0, ni, j, e1 + i0, e2 + i0, obs + i0);
  }
}

// checks that a view is that of a source of genotype likelihoods, without
// a cross-validation mask
void check_gl_view(const GenoView& geno, int n_cores) {
  if (!geno.is_gl())
    stop("These functions are for genotype likelihoods.");
  if (geno.cv_masked())
    stop("The SVD start does not take a cross-validation mask.");
  if (n_cores < 1) stop("'n_cores' must be at least 1.");
}

}  // namespace

// The frequency of allele2 at each locus, estimated from the genotype
// likelihoods by EM (emMAF): starting from 0.25, x is replaced by the mean
// of E[g | D, x] / 2 over the individuals with data, until it changes by
// less than tol (or after max_iter steps). A locus without data gets 0.5.
// [[Rcpp::export]]
NumericVector radmix_gl_freq_cpp(List data, int n_cores, double tol = 1e-6,
                                 int max_iter = 100) {
  const GenoView geno(data);
  check_gl_view(geno, n_cores);
  const size_t I = geno.nrow(), J = geno.ncol();
  NumericVector f(J);
  double* fp = f.begin();
#pragma omp parallel num_threads(n_cores)
{
  std::vector<double> e1(I), e2(I), obs(I);
#pragma omp for schedule(dynamic, 64)
  for (ptrdiff_t jj = 0; jj < (ptrdiff_t) J; jj++) {
    const size_t j = jj;
    read_locus(geno, j, e1.data(), e2.data(), obs.data());
    double n = 0.0;
    for (size_t i = 0; i < I; i++) n += obs[i];
    if (!(n > 0)) {
      fp[j] = 0.5;
      continue;
    }
    double x = 0.25;
    for (int it = 0; it < max_iter; it++) {
      double s = 0.0;
      for (size_t i = 0; i < I; i++)
        s += obs[i] * gl_dosage(e1[i], e2[i], x);
      const double x_new = s / (2.0 * n);
      const bool done = std::fabs(x_new - x) < tol;
      x = x_new;
      if (done) break;
    }
    fp[j] = x;
  }
}
  return f;
}

// The product of the standardised posterior dosages X (individuals x loci,
// see the top of the file) with b: X b (b: loci x l, result: individuals x
// l) or, with transpose, X' b (b: individuals x l, result: loci x l). f is
// the frequency of allele2 at each locus (radmix_gl_freq_cpp()).
// [[Rcpp::export]]
NumericMatrix radmix_gl_dosage_prod_cpp(List data, NumericVector f,
                                        NumericMatrix b, bool transpose,
                                        int n_cores) {
  const GenoView geno(data);
  check_gl_view(geno, n_cores);
  const size_t I = geno.nrow(), J = geno.ncol(), L = b.ncol();
  if ((size_t) f.size() != J)
    stop("'f' must have one value per locus.");
  if ((size_t) b.nrow() != (transpose ? I : J))
    stop("'b' must have one row per %s.",
         transpose ? "individual" : "locus");
  const double* fp = f.begin();
  NumericMatrix out(transpose ? J : I, L);
  double* op = out.begin();
  if (!transpose) {
    // X b: one task per slice of individuals, which accumulates the rows of
    // its individuals over all the loci; b is transposed first (l x J), so
    // that the l values of a locus are contiguous
    std::vector<double> bt(L * J);
    for (size_t j = 0; j < J; j++)
      for (size_t c = 0; c < L; c++) bt[c + L * j] = b(j, c);
    const ptrdiff_t n_slices = (ptrdiff_t) ((I + SLICE - 1) / SLICE);
#pragma omp parallel num_threads(n_cores)
{
    std::vector<double> e1(SLICE), e2(SLICE), obs(SLICE), acc(SLICE * L);
#pragma omp for schedule(dynamic, 1)
    for (ptrdiff_t sl = 0; sl < n_slices; sl++) {
      const size_t i0 = sl * SLICE, ni = std::min(SLICE, I - i0);
      std::fill(acc.begin(), acc.end(), 0.0);
      for (size_t j = 0; j < J; j++) {
        geno.gl_slice(i0, ni, j, e1.data(), e2.data(), obs.data());
        const double fj = fp[j], s = dosage_scale(fj);
        const double* bj = &bt[L * j];
        for (size_t ii = 0; ii < ni; ii++) {
          if (!(obs[ii] > 0)) continue;
          const double x = (gl_dosage(e1[ii], e2[ii], fj) - 2.0 * fj) / s;
          double* a = &acc[L * ii];
          for (size_t c = 0; c < L; c++) a[c] += x * bj[c];
        }
      }
      for (size_t ii = 0; ii < ni; ii++)
        for (size_t c = 0; c < L; c++)
          op[(i0 + ii) + I * c] = acc[L * ii + c];
    }
}
  } else {
    // X' b: one locus at a time, over all the individuals
    const double* bp = b.begin();
#pragma omp parallel num_threads(n_cores)
{
    std::vector<double> e1(I), e2(I), obs(I), x(I);
#pragma omp for schedule(dynamic, 64)
    for (ptrdiff_t jj = 0; jj < (ptrdiff_t) J; jj++) {
      const size_t j = jj;
      read_locus(geno, j, e1.data(), e2.data(), obs.data());
      const double fj = fp[j], s = dosage_scale(fj);
      for (size_t i = 0; i < I; i++)
        x[i] = obs[i] > 0 ? (gl_dosage(e1[i], e2[i], fj) - 2.0 * fj) / s
                          : 0.0;
      for (size_t c = 0; c < L; c++) {
        const double* bc = bp + I * c;
        double sum = 0.0;
        for (size_t i = 0; i < I; i++) sum += x[i] * bc[i];
        op[j + J * c] = sum;
      }
    }
}
  }
  return out;
}

#include <Rcpp.h>
#include "openadmixture_geno_access.h"
using namespace Rcpp;
using namespace openadmixture;

// Decode the selected subset into a dense matrix (NA = missing). Used only in
// tests to check that the accessor agrees with tidypopgen::show_genotypes().
// The fitting code never builds this dense matrix: the kernels decode each
// genotype on the fly from the memory-mapped file, so the memory used does
// not grow with the size of the dataset. Loci are the outer loop because the
// FBM is column-major (consecutive individuals of a locus are contiguous).
// The cross-validation arguments are as for openadmixture_geno_ptr() (cv_fold
// is 1-based), so that the tests can check which genotypes each fold hides.
// [[Rcpp::export]]
NumericMatrix openadmixture_decode_cpp(Environment BM,
                                       const IntegerVector& rowInd,
                                       const IntegerVector& colInd,
                                       const NumericVector& lut,
                                       const NumericVector& n_trials,
                                       const NumericVector& scale,
                                       int cv_mode = 0, int cv_fold = 1,
                                       int n_folds = 1, double cv_seed = 0) {
  GenoView geno(BM, rowInd, colInd, lut, n_trials, scale, cv_mode,
                cv_fold - 1, n_folds, (uint64_t) cv_seed);
  size_t n = geno.nrow(), m = geno.ncol();
  NumericMatrix out(n, m);
  double g;
  for (size_t j = 0; j < m; j++)
    for (size_t i = 0; i < n; i++)
      out(i, j) = geno.get(i, j, g) ? g : NA_REAL;
  return out;
}

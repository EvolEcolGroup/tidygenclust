// An example of a source of genotype likelihoods for gt_radmix(): the
// likelihoods of every individual and locus held in memory, as dense
// matrices. tidygenclust does not store genotype likelihoods itself: a
// package that does implements the same
// interface (inst/include/tidygenclust/radmix_gl_source.h) on its own
// storage, and returns the source from its method of gt_radmix_data(). This
// one is used by the tests (radmix_test_gl() in
// tests/testthat/helper-radmix-data.R), so that the kernels can be checked
// on genotype likelihoods without any other package, and it shows what a
// source has to provide: a context, created and checked here in R's main
// thread, and a function that only copies the likelihoods of a slice out of
// it (called by the kernels from any thread).

#include <Rcpp.h>
#include <tidygenclust/radmix_gl_source.h>
#include <cmath>
#include <vector>

// The context: the likelihoods e1 and e2 (relative to 0 copies of allele2)
// and the indicator of data of each individual (rows) and locus (columns),
// column-major as in R. They are copied out of R memory, so that the
// context does not depend on R objects that the garbage collector could
// free while the source is in use.
struct RadmixDenseGl {
  size_t n_indiv;
  std::vector<double> e1, e2, obs;
};

// The function of the source (see radmix_gl_source.h): the likelihoods of
// individuals i0 to i0 + ni - 1 at locus j, which are adjacent in each
// column. It has C linkage, as the type of the function pointer.
extern "C" void radmix_gl_dense_fill(const void* ctx, size_t i0, size_t ni,
                                     size_t j, double* e1, double* e2,
                                     double* obs) {
  const RadmixDenseGl* d = static_cast<const RadmixDenseGl*>(ctx);
  const size_t off = j * d->n_indiv + i0;
  for (size_t ii = 0; ii < ni; ii++) {
    e1[ii] = d->e1[off + ii];
    e2[ii] = d->e2[off + ii];
    obs[ii] = d->obs[off + ii];
  }
}

// The function and the context of a source of the likelihoods e1 and e2
// (individuals x loci, positive and finite) with the indicator of data obs
// (TRUE or FALSE), as the elements fill and ctx of the list that a method
// of gt_radmix_data() returns. The context is freed when the garbage
// collector frees its external pointer.
// [[Rcpp::export]]
Rcpp::List radmix_gl_dense_source_cpp(const Rcpp::NumericMatrix& e1,
                                      const Rcpp::NumericMatrix& e2,
                                      const Rcpp::LogicalMatrix& obs) {
  const size_t n = e1.nrow(), m = e1.ncol();
  if ((size_t) e2.nrow() != n || (size_t) e2.ncol() != m ||
      (size_t) obs.nrow() != n || (size_t) obs.ncol() != m)
    Rcpp::stop("'e1', 'e2' and 'obs' must have the same dimensions.");
  // the checks that the kernels rely on (see radmix_gl_source.h): positive,
  // finite likelihoods, and data or no data
  for (size_t x = 0; x < n * m; x++) {
    if (!(e1[x] > 0 && e2[x] > 0 && std::isfinite(e1[x]) &&
          std::isfinite(e2[x])))
      Rcpp::stop("'e1' and 'e2' must be positive and finite.");
    if (obs[x] == NA_LOGICAL) Rcpp::stop("'obs' must not be NA.");
  }
  Rcpp::XPtr<RadmixDenseGl> ctx(new RadmixDenseGl, true);
  ctx->n_indiv = n;
  ctx->e1.assign(e1.begin(), e1.end());
  ctx->e2.assign(e2.begin(), e2.end());
  ctx->obs.assign(obs.begin(), obs.end());
  // the function, tagged so that tidygenclust can recognise it
  Rcpp::RObject fill = R_MakeExternalPtrFn(
    reinterpret_cast<DL_FUNC>(&radmix_gl_dense_fill),
    Rf_install(TIDYGENCLUST_GL_FILL_TAG), R_NilValue);
  return Rcpp::List::create(Rcpp::Named("fill") = fill,
                            Rcpp::Named("ctx") = ctx);
}

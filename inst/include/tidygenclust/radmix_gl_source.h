// Genotype likelihoods from another package, for gt_radmix() (version 1).
//
// gt_radmix() fits the model of NGSadmix (Skotte et al. 2013) to genotype
// likelihoods, but tidygenclust does not store them: a package that does
// gives them to the kernels of gt_radmix() through a "source", a C function that fills the likelihoods of a block of
// individuals at one locus, and its context (the data that the function
// reads). The package provides a method of the R generic gt_radmix_data()
// for its class, which returns a list with the elements
//   kind = "gl_source";
//   version = 1L (TIDYGENCLUST_GL_SOURCE_VERSION);
//   fill: an external pointer to the function (made with
//     R_MakeExternalPtrFn()), whose tag is the symbol tidygenclust_gl_fill
//     (TIDYGENCLUST_GL_FILL_TAG), so that tidygenclust can tell it from
//     another external pointer;
//   ctx: an external pointer to the context, passed to the function as is;
//   n_indiv, n_loci: the numbers of individuals and loci, in the order of
//     the rows and of the loci (show_loci()) of the object;
// and any other elements that must stay alive during the fit (e.g. the
// file-backed matrix whose mapping the context points into): tidygenclust
// keeps the list, and hence them, for as long as it uses the source.
//
// This header has no dependency on R or on tidygenclust: such a package can
// copy the typedef below (and the two constants), rather than LinkingTo
// tidygenclust. tidygenclust checks the version and the tag at run time.
//
// The function. fill(ctx, i0, ni, j, e1, e2, obs) writes, for the
// individuals i0, ..., i0 + ni - 1 (0-based, in the order of the rows of
// the object; ni >= 1) at locus j (0-based, in the order of its loci), into
// e1[0..ni), e2[0..ni) and obs[0..ni):
//   e1[ii], e2[ii]: the likelihoods of the genotypes allele1/allele2 and
//     allele2/allele2, relative to that of allele1/allele1 (i.e. the
//     likelihoods of 0, 1 and 2 copies of allele2 are proportional to
//     1, e1, e2); both must be positive and finite (a likelihood of exactly
//     0 would make some estimates impossible: clamp it to a small ratio
//     instead), also for an individual without data, which can use 1 and 1;
//   obs[ii]: 1 if the individual has data at the locus, 0 if not (it then
//     contributes nothing to the log-likelihood).
// The estimates of P are the frequencies of allele2 (often the minor
// allele). The function is called many times per pass over the data,
// from several OpenMP threads at once, so it must be thread safe (read
// only), must not call R (no allocation from R, no Rcpp, no R errors), must
// not throw C++ exceptions, and should be fast (it is called once per
// block of individuals and locus, out of line, so a plain loop over the
// individuals is enough). Check everything (indices, types, sizes) when the
// context is created, in R's main thread.
#ifndef TIDYGENCLUST_RADMIX_GL_SOURCE_H
#define TIDYGENCLUST_RADMIX_GL_SOURCE_H

#include <stddef.h>

// the version of this interface, which the list must give as `version`
#define TIDYGENCLUST_GL_SOURCE_VERSION 1
// the name of the symbol that must tag the external pointer to the function
#define TIDYGENCLUST_GL_FILL_TAG "tidygenclust_gl_fill"

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*tidygenclust_gl_fill_fn)(const void* ctx, size_t i0,
                                        size_t ni, size_t j, double* e1,
                                        double* e2, double* obs);

#ifdef __cplusplus
}
#endif

#endif

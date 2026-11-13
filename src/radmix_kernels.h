#ifndef RADMIX_KERNELS_H
#define RADMIX_KERNELS_H

// The passes over the genotypes of gt_radmix() (radmix_kernels.cpp).
//
// The model
// ---------
// For individual i and locus j, the genotype g_ij (the number of copies of
// the alternate allele) is binomial,
//   g_ij ~ Binomial(n_i, qp_ij),   qp_ij = sum_k q_ki p_kj,
// where n_i is the number of allele copies (2 for diploids, 1 for
// pseudohaploids), q_ki the proportion of the genome of individual i that
// comes from ancestral population k, and p_kj the frequency of the
// alternate allele at locus j in population k. The log-likelihood, over the
// observed genotypes and without the binomial coefficients (which do not
// depend on q or p), is
//   L(q, p) = sum_ij g_ij log(qp_ij) + (n_i - g_ij) log(1 - qp_ij).
// Each column of q lies on the simplex (q_ki >= 0, sum_k q_ki = 1), and each
// p_kj in [0, 1]. Given p, L is a sum of one term per individual, and given
// q, one term per locus, so the passes work one column at a time: Q-side
// passes sum over the loci of each individual, P-side passes over the
// individuals of each locus.
//
// Conventions
// -----------
// q is K x I (column i holds the K ancestry proportions of individual i) and
// p is K x J (column j holds the K allele frequencies of locus j), both
// column-major, so the K values updated together are contiguous. Hessians
// are K x K per column (K x K x I or K x K x J), and are the NEGATIVE
// second derivatives of L (positive semi-definite).
//
// Common arguments:
//   geno     the genotypes (I individuals x J loci), see radmix_geno_access.h
//   q, p     the current estimates (K x I and K x J), read only
//   K        the number of ancestral populations
//   n_cores  the number of OpenMP threads (results are identical for any)
// Outputs are allocated by the caller, with the sizes given below, and are
// overwritten; they must not overlap q or p. Nothing is checked here: the
// callers (the R wrappers in radmix_exports.cpp and the fitting loop in
// radmix_fit.cpp) check the dimensions, and require I, J >= 1. None of these
// functions calls R, and their results do not depend on the number of
// threads.

#include <cstddef>
#include "radmix_geno_access.h"

namespace radmix {

// The log-likelihood L(q, p) (-Inf if a genotype is impossible under q and
// p, i.e. qp of exactly 0 with g > 0, or of exactly 1 with g < n_i).
double loglik(const GenoView& geno, const double* q, const double* p, int K,
              int n_cores);

// One EM step, updating both blocks from the same (q, p): q_new (K x I) and
// p_new (K x J). The result is not projected: each column of q_new is on
// the simplex, but p_new can be exactly 0 or 1 (see em()).
void em(const GenoView& geno, const double* q, const double* p, int K,
        double* q_new, double* p_new, int n_cores);

// The gradient (K x I) and negative Hessian (K x K x I) of L with respect to
// each column of q.
void grad_q(const GenoView& geno, const double* q, const double* p, int K,
            double* grad, double* hess, int n_cores);

// The gradient (K x J) and negative Hessian (K x K x J) of L with respect to
// each column of p.
void grad_p(const GenoView& geno, const double* q, const double* p, int K,
            double* grad, double* hess, int n_cores);

// The block update of q given p (q_new, K x I): each column takes the step
// that maximises the quadratic approximation of L subject to staying on the
// simplex, and is then projected onto {sum = 1, q >= 1e-5}.
void update_q(const GenoView& geno, const double* q, const double* p, int K,
              double* q_new, int n_cores);

// update_q() and the log-likelihood L(q, p), which is returned, in a single
// pass over the genotypes: the Q-side pass computes qp for every genotype,
// so the log-likelihood of its input costs little more (the fitting loop
// uses it to evaluate a point and start the next iteration from it at the
// same time). q_new is identical to that of update_q(), and the
// log-likelihood equals that of loglik() up to rounding (the same terms,
// added up in another order); neither depends on the number of threads.
double update_q_loglik(const GenoView& geno, const double* q,
                       const double* p, int K, double* q_new, int n_cores);

// The block update of p given q (p_new, K x J): the same, within the box
// [0, 1], then clamped to [1e-5, 1 - 1e-5].
void update_p(const GenoView& geno, const double* q, const double* p, int K,
              double* p_new, int n_cores);

// The binomial deviance of the visible genotypes (*dev) and their number
// (*count), the prediction error of cross-validation, computed with a
// held-out view of the genotypes (see radmix_geno_access.h).
void deviance(const GenoView& geno, const double* q, const double* p, int K,
              int n_cores, double* dev, double* count);

// The number of segments of loci into which the Q-side passes split their
// work for I individuals and J loci (1 = none). It depends only on the
// dimensions of the data, never on the number of threads.
size_t q_segments(size_t I, size_t J);

}  // namespace radmix

#endif

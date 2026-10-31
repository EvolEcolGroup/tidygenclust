#ifndef RADMIX_KERNELS_H
#define RADMIX_KERNELS_H

// Port of the CPU kernels of OpenADMIXTURE.jl (src/loops.jl and the inner
// algorithms in src/algorithms_inner.jl), reading genotypes through GenoView.
//
// All matrices are plain column-major arrays: q is K x I (column i holds the
// ancestry proportions of individual i), p is K x J (column j holds the alt
// allele frequencies of locus j), Hessians are K x K x I or K x K x J.
// Missing genotypes are skipped, and the constant 2 of the Julia code is
// replaced by the number of trials n_i of each individual (2 for diploids,
// 1 for pseudohaploids).
//
// None of these functions touches the R API, so they can be called from other
// C++ code (e.g. the outer loop) without going back to R. Results do not
// depend on the number of threads: each column is always summed in the same
// order, and there are no cross-thread reductions.
//
// The model, in brief (see radmix_kernels.cpp for the derivations):
//   qp_ij = sum_k q_ki p_kj                  (expected alt allele frequency)
//   g_ij ~ Binomial(n_i, qp_ij)
//   L(q, p) = sum_{i,j observed} g_ij log(qp_ij) + (n_i - g_ij) log(1 - qp_ij)
// (the binomial coefficient does not depend on q or p and is dropped).
// q is constrained to the simplex (q_ki >= 0, sum_k q_ki = 1) and p to the
// box [0, 1]. Because L separates over individuals once p is fixed, and over
// loci once q is fixed, every kernel below works one column at a time:
// Q-side kernels sum over loci for each individual, P-side kernels sum over
// individuals for each locus.
//
// Common arguments:
//   geno     genotype accessor (I individuals x J loci), see
//            radmix_geno_access.h
//   q, p     current estimates, K x I and K x J, column-major, read only
//   K        number of ancestral populations
//   n_cores  number of OpenMP threads (results are identical for any value)
// Output buffers are allocated by the caller with the sizes stated below and
// are completely overwritten. They must not overlap q or p: the kernels are
// not written to work in place (em(), for instance, reads q while it
// normalises q_new).

#include "radmix_geno_access.h"

// internal code of gt_radmix() (see the note on the namespace in
// radmix_geno_access.h)
namespace radmix {


// log-likelihood L(q, p) as defined above (a scalar, without the binomial
// constant); Julia: loglikelihood() in algorithms_inner.jl
double loglik(const GenoView& geno, const double* q, const double* p, int K,
              int n_cores);

// one FRAPPE EM step (q_new and p_new are computed from the same q and p)
// q_new is K x I, p_new is K x J; Julia: em_q!() and em_p!() in "base" mode
void em(const GenoView& geno, const double* q, const double* p, int K,
        double* q_new, double* p_new, int n_cores);

// gradient (K x I) and negative Hessian (K x K x I) of the log-likelihood
// with respect to q
// (in the homogeneous form used by Julia's update_q_loop!(); see
// radmix_kernels.cpp)
void grad_q(const GenoView& geno, const double* q, const double* p, int K,
            double* grad, double* hess, int n_cores);

// gradient (K x J) and negative Hessian (K x K x J) of the log-likelihood
// with respect to p
// (as in Julia's update_p_loop!())
void grad_p(const GenoView& geno, const double* q, const double* p, int K,
            double* grad, double* hess, int n_cores);

// one constrained Newton (QP) update of q (q_new is K x I), projected
// For each individual, the step delta maximises the quadratic approximation
// of L subject to q + delta lying on the simplex, then q + delta is
// projected onto {sum = 1, q >= eps}; Julia: update_q!()
void update_q(const GenoView& geno, const double* q, const double* p, int K,
              double* q_new, int n_cores);

// one constrained Newton (QP) update of p (p_new is K x J), projected
// Same as above for each locus, with the box constraint 0 <= p + delta <= 1
// and a final clamp to [eps, 1 - eps]; Julia: update_p!()
void update_p(const GenoView& geno, const double* q, const double* p, int K,
              double* p_new, int n_cores);

// binomial deviance of the visible genotypes (*dev) and their number
// (*count), for ADMIXTURE-style cross-validation with a held-out view (see
// radmix_geno_access.h); no Julia counterpart (OpenADMIXTURE.jl has no
// cross-validation)
void deviance(const GenoView& geno, const double* q, const double* p, int K,
              int n_cores, double* dev, double* count);

}  // namespace radmix

#endif

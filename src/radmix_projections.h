#ifndef RADMIX_PROJECTIONS_H
#define RADMIX_PROJECTIONS_H

// internal code of gt_radmix() (see the note on the namespace in
// radmix_geno_access.h)
namespace radmix {


// Port of OpenADMIXTURE.jl src/projections.jl.
//
// project_q: Euclidean projection of q[0..K-1] onto the shifted simplex
//   {sum(q) == 1, q >= eps}, in place.
// project_p: clamp p[0..K-1] to [eps, 1 - eps], in place.
// Both are thread-safe (no shared state, no R API).
// Unlike the Julia versions, both tolerate non-finite input (a deviation
// from OpenADMIXTURE.jl): project_q resets a column containing NaN or Inf
// to 1 / K, and project_p maps NaN to 0.5 (see radmix_projections.cpp).
//
// They are applied after every Newton step (src/radmix_kernels.cpp, on
// par + delta from the QP solver), after the quasi-Newton extrapolation
// (src/radmix_quasi_newton.cpp) and to the starting values (via
// radmix_project_q_cpp in R). The QP already keeps columns inside [0,
// 1], but the projections tidy up rounding error, any partial step left by a
// failed QP, and quasi-Newton points that can land anywhere.

// Lower bound for every entry of Q and P (and 1 - RADMIX_EPS the upper
// bound for P). 1e-5 matches the output bounds of ADMIXTURE and
// OpenADMIXTURE.jl, and keeps the log-likelihood finite: with q and p bounded
// away from 0 and 1, the per-genotype probability sum_k q_ik p_kj and its
// complement are strictly positive, so log() never sees 0.
const double RADMIX_EPS = 1e-5;

void project_q(double* q, int K, double eps = RADMIX_EPS);
void project_p(double* p, int K, double eps = RADMIX_EPS);

}  // namespace radmix

#endif

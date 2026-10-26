#ifndef OPENADMIXTURE_QP_SWEEP_H
#define OPENADMIXTURE_QP_SWEEP_H

// Port of OpenADMIXTURE.jl src/qp.jl: K. Lange's sweep-operator active-set
// quadratic program, used for the per-column Newton steps of Q and P.
//
// For one column (an individual for Q, a locus for P) with current value
// par[0..K-1], gradient of the log-likelihood grad[0..K-1] and NEGATIVE
// Hessian H (K x K, column-major, symmetric positive semi-definite), find
//   delta = argmin 0.5 * delta' H delta - grad' delta
// subject to 0 <= par + delta <= 1 and, if simplex is true,
//   sum(delta) == 1 - sum(par).
// On failure delta is left as far as the solver got (0 if it failed during
// the initial sweeps), exactly as in the Julia code.
//
// The objective is the second-order Taylor expansion of the log-likelihood
// around par, with its sign flipped so that we minimise; without the
// constraints its minimiser is the plain Newton step H^{-1} grad. The
// constraints keep the updated column a valid set of ancestry proportions
// (simplex) or allele frequencies (box).
//
// Callers: update_q (simplex = true, one call per individual) and
// update_p (simplex = false, one call per locus) in
// src/openadmixture_kernels.cpp, which add delta to the current column and then
// project the result with project_q / project_p
// (src/openadmixture_projections.h). Both callers ignore the return value and
// use whatever delta the solver produced.
//
// QPWork holds all scratch memory for a given K; create one per thread.
// Nothing here touches the R API, so it is safe inside OpenMP regions.

#include <vector>

// internal code of the OpenADMIXTURE implementation (see the note on the
// namespace in openadmixture_geno_access.h)
namespace openadmixture {


// Scratch space for qp_solve_column. Everything is sized once for the
// larger (simplex, Q) case so that the same object can serve both Q and P
// columns with no allocation inside the hot loop. The solver overwrites the
// contents on every call, so an instance must never be shared between
// threads: the kernels create one per OpenMP thread.
struct QPWork {
  // number of ancestral populations (the dimension of each QP)
  int K;
  // The tableau, stored column-major with leading dimension n, where
  // n = K + 2 for Q and n = K + 1 for P. Only the first n * n entries are
  // used, so for P the tableau is packed at the start of the buffer rather
  // than being the top-left block of a (K+2) x (K+2) matrix.
  std::vector<double> tab;     // tableau, (K+2) x (K+2) column-major max
  // The tableau diagonal before any sweep; the positive-definiteness test
  // compares each pivot against its original value.
  std::vector<double> diag;    // saved diagonal
  // swept[k] != 0 means parameter k is currently "free" (its row and
  // column of the tableau have been swept, so it moves with the Newton
  // step); swept[k] == 0 means it is held fixed on a bound. Only the first
  // K entries are used; the extra two keep the sizes uniform.
  std::vector<char> swept;     // swept flags
  // An orthonormal basis of R^K whose first column is ones(K) / sqrt(K),
  // used only to compute the penalty mu for the simplex case. It depends
  // only on K, so it is built once in the constructor.
  std::vector<double> V;       // K x K right singular vectors of ones(1, K)
  // Scratch vector: the saved pivot column during a sweep (length n, hence
  // K + 2 entries) and H * V[:, 1] while computing mu (length K).
  std::vector<double> tmp;     // scratch, K + 2 entries
  explicit QPWork(int K);
};

// Solve the column QP described above. H is K x K column-major, grad and
// par have length K, delta (length K) receives the step. simplex selects
// the Q problem (box plus sum constraint) rather than the P problem (box
// only).
// returns the number of active-set iterations used, or 0 on failure
// (Hessian not numerically positive definite, or no convergence within
// 1000 iterations; in the latter case delta may still be a good step)
int qp_solve_column(const double* H, const double* grad, const double* par,
                    int K, bool simplex, double* delta, QPWork& work);

}  // namespace openadmixture

#endif

#ifndef RADMIX_QP_SWEEP_H
#define RADMIX_QP_SWEEP_H

// The quadratic program of the Newton step of one column of q or p, solved
// exactly by an active-set method based on the sweep operator (Lange,
// "Numerical Analysis for Statisticians"); see radmix_qp_sweep.cpp.
//
// For one column (an individual for q, a locus for p) with current value
// par (K values), gradient of the log-likelihood grad (K values) and negative
// Hessian H (K x K, column-major, symmetric positive semi-definite), find
//   delta = argmin 0.5 delta' H delta - grad' delta
// subject to 0 <= par + delta <= 1 and, if simplex is true (a column of q),
//   sum(delta) == 1 - sum(par).
// The objective is minus the second-order Taylor expansion of the
// log-likelihood around par, so without the constraints its minimiser is the
// Newton step H^-1 grad; the constraints keep the column a valid set of
// ancestry proportions (simplex) or allele frequencies (box).
//
// Called by newton_column() (radmix_kernels.cpp) for every column of q and
// p, inside parallel regions: nothing here calls R, and all the state lives
// in the caller's QPWork. (radmix_qp_solve_cpp(), in radmix_exports.cpp,
// calls it from R for the tests.)

#include <vector>

namespace radmix {

// Scratch memory for qp_solve_column(), sized for the larger problem (a
// column of q), so that one object serves both blocks without allocations.
// It is overwritten on every call, so each thread needs its own.
struct QPWork {
  // the number of ancestral populations (the dimension of each problem)
  int K;
  // the tableau, column-major with leading dimension n (n = K + 2 for q,
  // K + 1 for p; only the first n * n entries are used)
  std::vector<double> tab;
  // the diagonal of the tableau before any sweep, against which the pivots
  // are compared
  std::vector<double> diag;
  // swept[k] != 0: parameter k is free (its row and column have been swept,
  // and it moves with the Newton step); 0: it is held on a bound
  std::vector<char> swept;
  // the saved pivot column of a sweep (n entries)
  std::vector<double> tmp;
  // the step of the column being updated (K entries), used by
  // newton_column()
  std::vector<double> delta;
  // the number of passes of the active-set loop in the last call, 0 if the
  // initial sweeps failed (used by the tests; unlike the return value of
  // qp_solve_column(), it is also set when the loop stops early)
  int passes = 0;
  // allocate the scratch for problems with K parameters
  explicit QPWork(int K);
};

// Solve the problem described above: H (K x K), grad and par (K values) are
// read only, delta (K values) receives the step, and work must have been
// built for at least K parameters.
// Returns the number of iterations of the active-set loop, or 0 if the
// initial sweeps fail (the tableau not numerically positive definite, i.e.
// H for p, or H + mu 11' for q, see radmix_qp_sweep.cpp: delta = 0), if the
// active set cycles at a vertex of the simplex (delta is then the
// solution), or if the loop does not settle within 1000 iterations (delta is
// then as far as it got, usually the solution). The kernels therefore use
// delta whatever the return value.
int qp_solve_column(const double* H, const double* grad, const double* par,
                    int K, bool simplex, double* delta, QPWork& work);

}  // namespace radmix

#endif

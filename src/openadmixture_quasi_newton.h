#ifndef OPENADMIXTURE_QUASI_NEWTON_H
#define OPENADMIXTURE_QUASI_NEWTON_H

// Quasi-Newton acceleration of Zhou, Alexander and Lange (2011); see
// openadmixture_quasi_newton.cpp for the method. These are plain C++ functions
// on raw pointers, called by the fitting loop (openadmixture_fit.cpp);
// openadmixture_quasi_newton.cpp also exports thin Rcpp wrappers around them
// for the tests.
//
// Parameters are conceptually stacked as x = (q, p): the K x I values of q
// followed by the K x J values of p, both column-major. They are always
// passed as the two separate blocks, never stacked in memory.

#include <cstddef>
#include <vector>

// internal code of the OpenADMIXTURE implementation (see the note on the
// namespace in openadmixture_geno_access.h)
namespace openadmixture {


// The last n_qn secant pairs, as the columns of U and V (n_par x n_qn each,
// column-major, n_par = K (I + J)). Together they take 2 n_qn copies of the
// parameters (6 with the default n_qn = 3), more than the four sets of
// estimates kept by the fitting loop.
struct QNState {
  size_t n_par, n_qn;
  std::vector<double> U, V;
  QNState(size_t n_par, size_t n_qn)
    : n_par(n_par), n_qn(n_qn), U(n_par * n_qn, 0.0), V(n_par * n_qn, 0.0) {}
};

// Store u = x1 - x0 and v = x2 - x1 in column col (0-based) of U and V, with
// x0 = (q, p), x1 = (q1, p1), x2 = (q2, p2); nq = K I and np = K J are the
// sizes of the q and p blocks.
void qn_push(QNState& st, size_t col, const double* q, const double* p,
             const double* q1, const double* p1, const double* q2,
             const double* p2, size_t nq, size_t np);

// The quasi-Newton step from x1 = (q1, p1), using the first n_used pairs,
// written to q_out (K x n_indiv) and p_out (K x n_loci) and projected onto
// the feasible set. Returns false (leaving the outputs undefined) if the
// n_used x n_used system is singular to machine precision (or the system
// or its solution is not finite). The outputs are written only after the
// last read of q and p.
bool qn_step(const QNState& st, size_t n_used, const double* q,
             const double* p, const double* q1, const double* p1, int K,
             size_t n_indiv, size_t n_loci, double* q_out, double* p_out,
             int n_cores);

}  // namespace openadmixture

#endif

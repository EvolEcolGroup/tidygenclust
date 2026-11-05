#ifndef RADMIX_QUASI_NEWTON_H
#define RADMIX_QUASI_NEWTON_H

// Quasi-Newton acceleration of the block updates (Zhou, Alexander and Lange
// 2011); see radmix_quasi_newton.cpp.
//
// The parameters x = (q, p) are never stacked into one vector: every
// function takes the q block (nq = K I values) and the p block (np = K J
// values) separately, and positions 0 to nq - 1 of a secant vector refer to
// the q block, positions nq to nq + np - 1 to the p block.

#include <cstddef>
#include <vector>

namespace radmix {

// The last n_qn secant pairs, as the columns of U and V (n_par x n_qn,
// column-major, n_par = nq + np). Column col holds the pair of one iteration
// of the fitting loop: u = x1 - x0 and v = x2 - x1. The fitting loop uses
// the columns as a circular buffer (iteration it writes column
// (it - 1) mod n_qn); the order of the columns does not matter to
// qn_step().
struct SecantPairs {
  // the number of parameters (rows) and of pairs (columns)
  size_t n_par, n_qn;
  // the secant vectors u (U) and v (V), n_par x n_qn each, zero until
  // stored
  std::vector<double> U, V;
  SecantPairs(size_t n_par, size_t n_qn)
    : n_par(n_par), n_qn(n_qn), U(n_par * n_qn, 0.0), V(n_par * n_qn, 0.0) {}
  // store the pair of the estimates x0 = (q, p), x1 = (q1, p1) and
  // x2 = (q2, p2) in column col (< n_qn), overwriting the pair stored there;
  // nq and np are the sizes of the q and p blocks (nq + np = n_par)
  void store(size_t col, const double* q, const double* p, const double* q1,
             const double* p1, const double* q2, const double* p2,
             size_t nq, size_t np);
};

// The quasi-Newton point x1 - V (U'(U - V))^-1 U'(x0 - x1), from the first m
// columns of U and V (n_par x m or larger, column-major), projected onto the
// feasible set, written into (q_out, p_out). x0 = (q, p) and x1 = (q1, p1),
// with q blocks of K x n_indiv values and p blocks of K x n_loci values
// (n_par = K (n_indiv + n_loci)); q_out and p_out may be q and p (x0 is
// only read before the output is written). n_cores is the number of OpenMP
// threads; the result does not depend on it. Returns false, without writing
// anything, if the m x m system is singular to machine precision or not
// finite.
bool qn_step(const double* U, const double* V, size_t n_par, size_t m,
             const double* q, const double* p, const double* q1,
             const double* p1, int K, size_t n_indiv, size_t n_loci,
             double* q_out, double* p_out, int n_cores);

}  // namespace radmix

#endif

#ifndef RADMIX_PROJECTIONS_H
#define RADMIX_PROJECTIONS_H

// Projections of a column of q or p onto its feasible set, applied after
// every Newton step (radmix_kernels.cpp), after the quasi-Newton step
// (radmix_quasi_newton.cpp) and to the starting values given by the user
// (radmix_project_cpp(), from R). The Newton steps already stay within
// [0, 1], but the projections remove rounding errors, and the quasi-Newton
// point, which ignores the constraints, can land anywhere.
//
// project_q: the Euclidean projection of q[0..K-1] onto the shifted simplex
//   {sum(q) == 1, q >= eps}, in place.
// project_p: clamps p[0..K-1] to [eps, 1 - eps], in place (K can be any
//   number of values).
// Both tolerate non-finite values (see radmix_projections.cpp), have no
// shared state and never call R, so they can run inside parallel regions.

namespace radmix {

// The lower bound of every entry of q and p (and 1 - RADMIX_EPS the upper
// bound of p), as in ADMIXTURE. With q and p bounded away from 0 and 1, qp
// and 1 - qp are strictly positive, so the log-likelihood is always finite.
const double RADMIX_EPS = 1e-5;

void project_q(double* q, int K, double eps = RADMIX_EPS);
void project_p(double* p, int K, double eps = RADMIX_EPS);

}  // namespace radmix

#endif

// Port of OpenADMIXTURE.jl src/projections.jl (project_q!, project_p!).
#include "openadmixture_projections.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

// internal code of the OpenADMIXTURE implementation (see the note on the
// namespace in openadmixture_geno_access.h)
namespace openadmixture {


// Projection onto {sum(q) == 1, q >= eps}: shift by eps, project onto the
// simplex of mass tau = 1 - K * eps with the sort-and-threshold algorithm,
// then shift back. Julia sorts an index permutation with its own quicksort;
// only the sorted values matter, so we sort a copy of the values instead.
//
// Why this works. Substituting x = q - eps turns the target set into
// {sum(x) == tau, x >= 0}. The closest point of that set to a vector b is
// x_k = max(b_k - t, 0) for the single threshold t at which these entries
// add up to tau (this follows from the KKT conditions: every entry is
// shifted down by the same Lagrange multiplier t and negative results are
// clipped to 0). So the whole problem reduces to finding t.
//
// With b sorted in decreasing order, b_(1) >= b_(2) >= ..., suppose exactly
// the m largest entries stay positive. Then sum_{i <= m} (b_(i) - t) = tau
// gives t_m = (b_(1) + ... + b_(m) - tau) / m, and this guess is right when
// the next entry would be clipped, i.e. t_m >= b_(m+1). Scanning m = 1, 2,
// ... and stopping at the first m that passes this test yields the correct
// t; if no m < K passes, all K entries stay positive and t = t_K.
//
// Worked example (eps = 0 for clarity, tau = 1): b = (0.9, 0.5, -0.2).
// m = 1: t = (0.9 - 1) / 1 = -0.1, which is below b_(2) = 0.5, keep going.
// m = 2: t = (1.4 - 1) / 2 = 0.2, which is >= b_(3) = -0.2, stop.
// Result: (0.9 - 0.2, 0.5 - 0.2, 0) = (0.7, 0.3, 0), summing to 1.
//
// Cost is O(K log K) for the sort; K is small, so this is negligible.
void project_q(double* q, int K, double eps) {
  // The sorted copy lives on the stack for the usual small K, avoiding a
  // heap allocation per individual; very large K falls back to the heap.
  const int stack_max = 64;
  double stack_buf[stack_max];
  std::vector<double> heap_buf;
  double* b = stack_buf;
  if (K > stack_max) {
    heap_buf.resize(K);
    b = heap_buf.data();
  }
  // Deviation from OpenADMIXTURE.jl: guard against non-finite input.
  // Julia's project_q! assumes finite values. Here a NaN would reach
  // std::sort with std::greater, which requires a strict weak ordering that
  // NaN breaks; that is undefined behaviour and, in practice, can make the
  // sort read outside the buffer. An infinite value would sort fine but turn
  // the threshold, and hence every entry, into NaN. Neither should happen
  // (the kernels discard non-finite Newton steps and qn_step() discards
  // non-finite quasi-Newton coefficients before projecting), so this is a last
  // line of defence: rather than propagate garbage, the column is reset to
  // the centre of the simplex (1 / K each), a valid point from which the
  // following iterations can recover.
  for (int k = 0; k < K; k++) {
    if (!std::isfinite(q[k])) {
      for (int k2 = 0; k2 < K; k2++) q[k2] = 1.0 / K;
      return;
    }
  }
  // mass left once every entry has been given its eps floor
  const double tau = 1.0 - K * eps;
  // shift to x = q - eps (in place) and take a copy to sort
  for (int k = 0; k < K; k++) {
    q[k] -= eps;
    b[k] = q[k];
  }
  // decreasing order, as the threshold scan expects
  std::sort(b, b + K, std::greater<double>());
  // find the threshold: the first i at which the running mean exceeds the
  // next largest value; otherwise use all K entries
  //
  // In the notation above, after adding b[i] the loop has m = i + 1 kept
  // entries, tsum = b_(1) + ... + b_(m) and tmax = t_m; it stops as soon as
  // t_m >= b_(m+1) (b[i + 1] in 0-based indexing).
  double tsum = 0, tmax = 0;
  bool found = false;
  for (int i = 0; i < K - 1; i++) {
    tsum += b[i];
    tmax = (tsum - tau) / (i + 1);
    if (tmax >= b[i + 1]) {
      found = true;
      break;
    }
  }
  // no early stop: all K entries are kept and t = t_K (this is also the
  // only branch taken for K == 1, giving q = 1). This reads b[K - 1], so
  // K >= 1 is required (K == 0 would read before the buffer).
  if (!found) tmax = (tsum + b[K - 1] - tau) / K;
  // apply the threshold, clip at 0 and undo the eps shift; the result sums
  // to tau + K * eps = 1 with every entry >= eps
  for (int k = 0; k < K; k++) q[k] = std::max(q[k] - tmax, 0.0) + eps;
}

// Allele frequencies are only box-constrained, and the box is a product of
// intervals, so the Euclidean projection is simply entry-wise clamping to
// [eps, 1 - eps]. K here is just the number of entries, so callers may pass
// a whole matrix at once.
//
// Deviation from OpenADMIXTURE.jl: NaN entries are set to 0.5. Clamping
// already maps -Inf to eps and +Inf to 1 - eps, but std::max and std::min
// let NaN through unchanged (every comparison with NaN is false), and a NaN
// frequency would make the log-likelihood NaN for the whole locus. As in
// project_q, this should never trigger, and 0.5 (the centre of the box) is
// simply a valid value to restart from.
void project_p(double* p, int K, double eps) {
  for (int k = 0; k < K; k++) {
    if (std::isnan(p[k])) {
      p[k] = 0.5;
      continue;
    }
    p[k] = std::min(std::max(p[k], eps), 1.0 - eps);
  }
}

}  // namespace openadmixture

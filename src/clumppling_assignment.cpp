// Cluster alignment for the native clumppling port.
//
// clumppling (Liu, Kopelman and Rosenberg 2024) aligns the clusters of two
// membership matrices P (N x K_P) and Q (N x K_Q, K_Q >= K_P) by solving an
// integer linear program (align_ILP() in clumppling/core.py, solved with
// GLPK through cvxpy): every column of Q is mapped to exactly one column of
// P, every column of P receives at least one column of Q, and the sum of the
// squared Euclidean distances (divided by N) between matched columns is
// minimised.
//
// We solve the same problem exactly without an ILP solver. When K_Q == K_P it
// is a linear assignment problem. When K_Q > K_P, we add K_Q - K_P "free"
// columns to the cost matrix, where the cost of assigning column q of Q is
// the smallest distance between q and any column of P: the K_P original
// columns guarantee that every column of P is covered, and a column of Q
// assigned to a free column is mapped to its closest column of P. Any
// solution of the ILP maps to an assignment with the same or lower cost, and
// vice versa, so the optima coincide. The square assignment problem is solved
// with the shortest augmenting path (Hungarian) algorithm in O(K^3).

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

namespace {

struct clumppling_alignment {
  double obj;
  std::vector<int> idx; // 0-based column of P for each column of Q
};

// Solve the square linear sum assignment problem (minimisation) for the
// n x n row-major cost matrix `a`; returns the column assigned to each row.
std::vector<int> clumppling_lsap(const std::vector<double>& a, int n) {
  const double inf = std::numeric_limits<double>::infinity();
  // 1-based potentials and matchings, with index 0 as a sentinel
  std::vector<double> u(n + 1, 0.0), v(n + 1, 0.0), minv(n + 1);
  std::vector<int> p(n + 1, 0), way(n + 1, 0);
  std::vector<char> used(n + 1);
  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::fill(minv.begin(), minv.end(), inf);
    std::fill(used.begin(), used.end(), 0);
    do {
      used[j0] = 1;
      const int i0 = p[j0];
      double delta = inf;
      int j1 = 0;
      for (int j = 1; j <= n; ++j) {
        if (!used[j]) {
          const double cur = a[(i0 - 1) * n + (j - 1)] - u[i0] - v[j];
          if (cur < minv[j]) {
            minv[j] = cur;
            way[j] = j0;
          }
          if (minv[j] < delta) {
            delta = minv[j];
            j1 = j;
          }
        }
      }
      if (j1 == 0) {
        Rcpp::stop("the assignment problem could not be solved (non-finite costs?)");
      }
      for (int j = 0; j <= n; ++j) {
        if (used[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        } else {
          minv[j] -= delta;
        }
      }
      j0 = j1;
    } while (p[j0] != 0);
    do {
      const int j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0 != 0);
  }
  std::vector<int> row2col(n);
  for (int j = 1; j <= n; ++j) {
    row2col[p[j] - 1] = j - 1;
  }
  return row2col;
}

// Squared Euclidean distance between column q of Q and column p of P,
// divided by the number of rows (as cdist(Q.T, P.T, 'sqeuclidean') / N).
// Matrices are column-major (as in R). Returns a k_q x k_p row-major matrix.
std::vector<double> clumppling_sqdist(const double* P, int n, int kp,
                                      const double* Q, int kq) {
  std::vector<double> d(static_cast<size_t>(kq) * kp);
  for (int q = 0; q < kq; ++q) {
    const double* qc = Q + static_cast<size_t>(q) * n;
    for (int p = 0; p < kp; ++p) {
      const double* pc = P + static_cast<size_t>(p) * n;
      double s = 0.0;
      for (int i = 0; i < n; ++i) {
        const double diff = qc[i] - pc[i];
        s += diff * diff;
      }
      d[static_cast<size_t>(q) * kp + p] = s / n;
    }
  }
  return d;
}

// Optimal alignment given the k_q x k_p distance matrix `d` (see the
// reduction described at the top of the file).
clumppling_alignment clumppling_align_from_dist(const std::vector<double>& d,
                                                int kp, int kq) {
  std::vector<double> a(static_cast<size_t>(kq) * kq);
  std::vector<int> closest(kq);
  for (int q = 0; q < kq; ++q) {
    const double* dq = &d[static_cast<size_t>(q) * kp];
    int best = 0;
    for (int p = 1; p < kp; ++p) {
      if (dq[p] < dq[best]) best = p;
    }
    closest[q] = best;
    for (int c = 0; c < kq; ++c) {
      a[static_cast<size_t>(q) * kq + c] = c < kp ? dq[c] : dq[best];
    }
  }
  const std::vector<int> row2col = clumppling_lsap(a, kq);
  clumppling_alignment res;
  res.idx.resize(kq);
  res.obj = 0.0;
  for (int q = 0; q < kq; ++q) {
    res.idx[q] = row2col[q] < kp ? row2col[q] : closest[q];
    res.obj += d[static_cast<size_t>(q) * kp + res.idx[q]];
  }
  return res;
}

clumppling_alignment clumppling_align_ilp(const double* P, int n, int kp,
                                          const double* Q, int kq) {
  return clumppling_align_from_dist(clumppling_sqdist(P, n, kp, Q, kq), kp, kq);
}

// As align_ILP_diff1() in clumppling: when Q has exactly one more cluster
// than P, try merging every pair of columns of Q (in the order given by
// itertools.combinations), align the merged matrix to P, and keep the first
// merge with the lowest objective.
clumppling_alignment clumppling_align_ilp_diff1(const double* P, int n, int kp,
                                                const double* Q, int kq) {
  const std::vector<double> d = clumppling_sqdist(P, n, kp, Q, kq);
  clumppling_alignment best;
  best.obj = std::numeric_limits<double>::infinity();
  std::vector<int> old_idx;
  old_idx.reserve(kq);
  std::vector<double> merged(n);
  std::vector<double> d_comb(static_cast<size_t>(kp) * kp);
  for (int a = 0; a < kq - 1; ++a) {
    for (int b = a + 1; b < kq; ++b) {
      old_idx.clear();
      for (int q = 0; q < kq; ++q) {
        if (q != a && q != b) old_idx.push_back(q);
      }
      // the unmerged columns keep their distances, the merged column (the
      // last one, as in clumppling) needs new ones
      for (size_t k = 0; k < old_idx.size(); ++k) {
        for (int p = 0; p < kp; ++p) {
          d_comb[k * kp + p] = d[static_cast<size_t>(old_idx[k]) * kp + p];
        }
      }
      const double* qa = Q + static_cast<size_t>(a) * n;
      const double* qb = Q + static_cast<size_t>(b) * n;
      for (int i = 0; i < n; ++i) merged[i] = qa[i] + qb[i];
      const std::vector<double> d_merged =
          clumppling_sqdist(P, n, kp, merged.data(), 1);
      for (int p = 0; p < kp; ++p) {
        d_comb[static_cast<size_t>(kp - 1) * kp + p] = d_merged[p];
      }
      const clumppling_alignment res =
          clumppling_align_from_dist(d_comb, kp, kp);
      if (res.obj < best.obj) {
        best.obj = res.obj;
        best.idx.assign(kq, 0);
        for (size_t k = 0; k < old_idx.size(); ++k) {
          best.idx[old_idx[k]] = res.idx[k];
        }
        best.idx[a] = res.idx[kp - 1];
        best.idx[b] = res.idx[kp - 1];
      }
    }
  }
  return best;
}

void clumppling_check_pair(const Rcpp::NumericMatrix& P,
                           const Rcpp::NumericMatrix& Q) {
  if (P.nrow() != Q.nrow()) {
    Rcpp::stop("P and Q must have the same number of rows (individuals)");
  }
  if (P.ncol() > Q.ncol()) {
    Rcpp::stop("Q must have no fewer columns (clusters) than P");
  }
}

Rcpp::List clumppling_alignment_to_r(const clumppling_alignment& res) {
  Rcpp::IntegerVector idx(res.idx.size());
  for (size_t q = 0; q < res.idx.size(); ++q) idx[q] = res.idx[q] + 1;
  return Rcpp::List::create(Rcpp::Named("obj") = res.obj,
                            Rcpp::Named("idx") = idx);
}

} // namespace

// Solve a square linear assignment problem (minimisation); returns the
// 1-based column assigned to each row. Used for testing.
// [[Rcpp::export]]
Rcpp::IntegerVector clumppling_lsap_cpp(Rcpp::NumericMatrix cost) {
  const int n = cost.nrow();
  if (cost.ncol() != n) Rcpp::stop("the cost matrix must be square");
  std::vector<double> a(static_cast<size_t>(n) * n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) a[static_cast<size_t>(i) * n + j] = cost(i, j);
  }
  const std::vector<int> row2col = clumppling_lsap(a, n);
  Rcpp::IntegerVector res(n);
  for (int i = 0; i < n; ++i) res[i] = row2col[i] + 1;
  return res;
}

// Align the columns of Q to those of P (align_ILP() in clumppling). Returns
// the optimal objective and, for each column of Q, the 1-based column of P it
// maps to.
// [[Rcpp::export]]
Rcpp::List clumppling_align_ilp_cpp(Rcpp::NumericMatrix P,
                                    Rcpp::NumericMatrix Q) {
  clumppling_check_pair(P, Q);
  return clumppling_alignment_to_r(clumppling_align_ilp(
      P.begin(), P.nrow(), P.ncol(), Q.begin(), Q.ncol()));
}

// Align Q to P when Q has exactly one more column than P, merging two of its
// columns (align_ILP_diff1() in clumppling).
// [[Rcpp::export]]
Rcpp::List clumppling_align_ilp_diff1_cpp(Rcpp::NumericMatrix P,
                                          Rcpp::NumericMatrix Q) {
  clumppling_check_pair(P, Q);
  if (Q.ncol() - P.ncol() != 1) {
    Rcpp::stop("Q must have exactly one more column (cluster) than P");
  }
  return clumppling_alignment_to_r(clumppling_align_ilp_diff1(
      P.begin(), P.nrow(), P.ncol(), Q.begin(), Q.ncol()));
}

// Align all pairs of replicates with the same K (align_within_k() in
// clumppling). For each pair i < j, Q_j is aligned to Q_i, and the cost is
// ||Q_i - aligned Q_j||_F^2 / (2N). Returns the symmetric cost matrix and a
// K x R x R array where [, i, j] maps the columns of Q_j onto those of Q_i
// (1-based), with the reverse alignment for i > j and the identity for i = j.
// [[Rcpp::export]]
Rcpp::List clumppling_align_within_k_cpp(Rcpp::List q_list) {
  const int r = q_list.size();
  std::vector<Rcpp::NumericMatrix> qs;
  qs.reserve(r);
  for (int i = 0; i < r; ++i) qs.push_back(Rcpp::as<Rcpp::NumericMatrix>(q_list[i]));
  const int n = r > 0 ? qs[0].nrow() : 0;
  const int k = r > 0 ? qs[0].ncol() : 0;
  for (int i = 1; i < r; ++i) {
    if (qs[i].nrow() != n || qs[i].ncol() != k) {
      Rcpp::stop("all Q matrices must have the same dimensions");
    }
  }
  Rcpp::NumericMatrix cost(r, r);
  Rcpp::IntegerVector pattern(static_cast<size_t>(k) * r * r);
  pattern.attr("dim") = Rcpp::IntegerVector::create(k, r, r);
  auto pat = [&](int c, int i, int j) -> int& {
    return pattern[c + static_cast<size_t>(k) * (i + static_cast<size_t>(r) * j)];
  };
  for (int i = 0; i < r; ++i) {
    for (int c = 0; c < k; ++c) pat(c, i, i) = c + 1;
  }
  for (int i = 0; i < r - 1; ++i) {
    const double* P = qs[i].begin();
    for (int j = i + 1; j < r; ++j) {
      const double* Q = qs[j].begin();
      const clumppling_alignment res = clumppling_align_ilp(P, n, k, Q, k);
      // cost of P against Q with its columns moved to the matched columns of
      // P (alignQ_wrtP() with merge, then cost_membership())
      double s = 0.0;
      for (int c = 0; c < k; ++c) {
        const double* pc = P + static_cast<size_t>(res.idx[c]) * n;
        const double* qc = Q + static_cast<size_t>(c) * n;
        for (int ii = 0; ii < n; ++ii) {
          const double diff = pc[ii] - qc[ii];
          s += diff * diff;
        }
      }
      // numpy computes the Frobenius norm and then squares it
      const double norm = std::sqrt(s);
      cost(i, j) = cost(j, i) = norm * norm / (2.0 * n);
      for (int c = 0; c < k; ++c) {
        pat(c, i, j) = res.idx[c] + 1;
        pat(res.idx[c], j, i) = c + 1;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("cost") = cost,
                            Rcpp::Named("pattern") = pattern);
}

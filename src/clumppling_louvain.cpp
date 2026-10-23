// Louvain community detection for the native clumppling port.
//
// clumppling detects modes with cdlib's louvain(), which calls
// best_partition() from python-louvain (Thomas Aynaud, BSD-3-Clause,
// https://github.com/taynaud/python-louvain). This is a C++ port of that
// implementation: the same local moving phase (nodes and candidate
// communities are visited in random order, a node only moves for a strictly
// positive gain, a level stops when modularity improves by less than 1e-7),
// the same aggregation into an induced graph, and the partition at the top
// level of the dendrogram is returned. Random orders are drawn from R's RNG
// rather than numpy's, so results for a given seed are not identical to the
// python ones when the community structure is ambiguous.
//
// Graphs are dense weight matrices (the number of nodes, i.e. replicates, is
// small); a zero weight means no edge, as in networkx.from_numpy_array(), and
// the diagonal holds self loops.

#include <Rcpp.h>
#include <vector>

namespace {

const double clumppling_louvain_min = 0.0000001;

struct clumppling_graph {
  int n;
  std::vector<double> w; // n x n symmetric weights, diagonal = self loops
  std::vector<double> degree; // weighted degree, self loops counted twice
  double total_weight; // sum of the weights of all edges (self loops once)

  double weight(int i, int j) const { return w[static_cast<size_t>(i) * n + j]; }

  void finalise() {
    degree.assign(n, 0.0);
    total_weight = 0.0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        degree[i] += weight(i, j);
        if (j >= i) total_weight += weight(i, j);
      }
      degree[i] += weight(i, i);
    }
  }
};

struct clumppling_status {
  std::vector<int> node2com;
  std::vector<double> degrees; // total degree of each community
  std::vector<double> gdegrees; // degree of each node
  std::vector<double> loops; // self loop of each node
  std::vector<double> internals; // internal weight of each community
  double total_weight;

  void init(const clumppling_graph& g) {
    node2com.resize(g.n);
    degrees = g.degree;
    gdegrees = g.degree;
    loops.resize(g.n);
    internals.resize(g.n);
    for (int i = 0; i < g.n; ++i) {
      node2com[i] = i;
      loops[i] = g.weight(i, i);
      internals[i] = loops[i];
    }
    total_weight = g.total_weight;
  }
};

// Fisher-Yates shuffle (as numpy's shuffle), with R's RNG
void clumppling_shuffle(std::vector<int>& x) {
  for (int i = static_cast<int>(x.size()) - 1; i > 0; --i) {
    const int j = static_cast<int>(R_unif_index(i + 1));
    std::swap(x[i], x[j]);
  }
}

double clumppling_modularity(const clumppling_status& s, double resolution) {
  const double links = s.total_weight;
  std::vector<char> present(s.node2com.size(), 0);
  for (int c : s.node2com) present[c] = 1;
  double result = 0.0;
  for (size_t c = 0; c < present.size(); ++c) {
    if (present[c] && links > 0) {
      const double frac = s.degrees[c] / (2.0 * links);
      result += s.internals[c] * resolution / links - frac * frac;
    }
  }
  return result;
}

// One level of local moving (__one_level() in python-louvain)
void clumppling_one_level(const clumppling_graph& g, clumppling_status& s,
                          double resolution) {
  const int n = g.n;
  std::vector<int> nodes(n);
  // weight towards each neighbouring community, and the communities in the
  // order in which they are found
  std::vector<double> neigh_weight(n, 0.0);
  std::vector<char> neigh_seen(n, 0);
  std::vector<int> neigh_coms;
  neigh_coms.reserve(n);

  bool modified = true;
  double new_mod = clumppling_modularity(s, resolution);
  while (modified) {
    const double cur_mod = new_mod;
    modified = false;
    for (int i = 0; i < n; ++i) nodes[i] = i;
    clumppling_shuffle(nodes);
    for (int node : nodes) {
      const int com_node = s.node2com[node];
      const double degc_totw = s.gdegrees[node] / (s.total_weight * 2.0);
      neigh_coms.clear();
      for (int nb = 0; nb < n; ++nb) {
        const double w = g.weight(node, nb);
        if (nb == node || w == 0.0) continue;
        const int c = s.node2com[nb];
        if (!neigh_seen[c]) {
          neigh_seen[c] = 1;
          neigh_weight[c] = 0.0;
          neigh_coms.push_back(c);
        }
        neigh_weight[c] += w;
      }
      const double w_own = neigh_seen[com_node] ? neigh_weight[com_node] : 0.0;
      const double remove_cost =
          -w_own + resolution * (s.degrees[com_node] - s.gdegrees[node]) * degc_totw;
      // remove the node from its community
      s.degrees[com_node] -= s.gdegrees[node];
      s.internals[com_node] -= w_own + s.loops[node];
      s.node2com[node] = -1;

      int best_com = com_node;
      double best_increase = 0.0;
      clumppling_shuffle(neigh_coms);
      for (int com : neigh_coms) {
        const double incr = remove_cost + neigh_weight[com] -
                            resolution * s.degrees[com] * degc_totw;
        if (incr > best_increase) {
          best_increase = incr;
          best_com = com;
        }
      }
      // insert it into the best community
      const double w_best = neigh_seen[best_com] ? neigh_weight[best_com] : 0.0;
      s.node2com[node] = best_com;
      s.degrees[best_com] += s.gdegrees[node];
      s.internals[best_com] += w_best + s.loops[node];
      if (best_com != com_node) modified = true;
      for (int c : neigh_coms) neigh_seen[c] = 0;
    }
    new_mod = clumppling_modularity(s, resolution);
    if (new_mod - cur_mod < clumppling_louvain_min) break;
  }
}

// Relabel communities as 0, 1, ... in order of first appearance; returns the
// number of communities
int clumppling_renumber(std::vector<int>& node2com) {
  std::vector<int> map(node2com.size(), -1);
  int m = 0;
  for (int& c : node2com) {
    if (map[c] < 0) map[c] = m++;
    c = map[c];
  }
  return m;
}

// Graph whose nodes are the communities of `partition` (induced_graph() in
// python-louvain): edge weights are summed, and edges within a community
// become self loops
clumppling_graph clumppling_induced_graph(const std::vector<int>& partition,
                                          int m, const clumppling_graph& g) {
  clumppling_graph ret;
  ret.n = m;
  ret.w.assign(static_cast<size_t>(m) * m, 0.0);
  for (int i = 0; i < g.n; ++i) {
    for (int j = i; j < g.n; ++j) {
      const double w = g.weight(i, j);
      if (w == 0.0) continue;
      const int ci = partition[i];
      const int cj = partition[j];
      ret.w[static_cast<size_t>(ci) * m + cj] += w;
      if (ci != cj) ret.w[static_cast<size_t>(cj) * m + ci] += w;
    }
  }
  ret.finalise();
  return ret;
}

} // namespace

// Louvain partition of the weighted graph with adjacency matrix `adj` (zero
// means no edge; the diagonal is ignored, as clumppling removes self loops).
// Returns 0-based community labels in order of first appearance.
// [[Rcpp::export]]
Rcpp::IntegerVector clumppling_louvain_cpp(Rcpp::NumericMatrix adj,
                                           double resolution = 1.0) {
  const int n = adj.nrow();
  if (adj.ncol() != n) Rcpp::stop("the adjacency matrix must be square");
  clumppling_graph g;
  g.n = n;
  g.w.assign(static_cast<size_t>(n) * n, 0.0);
  bool has_edges = false;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      const double w = adj(i, j);
      if (!(w >= 0.0)) Rcpp::stop("edge weights must be non-negative and finite");
      g.w[static_cast<size_t>(i) * n + j] = w;
      if (w != 0.0) has_edges = true;
    }
  }
  Rcpp::IntegerVector result(n);
  if (!has_edges) {
    for (int i = 0; i < n; ++i) result[i] = i;
    return result;
  }
  g.finalise();

  // the partition of the original nodes at the current top level
  std::vector<int> membership(n);
  for (int i = 0; i < n; ++i) membership[i] = i;

  clumppling_status status;
  status.init(g);
  clumppling_one_level(g, status, resolution);
  double mod = clumppling_modularity(status, resolution);
  std::vector<int> partition = status.node2com;
  int m = clumppling_renumber(partition);
  for (int i = 0; i < n; ++i) membership[i] = partition[membership[i]];
  clumppling_graph current = clumppling_induced_graph(partition, m, g);
  status.init(current);
  while (true) {
    clumppling_one_level(current, status, resolution);
    const double new_mod = clumppling_modularity(status, resolution);
    if (new_mod - mod < clumppling_louvain_min) break;
    partition = status.node2com;
    m = clumppling_renumber(partition);
    for (int i = 0; i < n; ++i) membership[i] = partition[membership[i]];
    mod = new_mod;
    current = clumppling_induced_graph(partition, m, current);
    status.init(current);
  }
  clumppling_renumber(membership);
  for (int i = 0; i < n; ++i) result[i] = membership[i];
  return result;
}

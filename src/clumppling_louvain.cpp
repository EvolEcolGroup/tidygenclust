// Louvain community detection for the native clumppling port.
//
// clumppling detects modes with cdlib's louvain(), which calls
// best_partition() from python-louvain (Thomas Aynaud, BSD-3-Clause,
// https://github.com/taynaud/python-louvain; see inst/COPYRIGHTS). This is a
// C++ port of that implementation: the same local moving phase (nodes and
// candidate communities are visited in random order, a node only moves for a
// strictly positive gain, a level stops when modularity improves by less than
// 1e-7),
// the same aggregation into an induced graph, and the partition at the top
// level of the dendrogram is returned. Random orders are drawn from R's RNG
// rather than numpy's, so results for a given seed are not identical to the
// python ones when the community structure is ambiguous.
//
// Graphs are dense weight matrices (the number of nodes, i.e. replicates, is
// small); a zero weight means no edge, as in networkx.from_numpy_array(), and
// the diagonal holds self loops.
//
// The algorithm (Blondel et al. 2008) alternates two phases:
// 1. local moving (clumppling_one_level()): starting with every node in its
//    own community, nodes are visited in random order and moved to the
//    neighbouring community with the largest gain in modularity, until no
//    node moves or a full pass improves modularity by less than 1e-7;
// 2. aggregation (clumppling_induced_graph()): every community becomes a
//    node of a new graph, the weights of the edges between two communities
//    are summed, and the edges within a community become a self loop.
// The two phases are repeated on the aggregated graph until a level improves
// modularity by less than 1e-7. Each level is a partition of the nodes of
// the previous one; their composition (`membership`) is the partition of
// the original nodes.
//
// Index conventions: at every level, nodes are 0, ..., n - 1, and a
// community is identified by the index of a node (each node starts in the
// community with its own index), so arrays indexed by community have n
// elements. Communities are renumbered 0, ..., m - 1 before aggregation, so
// that community c becomes node c of the aggregated graph.
//
// Modularity is computed as in python-louvain, sum over communities c of
// in_c * resolution / m - (tot_c / 2m)^2, where in_c is the weight of the
// edges within c (self loops included, each edge counted once), tot_c the
// total degree of the nodes in c, and m the total weight of the graph. Note
// that python-louvain applies the resolution to the first term here, but to
// the second one in the gains of local moves (as in the usual definition);
// this is kept for fidelity, and only matters if resolution != 1.

#include <Rcpp.h>
#include <cmath>
#include <vector>

namespace {

// the smallest improvement in modularity that continues the algorithm
// (__MIN in python-louvain)
const double clumppling_louvain_min = 0.0000001;

// An undirected weighted graph, stored as a dense matrix
struct clumppling_graph {
  int n;
  std::vector<double> w; // n x n symmetric weights, diagonal = self loops
  std::vector<double> degree; // weighted degree, self loops counted twice
  double total_weight; // sum of the weights of all edges (self loops once)

  double weight(int i, int j) const { return w[static_cast<size_t>(i) * n + j]; }

  // compute the degrees and total weight once the weights are set, as
  // networkx's degree(weight=) and size(weight=)
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

// The state of local moving (Status in python-louvain)
struct clumppling_status {
  std::vector<int> node2com; // community of each node (-1 while it moves)
  std::vector<double> degrees; // total degree of each community
  std::vector<double> gdegrees; // degree of each node
  std::vector<double> loops; // self loop of each node
  std::vector<double> internals; // internal weight of each community
  double total_weight;

  // every node in its own community (Status.init() without a partition)
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

// Fisher-Yates shuffle (as numpy's shuffle), with R's RNG: it draws
// size - 1 random indices, like numpy, but not the same ones. R_unif_index()
// needs the RNG state, which the RNGScope added by Rcpp::export provides.
void clumppling_shuffle(std::vector<int>& x) {
  for (int i = static_cast<int>(x.size()) - 1; i > 0; --i) {
    const int j = static_cast<int>(R_unif_index(i + 1));
    std::swap(x[i], x[j]);
  }
}

// Modularity of the current partition (__modularity() in python-louvain, see
// the top of the file), summed over the communities that have nodes
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
  // passes over all the nodes, until no node moves or a pass improves
  // modularity by less than clumppling_louvain_min
  while (modified) {
    const double cur_mod = new_mod;
    modified = false;
    // a new random order of the nodes at each pass; python shuffles the
    // nodes in the order of the graph, which is also 0, ..., n - 1
    for (int i = 0; i < n; ++i) nodes[i] = i;
    clumppling_shuffle(nodes);
    for (int node : nodes) {
      const int com_node = s.node2com[node];
      // k_i / 2m: the degree of the node over twice the total weight
      const double degc_totw = s.gdegrees[node] / (s.total_weight * 2.0);
      // weights from the node to each neighbouring community (__neighcom()),
      // accumulated in neigh_weight (a dense replacement for python's dict,
      // reset through neigh_seen after each node); self loops are excluded
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
      // change in modularity (times m) when the node leaves its community:
      // it loses the weight w_own of its edges to the other members, and the
      // null-model term resolution * (tot_own - k_i) * k_i / 2m
      const double w_own = neigh_seen[com_node] ? neigh_weight[com_node] : 0.0;
      const double remove_cost =
          -w_own + resolution * (s.degrees[com_node] - s.gdegrees[node]) * degc_totw;
      // remove the node from its community
      s.degrees[com_node] -= s.gdegrees[node];
      s.internals[com_node] -= w_own + s.loops[node];
      s.node2com[node] = -1;

      // net change in modularity (times m) for moving the node to each
      // neighbouring community com: remove_cost + w_com - resolution *
      // tot_com * k_i / 2m. The node stays in its community unless a move
      // gives a strictly positive gain; on ties the first community in the
      // (shuffled) order wins. Returning to its own community gives a gain
      // of zero up to rounding.
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
      // insert it into the best community (the weight of its edges to the
      // community becomes internal)
      const double w_best = neigh_seen[best_com] ? neigh_weight[best_com] : 0.0;
      s.node2com[node] = best_com;
      s.degrees[best_com] += s.gdegrees[node];
      s.internals[best_com] += w_best + s.loops[node];
      if (best_com != com_node) modified = true;
      for (int c : neigh_coms) neigh_seen[c] = 0;
    }
    new_mod = clumppling_modularity(s, resolution);
    // (written so that a NaN gain also stops the loop)
    if (!(new_mod - cur_mod >= clumppling_louvain_min)) break;
  }
}

// Relabel communities as 0, 1, ... in order of first appearance; returns the
// number of communities. python-louvain's __renumber() chooses different
// labels, but adds the communities to the induced graph in the same order
// (that of their first node), so the nodes of the next level are visited in
// the same order before shuffling; the labels themselves do not matter.
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
  // every edge (i, j), i <= j, once, self loops included, as graph.edges();
  // an edge within a community is added once to the self loop of the
  // community, an edge between communities to both symmetric entries
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
  if (!std::isfinite(resolution) || resolution <= 0.0)
    Rcpp::stop("the resolution must be a positive number");
  clumppling_graph g;
  g.n = n;
  g.w.assign(static_cast<size_t>(n) * n, 0.0);
  bool has_edges = false;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      const double w = adj(i, j);
      if (!(w >= 0.0) || !std::isfinite(w))
        Rcpp::stop("edge weights must be non-negative and finite");
      g.w[static_cast<size_t>(i) * n + j] = w;
      if (w != 0.0) has_edges = true;
    }
  }
  Rcpp::IntegerVector result(n);
  // as in python-louvain, a graph without edges gives singletons
  if (!has_edges) {
    for (int i = 0; i < n; ++i) result[i] = i;
    return result;
  }
  g.finalise();

  // the partition of the original nodes at the current top level
  std::vector<int> membership(n);
  for (int i = 0; i < n; ++i) membership[i] = i;

  // generate_dendrogram() in python-louvain: the first level is always kept
  // (even without any improvement), later levels only if they improve
  // modularity by at least clumppling_louvain_min. Rather than storing the
  // dendrogram and composing its levels at the end (partition_at_level()),
  // each kept level is composed into `membership` as it is found.
  clumppling_status status;
  status.init(g);
  clumppling_one_level(g, status, resolution);
  double mod = clumppling_modularity(status, resolution);
  std::vector<int> partition = status.node2com;
  int m = clumppling_renumber(partition);
  // original node i was in node membership[i] of this level's graph, which
  // is now in community partition[membership[i]]
  for (int i = 0; i < n; ++i) membership[i] = partition[membership[i]];
  clumppling_graph current = clumppling_induced_graph(partition, m, g);
  status.init(current);
  while (true) {
    clumppling_one_level(current, status, resolution);
    const double new_mod = clumppling_modularity(status, resolution);
    // the last level did not improve: its partition is discarded
    // (written so that a NaN gain also stops the loop)
    if (!(new_mod - mod >= clumppling_louvain_min)) break;
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

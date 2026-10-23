# Community detection for the native clumppling port, used to group
# replicates with the same K into modes (cd_default() and detect_communities()
# in clumppling/core.py, cost_to_adj() and test_comm_struc() in
# clumppling/detectMode/utils.py).

# community detection methods available in clumppling
clumppling_cd_methods <- c(
  "louvain", "leiden", "infomap", "markov_clustering", "label_propagation",
  "walktrap", "custom"
)

#' Detect communities of replicates from their cost matrix
#'
#' As `detect_communities()` in clumppling: one or two replicates, and
#' matrices where all costs are below `min_threshold` (or all above
#' `max_threshold`), are handled directly; otherwise community detection is
#' run on the adjacency matrix derived from the costs. Communities are
#' labelled by decreasing size (ties broken by order of first appearance).
#' @param cost_mat a symmetric matrix of alignment costs
#' @param test_comm boolean, whether to test for community structure first
#' @param method the community detection method (see `clumppling_cd_default()`)
#' @param res the resolution for the louvain method
#' @param min_threshold,max_threshold cost thresholds
#' @returns an integer vector of 0-based community labels
#' @keywords internal
#' @noRd
clumppling_detect_communities <- function(cost_mat, test_comm = FALSE,
                                          method = "louvain", res = 1,
                                          min_threshold = 1e-6,
                                          max_threshold = 1e-2) {
  n_nodes <- nrow(cost_mat)
  if (n_nodes == 0) {
    stop("Cost matrix is empty.")
  }
  if (n_nodes == 1) {
    return(0L)
  }
  off_diag <- cost_mat[row(cost_mat) != col(cost_mat)]
  cost_min <- min(off_diag)
  cost_max <- max(off_diag)
  if (n_nodes == 2) {
    return(if (cost_mat[1, 2] < min_threshold) c(0L, 0L) else c(0L, 1L))
  }
  if (cost_max < min_threshold) {
    return(rep(0L, n_nodes))
  }
  if (cost_min > max_threshold) {
    return(seq_len(n_nodes) - 1L)
  }
  if (cost_max == cost_min) {
    # the normalised adjacency matrix is undefined (0/0); clumppling's
    # louvain then finds no gain in modularity, and returns singletons
    return(seq_len(n_nodes) - 1L)
  }

  adj_mat <- clumppling_cost_to_adj(cost_mat, norm = TRUE)
  has_comm_struc <- if (test_comm) clumppling_test_comm_struc(adj_mat) else TRUE
  if (!has_comm_struc) {
    return(rep(0L, n_nodes))
  }
  communities <- clumppling_cd_default(adj_mat, method = method, res = res)
  if (length(communities) == 0) {
    return(seq_len(n_nodes) - 1L)
  }
  clumppling_relabel_communities(communities)
}

#' Relabel communities by decreasing size
#'
#' As done at the end of `detect_communities()` in clumppling: the largest
#' community becomes 0, the next 1, etc., with ties broken by the position
#' of the first member of each community.
#' @param communities a vector of community labels
#' @returns an integer vector of 0-based labels
#' @keywords internal
#' @noRd
clumppling_relabel_communities <- function(communities) {
  # unique() returns labels in order of first appearance
  labels <- unique(communities)
  sizes <- tabulate(match(communities, labels), nbins = length(labels))
  labels <- labels[order(-sizes, method = "radix")]
  match(communities, labels) - 1L
}

#' Convert a cost matrix into an adjacency matrix
#'
#' As `cost_to_adj()` in clumppling: adjacency is 1 - cost, where costs are
#' first rescaled to \[0, 1\] using the off-diagonal range if `norm` is TRUE.
#' The diagonal is set to 1.
#' @param cost_mat a symmetric cost matrix
#' @param norm boolean, whether to rescale the costs
#' @returns the adjacency matrix
#' @keywords internal
#' @noRd
clumppling_cost_to_adj <- function(cost_mat, norm = TRUE) {
  if (norm) {
    off_diag <- cost_mat[row(cost_mat) != col(cost_mat)]
    cost_min <- min(off_diag)
    cost_max <- max(off_diag)
    adj_mat <- 1 - (cost_mat - cost_min) / (cost_max - cost_min)
  } else {
    adj_mat <- 1 - cost_mat
  }
  diag(adj_mat) <- 1
  adj_mat
}

#' Test whether an adjacency matrix has community structure
#'
#' As `test_comm_struc()` in clumppling (with `alpha = 0.01`, the value used
#' by clumppling): the extreme eigenvalues of the standardised adjacency
#' matrix, and of its exponentiated version, are compared with the quantiles
#' of the Tracy-Widom distribution (beta = 1). The quantiles are those
#' computed by clumppling with the 'TracyWidom' python package on a grid from
#' -10 to 10 with step 0.02.
#' @param w a symmetric adjacency matrix
#' @returns TRUE if there is community structure
#' @keywords internal
#' @noRd
clumppling_test_comm_struc <- function(w) {
  ci_max <- 2.8000000000000007 # first grid value with cdf > 1 - 0.01 / 4
  ci_min <- -4.38 # last grid value with cdf < 0.01 / 4
  n <- nrow(w)
  off_diag <- row(w) != col(w)
  standardize <- function(x) {
    v <- x[off_diag]
    v_mean <- mean(v)
    out <- matrix(0, n, n)
    # numpy's std() divides by n
    out[off_diag] <- (v - v_mean) / sqrt(mean((v - v_mean)^2))
    out
  }
  t_mat <- standardize(w) / sqrt(n)
  w_exp <- matrix(0, n, n)
  w_exp[off_diag] <- exp(w[off_diag] * 0.5)
  te_mat <- standardize(w_exp) / sqrt(n)
  eig_t <- eigen(t_mat, symmetric = TRUE, only.values = TRUE)$values
  eig_te <- eigen(te_mat, symmetric = TRUE, only.values = TRUE)$values
  s <- (max(eig_t) < ci_max) + (min(eig_t) > ci_min) +
    (max(eig_te) < ci_max) + (min(eig_te) > ci_min)
  s != 4
}

#' Run a community detection method on an adjacency matrix
#'
#' As `cd_default()` in clumppling. 'louvain' is a C++ port of the
#' python-louvain implementation used by clumppling (via cdlib), and
#' 'markov_clustering' a port of cdlib's implementation of MCL. 'custom'
#' reproduces the example shipped with clumppling (MCL with a pruning
#' threshold of 0.2). 'leiden', 'infomap', 'label_propagation' and 'walktrap'
#' use the 'igraph' package, with the same choice of weighted or unweighted
#' graph as cdlib (only infomap uses the weights). `method` can also be a
#' function that takes the adjacency matrix and returns a vector of community
#' labels.
#' @param adj_mat the adjacency matrix (with ones on the diagonal)
#' @param method the name of the method, or a function
#' @param res the resolution for louvain
#' @returns a vector of community labels, one per node
#' @keywords internal
#' @noRd
clumppling_cd_default <- function(adj_mat, method = "louvain", res = 1) {
  n_nodes <- nrow(adj_mat)
  if (is.function(method)) {
    communities <- method(adj_mat)
    if (length(communities) != n_nodes || anyNA(communities)) {
      stop(
        "the custom community detection function must return one label ",
        "(not NA) for each replicate"
      )
    }
    return(communities)
  }
  # the graph has no self loops
  graph_adj <- adj_mat
  diag(graph_adj) <- 0
  if (method == "louvain") {
    return(clumppling_louvain_cpp(graph_adj, res))
  }
  if (method == "markov_clustering") {
    return(clumppling_mcl(graph_adj))
  }
  if (method == "custom") {
    return(clumppling_mcl(graph_adj, pruning_threshold = 0.2))
  }
  if (!method %in% clumppling_cd_methods) {
    stop(
      "Unknown method: ", method, ". Supported: ",
      paste(clumppling_cd_methods, collapse = ", ")
    )
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("the 'igraph' package is needed for cd_method = '", method, "'")
  }
  graph <- igraph::graph_from_adjacency_matrix(graph_adj,
    mode = "undirected", weighted = TRUE, diag = FALSE
  )
  # igraph uses the weight attribute whenever weights are not given, so we
  # remove it for the methods that cdlib runs on the unweighted graph
  unweighted <- igraph::delete_edge_attr(graph, "weight")
  comm <- switch(method,
    leiden = igraph::cluster_leiden(unweighted,
      objective_function = "modularity"
    ),
    infomap = igraph::cluster_infomap(graph,
      e.weights = igraph::E(graph)$weight
    ),
    label_propagation = igraph::cluster_label_prop(unweighted),
    walktrap = igraph::cluster_walktrap(unweighted)
  )
  as.integer(igraph::membership(comm)) - 1L
}

#' Markov clustering (MCL)
#'
#' A port of the implementation of MCL used by cdlib (and by the 'custom'
#' example in clumppling), from the 'markov_clustering' python package (Guy
#' Allard, MIT license). Nodes that belong to more than one cluster are
#' assigned to the last one, and nodes that belong to none get the label -1,
#' as in cdlib.
#' @param adj_mat the adjacency matrix of the graph
#' @param expansion,inflation,loop_value,iterations parameters of MCL
#'   (defaults as in cdlib)
#' @param pruning_threshold,pruning_frequency,convergence_check_frequency
#'   parameters of MCL (defaults as in cdlib)
#' @returns an integer vector of 0-based cluster labels
#' @keywords internal
#' @noRd
clumppling_mcl <- function(adj_mat, expansion = 2, inflation = 2,
                           loop_value = 1, iterations = 100,
                           pruning_threshold = 0.001, pruning_frequency = 1,
                           convergence_check_frequency = 1) {
  # sklearn's normalize(norm = "l1", axis = 0) leaves zero columns unchanged
  normalize <- function(m) {
    norms <- colSums(abs(m))
    norms[norms == 0] <- 1
    sweep(m, 2, norms, "/")
  }
  m <- adj_mat
  if (loop_value > 0) {
    diag(m) <- loop_value
  }
  m <- normalize(m)
  for (i in seq_len(iterations) - 1) {
    last_m <- m
    expanded <- m
    for (e in seq_len(expansion - 1)) {
      expanded <- expanded %*% m
    }
    m <- normalize(expanded^inflation)
    prune <- pruning_threshold > 0 &&
      i %% pruning_frequency == pruning_frequency - 1
    if (prune) {
      pruned <- m
      pruned[pruned < pruning_threshold] <- 0
      # keep the largest entry of each column
      col_max <- cbind(max.col(t(m), ties.method = "first"), seq_len(ncol(m)))
      pruned[col_max] <- m[col_max]
      m <- pruned
    }
    check <- i %% convergence_check_frequency ==
      convergence_check_frequency - 1
    if (check && all(abs(m - last_m) - 1e-5 * abs(last_m) <= 1e-8)) {
      break
    }
  }
  # clusters are the non-zero entries of the rows of the attractors
  attractors <- which(diag(m) != 0)
  clusters <- unique(lapply(attractors, function(a) which(m[a, ] != 0)))
  # sort clusters as python sorts tuples
  max_len <- max(lengths(clusters))
  padded <- do.call(rbind, lapply(clusters, function(x) {
    c(x, rep(0L, max_len - length(x)))
  }))
  clusters <- clusters[do.call(order, as.data.frame(padded))]
  communities <- rep(-1L, nrow(adj_mat))
  for (i in seq_along(clusters)) {
    communities[clusters[[i]]] <- i - 1L
  }
  communities
}

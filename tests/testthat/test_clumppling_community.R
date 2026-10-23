# Tests of the community detection (R/clumppling_community.R and
# src/clumppling_louvain.cpp)

# modularity of a partition of a weighted graph (self loops ignored)
modularity <- function(w, membership) {
  diag(w) <- 0
  m2 <- sum(w)
  k <- rowSums(w)
  same <- outer(membership, membership, "==")
  sum((w - outer(k, k) / m2)[same]) / m2
}

test_that("communities are relabelled by decreasing size", {
  relabel <- clumppling_relabel_communities
  expect_identical(relabel(c(5, 3, 3, 5, 3, 9)), c(1L, 0L, 0L, 1L, 0L, 2L))
  # ties are broken by order of first appearance
  expect_identical(relabel(c(2, 1, 1, 2)), c(0L, 1L, 1L, 0L))
  expect_identical(relabel(c(-1, 0, 0)), c(1L, 0L, 0L))
})

test_that("costs are converted to adjacencies as in clumppling", {
  cost <- matrix(c(0, 0.1, 0.3, 0.1, 0, 0.2, 0.3, 0.2, 0), 3)
  expect_equal(
    clumppling_cost_to_adj(cost, norm = TRUE),
    matrix(c(1, 1, 0, 1, 1, 0.5, 0, 0.5, 1), 3)
  )
  expect_equal(
    clumppling_cost_to_adj(cost, norm = FALSE),
    matrix(c(1, 0.9, 0.7, 0.9, 1, 0.8, 0.7, 0.8, 1), 3)
  )
})

test_that("special cases of mode detection match clumppling", {
  ref <- clumppling_test_py_reference()$detect_cases
  cost_mats <- list(
    matrix(0),
    matrix(c(0, 5e-7, 5e-7, 0), 2),
    matrix(c(0, 0.5, 0.5, 0), 2),
    (1 - diag(4)) * 5e-7,
    (1 - diag(4)) * 0.5,
    (1 - diag(4)) * 0.005
  )
  for (i in seq_along(cost_mats)) {
    expect_identical(
      clumppling_detect_communities(cost_mats[[i]]),
      as.integer(ref[[i]])
    )
  }
  expect_error(clumppling_detect_communities(matrix(0, 0, 0)), "empty")
})

test_that("louvain finds the same communities as python-louvain", {
  ref <- clumppling_test_py_reference()$graphs
  graphs <- clumppling_test_unit_inputs()$graphs
  for (i in seq_along(graphs)) {
    # the python results do not depend on the seed
    for (py_comm in ref[[i]]$louvain[-1]) {
      expect_identical(py_comm, ref[[i]]$louvain[[1]])
    }
    for (seed in 1:5) {
      withr::with_seed(seed, {
        comm <- clumppling_cd_default(graphs[[i]], "louvain")
      })
      expect_identical(
        clumppling_relabel_communities(comm),
        as.integer(ref[[i]]$louvain[[1]])
      )
    }
  }
})

test_that("louvain is a local optimum of modularity, and uses the seed", {
  withr::local_seed(10)
  # a graph without a clear structure
  w <- matrix(stats::runif(400), 20)
  w[lower.tri(w)] <- t(w)[lower.tri(w)]
  diag(w) <- 0
  comm <- withr::with_seed(1, clumppling_louvain_cpp(w))
  expect_identical(withr::with_seed(1, clumppling_louvain_cpp(w)), comm)
  q <- modularity(w, comm)
  expect_gt(q, modularity(w, rep(1, 20)))
  expect_gt(q, modularity(w, 1:20))
  # moving any single node to another community does not improve modularity
  for (node in 1:20) {
    for (com in setdiff(unique(comm), comm[node])) {
      moved <- comm
      moved[node] <- com
      expect_lte(modularity(w, moved), q + 1e-12)
    }
  }
  # the resolution changes the size of communities
  n_low <- length(unique(withr::with_seed(1, clumppling_louvain_cpp(w, 0.5))))
  n_high <- length(unique(withr::with_seed(1, clumppling_louvain_cpp(w, 2))))
  expect_lte(n_low, n_high)
  # graphs without edges give singletons
  expect_identical(clumppling_louvain_cpp(matrix(0, 3, 3)), 0:2)
  expect_error(clumppling_louvain_cpp(-w), "non-negative")
  expect_error(clumppling_louvain_cpp(w[, -1]), "square")
})

test_that("louvain agrees with igraph on modularity", {
  skip_if_not_installed("igraph")
  graphs <- clumppling_test_unit_inputs()$graphs
  for (w in graphs) {
    diag(w) <- 0
    graph <- igraph::graph_from_adjacency_matrix(w,
      mode = "undirected",
      weighted = TRUE
    )
    ig <- withr::with_seed(1, igraph::cluster_louvain(graph))
    comm <- withr::with_seed(1, clumppling_louvain_cpp(w))
    expect_equal(
      modularity(w, comm),
      modularity(w, as.integer(igraph::membership(ig)))
    )
  }
})

test_that("markov clustering matches cdlib's implementation", {
  ref <- clumppling_test_py_reference()$graphs
  graphs <- clumppling_test_unit_inputs()$graphs
  for (i in seq_along(graphs)) {
    expect_identical(
      clumppling_relabel_communities(
        clumppling_cd_default(graphs[[i]], "markov_clustering")
      ),
      as.integer(ref[[i]]$markov_clustering)
    )
    expect_identical(
      clumppling_relabel_communities(
        clumppling_cd_default(graphs[[i]], "custom")
      ),
      as.integer(ref[[i]]$custom)
    )
  }
  # two disconnected cliques
  w <- matrix(0, 6, 6)
  w[1:3, 1:3] <- 1
  w[4:6, 4:6] <- 1
  expect_identical(clumppling_mcl(w), rep(0:1, each = 3))
})

test_that("the test for community structure matches clumppling", {
  ref <- clumppling_test_py_reference()
  graphs <- clumppling_test_unit_inputs()$graphs
  for (i in seq_along(graphs)) {
    expect_identical(
      clumppling_test_comm_struc(graphs[[i]]),
      ref$graphs[[i]]$test_comm
    )
  }
  for (k_ref in ref$capeverde_within_k) {
    expect_identical(
      clumppling_test_comm_struc(clumppling_cost_to_adj(k_ref$cost_mat)),
      k_ref$test_comm
    )
  }
  # a matrix without structure
  withr::local_seed(1)
  w <- matrix(stats::runif(900, 0.45, 0.55), 30)
  w[lower.tri(w)] <- t(w)[lower.tri(w)]
  diag(w) <- 1
  expect_false(clumppling_test_comm_struc(w))
  # without structure, a single mode is returned
  cost <- 0.005 + (1 - w) / 100
  diag(cost) <- 0
  expect_identical(
    clumppling_detect_communities(cost, test_comm = TRUE),
    rep(0L, 30)
  )
})

test_that("capeverde modes match clumppling for all methods", {
  ref <- clumppling_test_py_reference()$capeverde_within_k
  methods <- c(
    "louvain", "markov_clustering", "infomap", "label_propagation",
    "walktrap"
  )
  for (k_ref in ref) {
    for (method in methods) {
      if (!method %in% c("louvain", "markov_clustering")) {
        skip_if_not_installed("igraph")
      }
      for (test_comm in c(FALSE, TRUE)) {
        py_comm <- k_ref$communities[[paste0(
          method, "_",
          if (test_comm) "True" else "False"
        )]]
        comm <- withr::with_seed(1, clumppling_detect_communities(
          k_ref$cost_mat,
          test_comm = test_comm, method = method
        ))
        expect_identical(comm, as.integer(py_comm), info = method)
      }
    }
  }
})

test_that("community detection methods can be chosen or supplied", {
  graphs <- clumppling_test_unit_inputs()$graphs
  w <- graphs[[1]]
  # a custom function
  halves <- function(adj) rep(1:2, length.out = nrow(adj))
  expect_identical(clumppling_cd_default(w, halves), rep(1:2, length.out = 30))
  expect_error(
    clumppling_cd_default(w, function(adj) 1:2),
    "one label"
  )
  expect_error(clumppling_cd_default(w, "kmeans"), "Unknown method: kmeans")
  skip_if_not_installed("igraph")
  # leiden is not available in clumppling's default environment; as in cdlib,
  # it ignores edge weights
  comm <- withr::with_seed(1, clumppling_cd_default(w, "leiden"))
  expect_length(comm, 30)
  expect_false(anyNA(comm))
  w_unweighted <- (w > 0) * 1
  expect_identical(
    withr::with_seed(1, clumppling_cd_default(w_unweighted, "leiden")),
    comm
  )
})

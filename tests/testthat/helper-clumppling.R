# Helpers for the tests of the native clumppling port. They are also sourced
# by data-raw/clumppling/regenerate_python_fixtures.R, so that the python
# reference results are computed on exactly the same inputs.

# path to the example data; when sourced from data-raw, the option
# tidygenclust.pkg_root points to the source package
clumppling_test_extdata <- function(file) {
  pkg_root <- getOption("tidygenclust.pkg_root")
  if (is.null(pkg_root)) {
    system.file("extdata", file, package = "tidygenclust")
  } else {
    file.path(pkg_root, "inst", "extdata", file)
  }
}

# the example gt_admix objects, created by
# data-raw/clumppling/create_example_gt_admix.R: capeverde has 200 replicates
# with K from 2 to 5, and chicken 60 replicates with K of 17, 19 and 21
clumppling_test_capeverde <- function() {
  readRDS(clumppling_test_extdata("capeverde_admix.rds"))
}

clumppling_test_chicken <- function() {
  readRDS(clumppling_test_extdata("chicken_admix.rds"))
}

# a gt_admix object with the given Q matrices (e.g. a subset of the
# replicates of another gt_admix object)
clumppling_test_admix <- function(q_list) {
  structure(
    list(k = vapply(q_list, ncol, integer(1)), Q = q_list),
    class = c("gt_admix", "list")
  )
}

# 10 replicates of capeverde for each K from 2 to 5, with shuffled
# individuals, and two runs with K = 1 (which are ignored)
clumppling_test_shuffled <- function() {
  q_list <- clumppling_test_capeverde()$Q[c(41:50, 91:100, 141:150, 191:200)]
  ind_order <- withr::with_seed(99, sample(399))
  q_list <- lapply(q_list, function(q) q[ind_order, ])
  clumppling_test_admix(
    c(list(matrix(1, 399, 1)), q_list, list(matrix(1, 399, 1)))
  )
}

# a compact summary of a Q matrix, used to compare large matrices with the
# stored python results: column sums, column sums of squares, and weighted
# column sums (which depend on the order of individuals)
clumppling_test_fingerprint <- function(q) {
  q <- as.matrix(q)
  w <- seq_len(nrow(q)) / nrow(q)
  rbind(colSums(q), colSums(q^2), colSums(q * w))
}

# simulate a random Q matrix (rows sum to 1); with `alpha` small,
# memberships are concentrated on a few clusters
clumppling_test_random_q <- function(n, k, alpha = 1) {
  q <- matrix(stats::rgamma(n * k, shape = alpha), n, k)
  q / rowSums(q)
}

# simulate a weighted graph with planted communities of the given sizes
# (high weights within communities, low weights between them, and a few
# missing edges)
clumppling_test_planted_graph <- function(sizes, within = c(0.7, 1),
                                          between = c(0, 0.2),
                                          p_missing = 0.02) {
  n <- sum(sizes)
  membership <- rep(seq_along(sizes), sizes)
  same <- outer(membership, membership, "==")
  w <- matrix(0, n, n)
  w[same] <- stats::runif(sum(same), within[1], within[2])
  w[!same] <- stats::runif(sum(!same), between[1], between[2])
  w[stats::runif(n * n) < p_missing] <- 0
  w[lower.tri(w)] <- t(w)[lower.tri(w)]
  diag(w) <- 1
  w
}

# the inputs of the unit-level comparisons with python, generated with a
# fixed seed
clumppling_test_unit_inputs <- function() {
  withr::with_seed(20261023, {
    # pairs of Q matrices to align, with K_P <= K_Q
    k_pairs <- list(
      c(2, 2), c(3, 3), c(5, 5), c(8, 8), c(2, 3), c(4, 5), c(7, 8),
      c(3, 5), c(4, 7), c(2, 6)
    )
    ilp <- list()
    for (k in k_pairs) {
      # unstructured matrices
      ilp[[length(ilp) + 1]] <- list(
        P = clumppling_test_random_q(40, k[1]),
        Q = clumppling_test_random_q(40, k[2])
      )
      # Q is a noisy version of P, with some clusters split in two, and
      # permuted clusters
      p <- clumppling_test_random_q(60, k[1], alpha = 0.3)
      q <- p
      for (j in sample(k[1], k[2] - k[1], replace = TRUE)) {
        frac <- stats::runif(60, 0.2, 0.8)
        q <- cbind(q, q[, j] * frac)
        q[, j] <- q[, j] * (1 - frac)
      }
      q <- q[, sample(k[2]), drop = FALSE] +
        matrix(stats::runif(60 * k[2], 0, 0.05), 60, k[2])
      ilp[[length(ilp) + 1]] <- list(P = p, Q = q / rowSums(q))
    }
    graphs <- list(
      clumppling_test_planted_graph(c(15, 10, 5)),
      clumppling_test_planted_graph(c(8, 8, 8, 8)),
      clumppling_test_planted_graph(c(12, 4)),
      clumppling_test_planted_graph(c(6, 6), between = c(0, 0.4))
    )
  })
  list(ilp = ilp, graphs = graphs)
}

# compare a gt_clumppling object with a python reference; if `fingerprint` is
# TRUE, the reference stores fingerprints of the aligned modes
expect_same_clumppling <- function(res, ref, fingerprint = FALSE,
                                   tolerance = 1e-10) {
  testthat::expect_s3_class(res, "gt_clumppling")
  testthat::expect_named(res, names(ref))
  testthat::expect_identical(res$K_range, ref$K_range)
  testthat::expect_identical(res$N, ref$N)
  testthat::expect_identical(res$mode_replicates, ref$mode_replicates)
  testthat::expect_identical(names(res$aligned_modes), names(ref$aligned_modes))
  testthat::expect_identical(
    res$mode_stats[c("Mode", "Representative", "Size")],
    ref$mode_stats[c("Mode", "Representative", "Size")]
  )
  testthat::expect_equal(res$mode_stats, ref$mode_stats, tolerance = tolerance)
  testthat::expect_identical(names(res$cost_acrossK), names(ref$cost_acrossK))
  testthat::expect_equal(res$cost_acrossK, ref$cost_acrossK,
    tolerance = tolerance
  )
  aligned <- res$aligned_modes
  if (fingerprint) {
    aligned <- lapply(aligned, clumppling_test_fingerprint)
  }
  testthat::expect_equal(aligned, ref$aligned_modes, tolerance = tolerance)
}

# the stored python reference results
clumppling_test_py_reference <- function() {
  readRDS(testthat::test_path("fixtures", "clumppling_python_reference.rds"))
}

# Tests of the alignment of modes across K (R/clumppling_across_k.R)

test_that("columns are moved, merged or added as in alignQ_wrtP()", {
  p <- matrix(c(0.5, 0.2, 0.5, 0.8), 2)
  q <- matrix(c(0.1, 0.1, 0.4, 0.1, 0.5, 0.8), 2)
  # the first two columns of Q both map to the first column of P
  merged <- clumppling_align_q_wrt_p(p, q, c(1L, 1L, 2L), merge = TRUE)
  expect_equal(merged$q, matrix(c(0.5, 0.2, 0.5, 0.8), 2))
  expect_identical(merged$pattern, c(1L, 1L, 2L))
  # without merging, the closest one stays, and the other goes to column 3
  not_merged <- clumppling_align_q_wrt_p(p, q, c(1L, 1L, 2L), merge = FALSE)
  expect_identical(not_merged$pattern, c(3L, 1L, 2L))
  expect_identical(not_merged$q, q[, c(2, 3, 1)])
  # a permutation
  expect_identical(clumppling_move_columns(q, c(3L, 1L, 2L)), q[, c(2, 3, 1)])
})

test_that("costs between modes are computed as in clumppling", {
  p <- matrix(c(0.5, 0.2, 0.5, 0.8), 2)
  q <- matrix(c(0.4, 0.3, 0.6, 0.7), 2)
  expect_equal(clumppling_cost_membership(p, q), 0.04 / 4)
  expect_error(clumppling_cost_membership(p, q[, 1, drop = FALSE]), "same size")
  q3 <- cbind(q[, 1], q[, 2] / 2, q[, 2] / 2)
  # only the first column of P is matched to a single column of Q
  expect_equal(
    clumppling_cost_membership_sep(p, q3, c(1L, 2L, 2L)),
    sum((p[, 1] - q3[, 1])^2) / 4
  )
  expect_equal(clumppling_cost_to_similarity(0.04), 0.8)
  expect_identical(clumppling_pattern_to_str(c(2L, 1L, 3L)), "2 1 3")
})

test_that("modes are aligned across K and reordered consistently", {
  withr::local_seed(1)
  # K = 2 and K = 3 modes, where the third cluster splits the second one
  q2 <- clumppling_test_random_q(40, 2, alpha = 0.5)
  # the first part of the split cluster is the larger one
  split <- stats::runif(40, 0.6, 0.9)
  q3 <- cbind(q2[, 2] * split, q2[, 1], q2[, 2] * (1 - split))
  q3_alt <- clumppling_test_random_q(40, 3)
  q_modes <- list(list(q2), list(q3, q3_alt))
  mode_names <- list("K2M1", c("K3M1", "K3M2"))
  res <- clumppling_align_across_k(2:3, q_modes, mode_names)
  expect_identical(names(res$cost_acrossK), c(
    "K2M1-K3M1", "K2M1-K3M2", "K3M1-K3M1", "K3M2-K3M1", "K3M1-K3M2",
    "K3M2-K3M2", "K2M1-K2M1"
  ))
  expect_identical(names(res$alignment_acrossK), names(res$cost_acrossK))
  expect_equal(res$cost_acrossK[["K2M1-K3M1"]], 0)
  expect_identical(res$alignment_acrossK[["K2M1-K3M1"]], c(2L, 1L, 2L))
  expect_identical(res$alignment_acrossK[["K3M1-K3M1"]], 1:3)
  expect_identical(res$best_acrossK[["Best Pair"]], "K2M1-K3M1")
  expect_identical(res$major_acrossK[["Major Pair"]], "K2M1-K3M1")
  expect_named(res$best_acrossK, c(
    "Best Pair", "Cost", "Performance", "Separate-Cluster Cost",
    "Separate-Cluster Performance"
  ))
  # the first cluster of K = 3 matches the first cluster of K = 2
  expect_equal(res$best_acrossK[["Separate-Cluster Cost"]], 0)

  aligned <- clumppling_reorder_q_across_k(2:3, q_modes, mode_names,
    res$alignment_acrossK,
    anchor_pairs = res$best_acrossK[["Best Pair"]]
  )
  expect_named(aligned, c("K2M1", "K3M1", "K3M2"))
  expect_identical(aligned$K2M1, q2)
  # the cluster that is not split comes first, then the split cluster, with
  # the new cluster last
  expect_identical(aligned$K3M1, q3[, c(2, 1, 3)])
  # the other mode is a permutation of the original one
  by_sum <- function(q) q[, order(colSums(q))]
  expect_identical(by_sum(aligned$K3M2), by_sum(q3_alt))
})

test_that("clusters are added at the end when K is not consecutive", {
  withr::local_seed(2)
  q2 <- clumppling_test_random_q(40, 2, alpha = 0.5)
  q4 <- cbind(q2[, 2] * 0.5, q2[, 1] * 0.3, q2[, 2] * 0.5, q2[, 1] * 0.7)
  q_modes <- list(list(q2), list(q4))
  mode_names <- list("K2M1", "K4M1")
  res <- clumppling_align_across_k(c(2L, 4L), q_modes, mode_names)
  expect_identical(res$alignment_acrossK[["K2M1-K4M1"]], c(2L, 1L, 2L, 1L))
  aligned <- clumppling_reorder_q_across_k(c(2L, 4L), q_modes, mode_names,
    res$alignment_acrossK,
    anchor_pairs = res$best_acrossK[["Best Pair"]]
  )
  # the columns closest to the K = 2 clusters keep their place
  expect_identical(aligned$K4M1, q4[, c(4, 1, 2, 3)])
})

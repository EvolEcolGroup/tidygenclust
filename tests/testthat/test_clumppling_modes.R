# Tests of the extraction of modes (R/clumppling_modes.R)

test_that("community labels are converted into modes", {
  expect_identical(
    clumppling_labels_to_modes(c(0L, 1L, 0L, 2L)),
    list(c(1L, 3L), 2L, 4L)
  )
  expect_error(
    clumppling_labels_to_modes(c(0L, 2L)),
    "zero-indexed and consecutive"
  )
})

test_that("the representative replicate has the smallest total cost", {
  cost <- matrix(c(
    0, 1, 2, 9,
    1, 0, 1, 9,
    2, 1, 0, 9,
    9, 9, 9, 0
  ), 4)
  expect_identical(
    clumppling_find_repr_modes(list(1:3, 4L), cost),
    c(2L, 4L)
  )
  # ties go to the first replicate
  expect_identical(clumppling_find_repr_modes(list(c(2L, 3L)), cost), 2L)
})

test_that("modes are summarised as in clumppling", {
  withr::local_seed(1)
  base <- clumppling_test_random_q(30, 3, alpha = 0.5)
  noisy <- function(q) {
    q <- q + matrix(stats::runif(length(q), 0, 0.02), nrow(q))
    q / rowSums(q)
  }
  q_list <- list(
    noisy(base), noisy(base)[, c(2, 1, 3)], clumppling_test_random_q(30, 3),
    noisy(base)[, c(3, 2, 1)]
  )
  q_names <- paste0(1:4, "_K3R", 1:4)
  aligned <- clumppling_align_within_k(q_list)
  modes <- list(c(1L, 2L, 4L), 3L)
  res <- clumppling_extract_modes(modes, aligned$cost, q_names,
    q_list, aligned$pattern,
    label_prefix = "K3M"
  )
  repr <- clumppling_find_repr_modes(modes, aligned$cost)
  expect_identical(res$repr_modes, repr)
  expect_identical(res$mode_stats$Mode, c("K3M1", "K3M2"))
  expect_identical(res$mode_stats$Representative, q_names[repr])
  expect_identical(res$mode_stats$Size, c(3, 1))
  sub_cost <- aligned$cost[modes[[1]], modes[[1]]]
  expect_equal(res$mode_stats$Cost, c(mean(sub_cost[lower.tri(sub_cost)]), 0))
  expect_equal(
    res$mode_stats$Performance,
    c(mean(1 - sqrt(sub_cost[lower.tri(sub_cost)])), 1)
  )
  # the representative and the average of the aligned replicates
  expect_identical(res$rep_q_modes, q_list[repr])
  aligned_q <- lapply(modes[[1]], function(r) {
    q_list[[r]][, order(aligned$pattern[, repr[1], r])]
  })
  expect_equal(res$avg_q_modes[[1]], Reduce(`+`, aligned_q) / 3)
  expect_identical(res$avg_q_modes[[2]], q_list[[3]])
  # alignment of each replicate to the representative of its mode
  expect_identical(res$align_info$Mode, c("K3M1", "K3M1", "K3M1", "K3M2"))
  expect_identical(res$align_info$Replicate, q_names[c(1, 2, 4, 3)])
  expect_identical(res$align_info$Alignment[4], "1 2 3")
  expect_identical(
    res$align_info$Alignment[match(repr[1], modes[[1]])],
    "1 2 3"
  )
})

test_that("alignment within K and modes of capeverde match clumppling", {
  ref <- clumppling_test_py_reference()$capeverde_within_k
  q_data <- clumppling_prepare_q_list(clumppling_test_capeverde()$Q)
  for (k in q_data$k_range) {
    k_ref <- ref[[paste0("K", k)]]
    idx <- which(q_data$k == k)
    expect_identical(q_data$q_names[idx], unlist(k_ref$q_names))
    aligned <- clumppling_align_within_k(q_data$q_mats[idx])
    expect_equal(aligned$cost, k_ref$cost_mat, tolerance = 1e-12)
    # numpy integers are returned to R as doubles
    py_pattern <- array(as.integer(k_ref$pattern), dim(k_ref$pattern))
    expect_identical(aligned$pattern, py_pattern)
  }
})

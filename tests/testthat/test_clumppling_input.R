# Tests of the input handling (R/clumppling_input.R)

test_that("Q matrices are validated and rescaled as in clumppling", {
  q <- matrix(c(0.2, 0.5, 0.8, 0.5), 2)
  expect_true(clumppling_validate_q(q))
  expect_error(clumppling_validate_q(q * 2), "outside \\[0, 1\\]")
  expect_error(clumppling_validate_q(-q), "outside \\[0, 1\\]")
  expect_error(clumppling_validate_q(q * NA), "outside \\[0, 1\\]")
  expect_error(
    clumppling_validate_q(matrix(c(0.2, 0.5, 0.7, 0.5), 2)),
    "Rows 1 do not sum to 1"
  )
  # numpy's allclose() adds a relative tolerance of 1e-5
  shift <- function(x) {
    q[1, 1] <- q[1, 1] + x
    q
  }
  expect_true(clumppling_validate_q(shift(1e-5)))
  expect_error(clumppling_validate_q(shift(1e-4)), "Rows 1 do not sum to 1")
  expect_true(clumppling_validate_q(shift(1e-3), tolerance = 1e-2))
  # rescaling adds 1e-9 to all entries
  q0 <- matrix(c(0, 1, 1, 0), 2)
  expect_equal(clumppling_rescale_q(q0), (q0 + 1e-9) / (1 + 2e-9))
  expect_equal(rowSums(clumppling_rescale_q(q + 1e-7)), c(1, 1))
})

test_that("lists of Q matrices are grouped by K and named", {
  withr::local_seed(1)
  q2a <- clumppling_test_random_q(10, 2)
  q2b <- clumppling_test_random_q(10, 2)
  q3 <- clumppling_test_random_q(10, 3)
  q1 <- matrix(1, 10, 1)
  expect_message(
    res <- clumppling_prepare_q_list(list(q3, q2a, q1, as.data.frame(q2b))),
    "Ignored Q matrix 3 with K=1"
  )
  expect_identical(res$q_names, c("1_K2R1", "2_K2R2", "3_K3R1"))
  expect_identical(res$k_range, 2:3)
  expect_identical(res$k, c(2L, 2L, 3L))
  expect_equal(res$q_mats[[1]], clumppling_rescale_q(q2a))
  expect_equal(res$q_mats[[2]], clumppling_rescale_q(q2b))
  expect_equal(res$q_mats[[3]], clumppling_rescale_q(q3))
  expect_error(
    clumppling_prepare_q_list(list(q2a, q3[-1, ])),
    "Q matrix 2 has 9 individuals, expected 10"
  )
  expect_error(
    clumppling_prepare_q_list(list(q2a, q2b * 2)),
    "Q matrix 2: Matrix contains values outside"
  )
})

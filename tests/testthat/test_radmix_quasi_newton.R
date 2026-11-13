# The quasi-Newton step of Zhou, Alexander and Lange (2011), computed in C++
# from the q and p blocks (src/radmix_quasi_newton.cpp)

# reference: the formula on the stacked vectors x = c(q, p), with the secant
# pairs in the columns of u_mat and v_mat and the step taken from x0 and x1:
# x_qn = x1 - V (U'(U - V))^-1 U'(x0 - x1)
qn_reference <- function(u_mat, v_mat, x0, x1) {
  coef <- solve(
    crossprod(u_mat, u_mat - v_mat),
    crossprod(u_mat, x0 - x1)
  )
  drop(x1 - v_mat %*% coef)
}

# the secant pair of three successive states, as a column of U and of V
secant_pair <- function(s0, s1, s2) {
  list(
    u = c(s1$q, s1$p) - c(s0$q, s0$p),
    v = c(s2$q, s2$p) - c(s1$q, s1$p)
  )
}

test_that("the quasi-Newton step matches the formula", {
  set.seed(1)
  k <- 3
  n_indiv <- 20
  # more parameters than a chunk (65536), to exercise the chunked sums
  n_loci <- 30000
  for (m in 1:3) {
    states <- replicate(m + 2,
      radmix_ref_start(n_indiv, n_loci, k, seed = NULL),
      simplify = FALSE
    )
    pairs <- lapply(seq_len(m), function(i) {
      secant_pair(states[[i]], states[[i + 1]], states[[i + 2]])
    })
    u_mat <- sapply(pairs, `[[`, "u")
    v_mat <- sapply(pairs, `[[`, "v")
    if (m == 1) {
      u_mat <- matrix(u_mat)
      v_mat <- matrix(v_mat)
    }
    s0 <- states[[m]]
    s1 <- states[[m + 1]]
    res <- radmix_qn_step_cpp(u_mat, v_mat, s0$q, s0$p, s1$q, s1$p, 2)
    x_ref <- qn_reference(u_mat, v_mat, c(s0$q, s0$p), c(s1$q, s1$p))
    # split the stacked vector back into q and p, and project each block
    # onto its feasible set, as the C++ step does
    q_ref <- radmix_project_cpp(matrix(x_ref[seq_len(k * n_indiv)], k), TRUE)
    p_ref <- radmix_project_cpp(
      matrix(x_ref[k * n_indiv + seq_len(k * n_loci)], k), FALSE
    )
    expect_equal(res$q, q_ref, tolerance = 1e-8)
    expect_equal(res$p, p_ref, tolerance = 1e-8)
    # the same with any number of threads
    expect_identical(
      radmix_qn_step_cpp(u_mat, v_mat, s0$q, s0$p, s1$q, s1$p, 1),
      res
    )
  }
})

test_that("a singular system skips the step", {
  set.seed(2)
  k <- 2
  st <- radmix_ref_start(5, 10, k, seed = NULL)
  s1 <- list(q = st$q * 0.9, p = st$p * 0.9)
  s2 <- list(q = st$q * 0.8, p = st$p * 0.8)
  pair <- secant_pair(st, s1, s2)
  # identical secant pairs make U'(U - V) singular
  u_mat <- cbind(pair$u, pair$u)
  v_mat <- cbind(pair$v, pair$v)
  expect_null(radmix_qn_step_cpp(u_mat, v_mat, st$q, st$p, s1$q, s1$p, 1))
  # all-zero pairs (e.g. with nothing left to update)
  zero <- matrix(0, k * 15, 1)
  expect_null(radmix_qn_step_cpp(zero, zero, st$q, st$p, st$q, st$p, 1))
})

test_that("invalid arguments are errors", {
  st <- radmix_ref_start(5, 10, 2, seed = NULL)
  u_mat <- matrix(0.1, 2 * 15, 2)
  expect_error(
    radmix_qn_step_cpp(u_mat[-1, ], u_mat, st$q, st$p, st$q, st$p, 1),
    "one row per parameter"
  )
  expect_error(
    radmix_qn_step_cpp(u_mat, u_mat, st$q, st$p, st$q, st$p[, -1], 1),
    "dimensions"
  )
  expect_error(
    radmix_qn_step_cpp(u_mat[, 0], u_mat[, 0], st$q, st$p, st$q, st$p, 1),
    "one row per parameter"
  )
  expect_error(
    radmix_qn_step_cpp(u_mat, u_mat, st$q, st$p, st$q, st$p, 0),
    "n_cores"
  )
})

test_that("the quasi-Newton step checks the number of rows", {
  # q with 2 rows and p with 3: the 30 values of p would be written into a
  # 2 x 10 output
  q <- matrix(0.5, 2, 10)
  p <- matrix(0.5, 3, 10)
  u_mat <- matrix(0.1, 50, 1)
  expect_error(
    radmix_qn_step_cpp(u_mat, u_mat, q, p, q, p, 1),
    "same number of rows"
  )
})

# The quasi-Newton step of Zhou, Alexander and Lange (2011), computed in C++
# from the q and p blocks (src/openadmixture_quasi_newton.cpp)

# reference: the formula on the stacked vectors x = c(q, p), with the secant
# pairs in `pairs` (each a list of x0, x1 and x2) and the step taken from the
# latest x0 and x1
qn_reference <- function(pairs, x0, x1) {
  u_mat <- sapply(pairs, function(x) x$x1 - x$x0)
  v_mat <- sapply(pairs, function(x) x$x2 - x$x1)
  coef <- solve(
    crossprod(u_mat, u_mat - v_mat),
    crossprod(u_mat, x0 - x1)
  )
  drop(x1 - v_mat %*% coef)
}

random_state <- function(k, n_indiv, n_loci) {
  q <- matrix(stats::runif(k * n_indiv), k)
  q <- sweep(q, 2, colSums(q), "/")
  p <- matrix(stats::runif(k * n_loci, 0.05, 0.95), k)
  list(q = q, p = p)
}

test_that("the quasi-Newton step matches the formula", {
  set.seed(1)
  k <- 3
  n_indiv <- 20
  # more parameters than a chunk (65536), to exercise the chunked sums
  n_loci <- 30000
  state <- openadmixture_qn_new_cpp(k * (n_indiv + n_loci), 3)
  xs <- list()
  for (iter in 1:4) {
    s0 <- random_state(k, n_indiv, n_loci)
    s1 <- random_state(k, n_indiv, n_loci)
    s2 <- random_state(k, n_indiv, n_loci)
    col <- (iter - 1) %% 3 + 1
    openadmixture_qn_push_cpp(state, col, s0$q, s0$p, s1$q, s1$p, s2$q, s2$p)
    xs[[col]] <- list(
      x0 = c(s0$q, s0$p), x1 = c(s1$q, s1$p), x2 = c(s2$q, s2$p)
    )
    n_used <- min(iter, 3)
    res <- openadmixture_qn_step_cpp(state, n_used, s0$q, s0$p, s1$q, s1$p, 2)
    # the pairs stored so far (the order of the columns does not matter)
    x_ref <- qn_reference(xs, xs[[col]]$x0, xs[[col]]$x1)
    q_ref <- openadmixture_project_q_cpp(matrix(x_ref[seq_len(k * n_indiv)], k))
    p_ref <- openadmixture_project_p_cpp(
      matrix(x_ref[k * n_indiv + seq_len(k * n_loci)], k)
    )
    expect_equal(res$q, q_ref, tolerance = 1e-8)
    expect_equal(res$p, p_ref, tolerance = 1e-8)
    # the same with any number of threads
    expect_identical(
      openadmixture_qn_step_cpp(state, n_used, s0$q, s0$p, s1$q, s1$p, 1),
      res
    )
  }
})

test_that("a singular system skips the step", {
  set.seed(2)
  k <- 2
  st <- random_state(k, 5, 10)
  state <- openadmixture_qn_new_cpp(k * 15, 3)
  # identical secant pairs make U'(U - V) singular
  for (col in 1:2) {
    openadmixture_qn_push_cpp(
      state, col, st$q, st$p, st$q * 0.9, st$p * 0.9,
      st$q * 0.8, st$p * 0.8
    )
  }
  expect_null(openadmixture_qn_step_cpp(
    state, 2, st$q, st$p, st$q * 0.9,
    st$p * 0.9, 1
  ))
  # all-zero pairs (e.g. with nothing left to update)
  zero_state <- openadmixture_qn_new_cpp(k * 15, 3)
  openadmixture_qn_push_cpp(zero_state, 1, st$q, st$p, st$q, st$p, st$q, st$p)
  expect_null(
    openadmixture_qn_step_cpp(zero_state, 1, st$q, st$p, st$q, st$p, 1)
  )
})

test_that("invalid arguments are errors", {
  st <- random_state(2, 5, 10)
  state <- openadmixture_qn_new_cpp(2 * 15, 3)
  expect_error(
    openadmixture_qn_push_cpp(state, 4, st$q, st$p, st$q, st$p, st$q, st$p),
    "out of range"
  )
  expect_error(
    openadmixture_qn_push_cpp(
      state, 1, st$q, st$p[, -1], st$q, st$p, st$q, st$p
    ),
    "dimensions"
  )
  expect_error(
    openadmixture_qn_step_cpp(state, 0, st$q, st$p, st$q, st$p, 1),
    "out of range"
  )
})

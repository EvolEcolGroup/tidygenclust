# reference solution with quadprog: minimise 0.5 d' H d - grad' d subject to
# 0 <= par + d <= 1 and, optionally, sum(d) == 1 - sum(par)
qp_reference <- function(h, grad, par, simplex) {
  k <- length(grad)
  # solve.QP() imposes t(amat) %*% d >= bvec (the first meq as equalities):
  # here d >= -par and -d >= par - 1, i.e. 0 <= par + d <= 1, preceded by
  # the equality sum(d) == 1 - sum(par) for the simplex
  amat <- cbind(diag(k), -diag(k))
  bvec <- c(-par, par - 1)
  meq <- 0
  if (simplex) {
    amat <- cbind(rep(1, k), amat)
    bvec <- c(1 - sum(par), bvec)
    meq <- 1
  }
  quadprog::solve.QP(h, grad, amat, bvec, meq = meq)$solution
}

# a random problem: a positive definite h, a starting point inside the
# feasible set (on the simplex, with entries of at least 0.05) and a gradient
# whose scale controls how often the bounds become active
random_qp <- function(k, grad_scale, simplex) {
  x <- matrix(rnorm((k + 3) * k), k + 3, k)
  h <- crossprod(x) + diag(0.01, k)
  par <- if (simplex) {
    p <- rexp(k)
    p <- p / sum(p)
    0.05 + p * (1 - 0.05 * k)
  } else {
    runif(k, 0.05, 0.95)
  }
  list(h = h, grad = rnorm(k, sd = grad_scale), par = par)
}

test_that("QP matches quadprog with and without the sum constraint", {
  set.seed(1)
  for (simplex in c(FALSE, TRUE)) {
    n_active <- 0
    for (rep in 1:300) {
      k <- sample(2:8, 1)
      prob <- random_qp(k, grad_scale = sample(c(0.1, 1, 10), 1), simplex)
      res <- openadmixture_qp_solve_cpp(prob$h, prob$grad, prob$par, simplex)
      ref <- qp_reference(prob$h, prob$grad, prob$par, simplex)
      info <- paste("simplex =", simplex, "rep =", rep)
      expect_equal(res$delta, ref, tolerance = 1e-6, info = info)
      est <- prob$par + res$delta
      # as in Julia, a vertex of the simplex (one entry at 1) is reached but
      # the active set then cycles until the iteration limit, returning 0
      vertex <- simplex && max(est) > 1 - 1e-8
      if (!vertex) expect_gt(res$iter, 0)
      expect_true(all(est >= -1e-10 & est <= 1 + 1e-10), info = info)
      if (simplex) expect_equal(sum(est), 1, tolerance = 1e-10, info = info)
      n_active <- n_active + any(est < 1e-8 | est > 1 - 1e-8)
    }
    # make sure the bounds were actually exercised
    expect_gt(n_active, 30)
  }
})

test_that("a simplex vertex solution is found but reported as iter 0", {
  h <- matrix(c(3.25, 1.13, 2.94, 1.13, 6.06, -2.28, 2.94, -2.28, 9.41), 3, 3)
  grad <- c(10, -1.6, -6.6)
  par <- c(0.36, 0.18, 0.46)
  res <- openadmixture_qp_solve_cpp(h, grad, par, simplex = TRUE)
  expect_equal(par + res$delta, c(1, 0, 0))
  expect_identical(res$iter, 0L)
})

test_that("an interior solution is the unconstrained Newton step", {
  h <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  grad <- c(0.01, -0.02)
  res <- openadmixture_qp_solve_cpp(h, grad, c(0.5, 0.5), simplex = FALSE)
  expect_equal(res$iter, 1)
  expect_equal(res$delta, solve(h, grad))
})

test_that("zero or indefinite Hessians fail with a zero step", {
  grad <- c(1, -1, 0.5)
  par <- c(0.2, 0.3, 0.5)
  bad_h <- list(
    matrix(0, 3, 3),
    diag(c(1, -1, 1)),
    matrix(c(1, 2, 0, 2, 1, 0, 0, 0, 1), 3, 3)
  )
  for (h in bad_h) {
    res <- openadmixture_qp_solve_cpp(h, grad, par, simplex = FALSE)
    expect_identical(res$iter, 0L)
    expect_identical(res$delta, rep(0, 3))
  }
  # the mu * 11' shift cannot rescue a zero Hessian either
  res <- openadmixture_qp_solve_cpp(matrix(0, 3, 3), grad, par, simplex = TRUE)
  expect_identical(res$iter, 0L)
  expect_identical(res$delta, rep(0, 3))
})

test_that("mismatched dimensions are rejected", {
  expect_error(
    openadmixture_qp_solve_cpp(diag(2), c(1, 1, 1), c(0, 0, 0), FALSE)
  )
})

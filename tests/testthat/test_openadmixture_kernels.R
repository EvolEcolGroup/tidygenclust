sim <- openadmixture_simulate_admixed(n = 60, m = 200, k = 3, missing = 0.05)

test_that("kernels match the reference on simulated data with missing values", {
  inp <- openadmixture_ref_inputs(openadmixture_test_gt(sim$genotypes))
  expect_true(anyNA(inp$g))
  openadmixture_expect_ref_match(inp, k = 3)
  openadmixture_expect_ref_match(inp, k = 1)
})

test_that("kernels match the reference on unsorted rows and columns", {
  test_gt <- openadmixture_test_gt(sim$genotypes)
  rows <- c(5, 2, 9, 1, 33, 17, 60, 41, 12, 8, 25)
  cols <- c(150, 3, 77, 1, 200, 64, 65, 12, 199, 100, 2, 130)
  # a literal vector, as tidyselect deprecates external vectors
  sub_gt <- tidypopgen::select_loci(
    test_gt[rows, ],
    c(150, 3, 77, 1, 200, 64, 65, 12, 199, 100, 2, 130)
  )
  inp <- openadmixture_ref_inputs(sub_gt)
  expect_equal(inp$g, sim$genotypes[rows, cols], ignore_attr = TRUE)
  openadmixture_expect_ref_match(inp, k = 3)
})

test_that("kernels match the reference with pseudohaploids", {
  hap_genotypes <- sim$genotypes
  # make the first 15 individuals fully homozygous
  hom <- 1:15
  hap_genotypes[hom, ] <- 2 * (hap_genotypes[hom, ] >= 1)
  pseudo_gt <- tidypopgen::gt_pseudohaploid(
    openadmixture_test_gt(hap_genotypes)
  )
  inp <- openadmixture_ref_inputs(pseudo_gt)
  expect_true(any(inp$n == 1))
  expect_true(any(inp$n == 2))
  openadmixture_expect_ref_match(inp, k = 3)
})

test_that("kernels handle tiles at the edges and tiny data", {
  big <- openadmixture_simulate_admixed(n = 130, m = 150, k = 2, seed = 3)
  inp <- openadmixture_ref_inputs(openadmixture_test_gt(big$genotypes))
  openadmixture_expect_ref_match(inp, k = 4)
  tiny <- openadmixture_simulate_admixed(
    n = 3, m = 2, k = 2, missing = 0, seed = 4
  )
  inp <- openadmixture_ref_inputs(openadmixture_test_gt(tiny$genotypes))
  openadmixture_expect_ref_match(inp, k = 2)
})

test_that("loglik handles extreme frequencies and non-integer dosages", {
  # the last two loci only carry the reference and only the alternate
  # allele, so that there are loci with exact frequencies of 0 and 1 below
  # (the simulation does not always produce them)
  genotypes <- sim$genotypes
  genotypes[, 199] <- ifelse(is.na(genotypes[, 199]), NA, 0)
  genotypes[, 200] <- ifelse(is.na(genotypes[, 200]), NA, 2)
  inp <- openadmixture_ref_inputs(openadmixture_test_gt(genotypes))
  st <- openadmixture_ref_start(60, 200, 3)
  p <- st$p
  # tiny qp where some genotypes are 2 (finite but below the batching limit)
  j_tiny <- which(colSums(inp$g == 2, na.rm = TRUE) > 0)[1:3]
  p[, j_tiny] <- 1e-16
  # exact 0 and 1 where they do not make the likelihood infinite: loci where
  # no observed genotype carries the alternate (j_zero) or the reference
  # (j_one) allele
  j_zero <- which(colSums(inp$g > 0, na.rm = TRUE) == 0)
  j_one <- which(colSums(inp$g < 2, na.rm = TRUE) == 0)
  expect_true(199 %in% j_zero && 200 %in% j_one)
  expect_false(any(c(j_zero, j_one) %in% j_tiny))
  p[, j_zero] <- 0
  p[, j_one] <- 1
  ll <- openadmixture_loglik_cpp(inp$geno, st$q, p, 2)
  expect_true(is.finite(ll))
  expect_equal(ll, openadmixture_ref_loglik(inp$g, inp$n, st$q, p),
    tolerance = 1e-10
  )
  # a qp of exactly 0 for an observed alt allele gives -Inf, not NaN
  p[, 1] <- 0
  expect_true(any(inp$g[, 1] > 0, na.rm = TRUE))
  expect_identical(openadmixture_loglik_cpp(inp$geno, st$q, p, 2), -Inf)
  # dosages that are not hard calls use the direct formula
  # (decoding byte 1 as 0.7 instead of 1 turns heterozygotes into
  # non-integer dosages)
  dat <- inp$dat
  dat$lut[1:3] <- c(0, 0.7, 2)
  geno <- dat
  g <- openadmixture_decode(dat)
  expect_true(any(g == 0.7, na.rm = TRUE))
  expect_equal(
    openadmixture_loglik_cpp(geno, st$q, st$p, 2),
    openadmixture_ref_loglik(g, inp$n, st$q, st$p),
    tolerance = 1e-10
  )
})

test_that("gradients and Hessians match numerical derivatives", {
  small <- openadmixture_simulate_admixed(
    n = 7, m = 12, k = 2, missing = 0.1, seed = 5
  )
  hom <- 1:2
  small$genotypes[hom, ] <- 2 * (small$genotypes[hom, ] >= 1)
  inp <- openadmixture_ref_inputs(tidypopgen::gt_pseudohaploid(
    openadmixture_test_gt(small$genotypes)
  ))
  st <- openadmixture_ref_start(7, 12, 3)
  q <- st$q
  p <- st$p
  gq <- openadmixture_grad_q_cpp(inp$geno, q, p, 1)
  gp <- openadmixture_grad_p_cpp(inp$geno, q, p, 1)
  # Like OpenADMIXTURE, the Q block differentiates the log-likelihood written
  # as a homogeneous function of q (1 - qp = sum_k q_k (1 - p_k) on the
  # simplex), which only changes the gradient by a constant across k
  ll_q <- function(x) {
    openadmixture_ref_loglik_homog(inp$g, inp$n, matrix(x, nrow(q)), p)
  }
  ll_p <- function(x) {
    openadmixture_ref_loglik(inp$g, inp$n, q, matrix(x, nrow(p)))
  }
  expect_equal(as.vector(gq$grad), numDeriv::grad(ll_q, as.vector(q)),
    tolerance = 1e-6
  )
  expect_equal(as.vector(gp$grad), numDeriv::grad(ll_p, as.vector(p)),
    tolerance = 1e-6
  )
  # on the simplex, the two forms have the same directional derivatives
  true_grad <- matrix(numDeriv::grad(
    function(x) openadmixture_ref_loglik(inp$g, inp$n, matrix(x, nrow(q)), p),
    as.vector(q)
  ), nrow(q))
  shift <- gq$grad - true_grad
  expect_equal(shift, matrix(shift[1, ], nrow(q), ncol(q), byrow = TRUE),
    tolerance = 1e-6
  )
  # the kernels return the negative Hessian of each column
  ll_qi <- function(x, i) {
    q[, i] <- x
    openadmixture_ref_loglik_homog(inp$g, inp$n, q, p)
  }
  ll_pj <- function(x, j) {
    p[, j] <- x
    openadmixture_ref_loglik(inp$g, inp$n, q, p)
  }
  for (i in c(1, 4)) {
    expect_equal(gq$hess[, , i], -numDeriv::hessian(ll_qi, q[, i], i = i),
      tolerance = 1e-5
    )
  }
  for (j in c(2, 12)) {
    expect_equal(gp$hess[, , j], -numDeriv::hessian(ll_pj, p[, j], j = j),
      tolerance = 1e-5
    )
  }
})

test_that("results do not depend on the number of cores", {
  big <- openadmixture_simulate_admixed(n = 130, m = 150, k = 3, seed = 6)
  inp <- openadmixture_ref_inputs(openadmixture_test_gt(big$genotypes))
  st <- openadmixture_ref_start(130, 150, 3)
  kernels <- list(
    openadmixture_loglik_cpp, openadmixture_em_cpp,
    openadmixture_grad_q_cpp, openadmixture_grad_p_cpp,
    openadmixture_update_q_cpp, openadmixture_update_p_cpp
  )
  # identical, not just equal: the threads must not change the order in
  # which the sums are accumulated
  for (kernel in kernels) {
    res1 <- kernel(inp$geno, st$q, st$p, 1)
    expect_identical(kernel(inp$geno, st$q, st$p, 2), res1)
    expect_identical(kernel(inp$geno, st$q, st$p, 4), res1)
  }
})

test_that("kernels do not modify their inputs", {
  inp <- openadmixture_ref_inputs(openadmixture_test_gt(sim$genotypes))
  st <- openadmixture_ref_start(60, 200, 3)
  # adding 0 makes new copies, so that a kernel writing into its inputs
  # would make them differ from st
  q <- st$q + 0
  p <- st$p + 0
  invisible(openadmixture_em_cpp(inp$geno, q, p, 2))
  invisible(openadmixture_update_q_cpp(inp$geno, q, p, 2))
  invisible(openadmixture_update_p_cpp(inp$geno, q, p, 2))
  expect_identical(q, st$q)
  expect_identical(p, st$p)
})

test_that("dimension mismatches and invalid datasets are errors", {
  inp <- openadmixture_ref_inputs(openadmixture_test_gt(sim$genotypes))
  st <- openadmixture_ref_start(60, 200, 3)
  expect_error(
    openadmixture_loglik_cpp(inp$geno, st$q[, -1], st$p, 1),
    "one column per individual"
  )
  expect_error(
    openadmixture_em_cpp(inp$geno, st$q, st$p[, -1], 1),
    "one column per locus"
  )
  expect_error(
    openadmixture_grad_q_cpp(inp$geno, st$q, st$p[-1, ], 1),
    "same number of rows"
  )
  expect_error(
    openadmixture_grad_p_cpp(inp$geno, st$q, st$p, 0),
    "n_cores"
  )
  expect_error(
    openadmixture_update_q_cpp(inp$geno, st$q[0, ], st$p[0, ], 1),
    "at least one row"
  )
  expect_error(
    openadmixture_update_p_cpp(list(), st$q, st$p, 1),
    "prepared by openadmixture_prepare"
  )
  bad <- inp$dat
  bad$lut <- bad$lut[-1]
  expect_error(
    openadmixture_loglik_cpp(bad, st$q, st$p, 1),
    "'lut' must have length 256"
  )
})

test_that("update_q returns columns on the simplex", {
  skip_if(!openadmixture_has_real_qp(), "QP solver not available yet")
  inp <- openadmixture_ref_inputs(openadmixture_test_gt(sim$genotypes))
  st <- openadmixture_ref_start(60, 200, 3)
  q_new <- openadmixture_update_q_cpp(inp$geno, st$q, st$p, 2)
  expect_equal(dim(q_new), dim(st$q))
  expect_equal(colSums(q_new), rep(1, 60), tolerance = 1e-8)
  expect_true(all(q_new >= 0))
  p_new <- openadmixture_update_p_cpp(inp$geno, st$q, st$p, 2)
  expect_true(all(p_new > 0 & p_new < 1))
})

test_that("update_q followed by update_p increases the log-likelihood", {
  skip_if(!openadmixture_has_real_qp(), "QP solver not available yet")
  inp <- openadmixture_ref_inputs(openadmixture_test_gt(sim$genotypes))
  st <- openadmixture_ref_start(60, 200, 3)
  q <- st$q
  p <- st$p
  # prime with a few EM steps, as in OpenADMIXTURE
  for (it in 1:5) {
    em <- openadmixture_em_cpp(inp$geno, q, p, 2)
    q <- em$q
    p <- em$p
  }
  # with missing data, the EM step (like OpenADMIXTURE's) leaves the columns
  # of q summing to less than 1, and the Q block quadratic model is only
  # valid on the simplex
  q <- openadmixture_project_q_cpp(q)
  ll0 <- openadmixture_loglik_cpp(inp$geno, q, p, 2)
  q <- openadmixture_update_q_cpp(inp$geno, q, p, 2)
  ll1 <- openadmixture_loglik_cpp(inp$geno, q, p, 2)
  p <- openadmixture_update_p_cpp(inp$geno, q, p, 2)
  ll2 <- openadmixture_loglik_cpp(inp$geno, q, p, 2)
  expect_gt(ll1, ll0)
  expect_gt(ll2, ll1)
  expect_equal(ll2, openadmixture_ref_loglik(inp$g, inp$n, q, p),
    tolerance = 1e-10
  )
})

# Deviation from OpenADMIXTURE.jl: non-finite Newton steps are discarded
# (see discard_nonfinite_step() in src/openadmixture_kernels.cpp)
test_that("update kernels keep a column whose Newton step is not finite", {
  geno <- rbind(c(2, 0, 1), c(0, 1, 2), c(1, 2, 0), c(2, 2, 1))
  dat <- openadmixture_prepare(openadmixture_test_gt(geno))
  q <- matrix(c(0.5, 0.5, 0.3, 0.7, 0.6, 0.4, 0.2, 0.8), 2, 4)
  p <- matrix(c(0.4, 0.6, 0.3, 0.2, 0.7, 0.5), 2, 3)
  # individual 1 is entirely from population 1, whose allele frequency at
  # locus 1 is exactly 0, yet it carries two copies: qp = 0 and g > 0, so
  # its gradient is infinite
  q[, 1] <- c(1, 0)
  p[1, 1] <- 0
  q_new <- openadmixture_update_q_cpp(dat, q, p, 1)
  expect_false(anyNA(q_new))
  # the step for individual 1 is discarded: its column is only projected
  expect_equal(q_new[, 1], c(1 - 1e-5, 1e-5))
  # the other individuals are still updated
  expect_false(isTRUE(all.equal(q_new[, -1], q[, -1])))
  # likewise for locus 1 in the P update
  p_new <- openadmixture_update_p_cpp(dat, q, p, 1)
  expect_false(anyNA(p_new))
  expect_equal(p_new[, 1], c(1e-5, 0.6))
})

# With few individuals, the Q-side kernels also split the loci into segments
# (two-dimensional blocking, see accum_q_segmented() in
# src/openadmixture_kernels.cpp)
test_that("kernels are correct when the loci are split into segments", {
  # 70 individuals (2 blocks) and 2000 loci: 7 segments
  expect_identical(openadmixture_q_segments_cpp(70, 2000), 7L)
  # many individuals, or few loci, use a single pass
  expect_identical(openadmixture_q_segments_cpp(5000, 2000), 1L)
  expect_identical(openadmixture_q_segments_cpp(70, 300), 1L)
  seg <- openadmixture_simulate_admixed(
    n = 70, m = 2000, k = 3, missing = 0.03, seed = 8
  )
  inp <- openadmixture_ref_inputs(openadmixture_test_gt(seg$genotypes))
  openadmixture_expect_ref_match(inp, k = 3)
  st <- openadmixture_ref_start(70, 2000, 3)
  # the Q update is the QP step of each individual, as solved from the
  # reference gradient and Hessian, then projected
  ref <- openadmixture_ref_grad_q(inp$g, inp$n, st$q, st$p)
  q_ref <- st$q
  for (i in seq_len(ncol(st$q))) {
    step <- openadmixture_qp_solve_cpp(
      ref$hess[, , i], ref$grad[, i], st$q[, i], TRUE
    )
    q_ref[, i] <- st$q[, i] + step$delta
  }
  q_ref <- openadmixture_project_q_cpp(q_ref)
  expect_equal(
    openadmixture_update_q_cpp(inp$geno, st$q, st$p, 1), q_ref,
    tolerance = 1e-8
  )
  # results do not depend on the number of cores
  seg_kernels <- list(
    openadmixture_em_cpp, openadmixture_grad_q_cpp, openadmixture_update_q_cpp
  )
  for (kernel in seg_kernels) {
    res1 <- kernel(inp$geno, st$q, st$p, 1)
    expect_identical(kernel(inp$geno, st$q, st$p, 3), res1)
    expect_identical(kernel(inp$geno, st$q, st$p, 4), res1)
  }
})

test_that("empty inputs give errors rather than crashes", {
  # these used to divide by zero (killing R) or read out of bounds
  expect_identical(openadmixture_q_segments_cpp(0, 1000), 1L)
  expect_error(
    openadmixture_project_q_cpp(matrix(0, 0, 3)),
    "at least one row"
  )
  gt <- openadmixture_test_gt(sim$genotypes)
  dat <- openadmixture_prepare(gt)
  # a dataset without individuals (as could be built by hand)
  empty <- dat
  empty$rows <- integer(0)
  empty$n_trials <- numeric(0)
  empty$scale <- numeric(0)
  p <- matrix(0.5, 3, dat$n_loci)
  expect_error(
    openadmixture_em_cpp(empty, matrix(0, 3, 0), p, 1),
    "at least one individual"
  )
  # projecting an empty gen_tibble onto a fixed P (allowed with any k)
  expect_error(
    gt_openadmixture(gt[0, ], k = 3, init_p = t(sim$p), fix_p = TRUE),
    "at least one individual and one locus"
  )
})

# Comparison with OpenADMIXTURE.jl
#
# Starting from the same values of Q and P, the C++ kernels of gt_radmix()
# should reproduce the internal functions of OpenADMIXTURE.jl to numerical
# precision: the log-likelihood, one EM step, one round of Newton updates of
# Q and P, and the full quasi-Newton algorithm. The data are simulated, with
# missing genotypes, and include more individuals and loci than the 64 x 64
# tiles used by both implementations.
#
# OpenADMIXTURE.jl is no longer actively developed, so rather than running
# Julia every time (which would make the tests depend on it), we compare with
# its results stored in fixtures/openadmixture_julia_reference.rds, created by
# data-raw/radmix/regenerate_julia_fixtures.R. The fixture also holds the
# simulated genotypes and the starting values. P matrices in the fixture have
# already been converted to frequencies of allele_alt (SnpArrays, used by
# OpenADMIXTURE.jl, counts the other allele).
#
# See data-raw/radmix/README.md for how the fixtures were created and how
# to regenerate them.

fixture_file <- test_path("fixtures", "openadmixture_julia_reference.rds")
skip_if_not(file.exists(fixture_file), "OpenADMIXTURE.jl fixtures not found")
fixture <- readRDS(fixture_file)

sim_gt <- radmix_test_gt(fixture$genotypes)
dat <- radmix_prepare(sim_gt)
q0 <- fixture$q0
p0 <- fixture$p0
k <- nrow(q0)

test_that("the log-likelihood matches OpenADMIXTURE.jl", {
  expect_equal(
    radmix_loglik_cpp(dat, q0, p0, 2),
    fixture$loglik,
    tolerance = 1e-12
  )
})

test_that("an EM step matches OpenADMIXTURE.jl", {
  em <- radmix_em_cpp(dat, q0, p0, 2)
  # P is updated exactly as in Julia
  expect_equal(em$p, fixture$em_p, tolerance = 1e-12)
  # Q differs by design when genotypes are missing: Julia divides by 2J, so
  # columns of Q no longer sum to one, whereas we divide by twice the number
  # of observed loci of each individual. The two differ only by the
  # normalisation of each column.
  expect_true(all(colSums(fixture$em_q) < 1 - 1e-3))
  expect_equal(em$q, sweep(fixture$em_q, 2, colSums(fixture$em_q), "/"),
    tolerance = 1e-12
  )
})

test_that("Newton updates of Q and P match OpenADMIXTURE.jl", {
  q1 <- radmix_update_q_cpp(dat, q0, p0, 2)
  expect_equal(q1, fixture$update_q, tolerance = 1e-10)
  # as in Julia, P is updated given the new Q
  p1 <- radmix_update_p_cpp(dat, q1, p0, 2)
  expect_equal(p1, fixture$update_p, tolerance = 1e-10)
})

test_that("the full algorithm follows the same path as OpenADMIXTURE.jl", {
  # (q0 and p0 are K x I and K x J, the starting values are I x K and J x K;
  # with both given, the EM priming is skipped, as in Julia)
  fit <- radmix_fit(dat,
    k = k, init_q = t(q0), init_p = t(p0),
    n_cores = 2
  )
  # same number of iterations, and the same log-likelihood at each of them
  # (loglik_history starts with that of the starting values, which the
  # fixture does not include)
  expect_true(fit$converged)
  expect_identical(fit$n_iter, length(fixture$run_loglik))
  expect_equal(fit$loglik_history[-1], fixture$run_loglik, tolerance = 1e-10)
  # same estimates
  expect_equal(fit$q, t(fixture$run_q), tolerance = 1e-8)
  expect_equal(fit$p, t(fixture$run_p), tolerance = 1e-8)
})

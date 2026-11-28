# Comparison with NGSadmix on simulated genotype likelihoods
#
# With genotype likelihoods (given by another package through
# gt_radmix_data()), gt_radmix() maximises the likelihood of NGSadmix
# (Skotte et al. 2013) with the algorithm of ADMIXTURE. The estimates of
# NGSadmix (from ANGSD 0.940) for simulated likelihoods (60 individuals,
# 600 loci from 3 ancestral populations, reads at a mean depth of 3) are
# stored, for K = 2 and 3 (the best of 10 seeds), in
# fixtures/ngsadmix_gl_reference.rds, with the likelihoods themselves (as
# 16-bit codes, decoded below), so that the test does not need NGSadmix. We
# check that:
# - our log-likelihood of NGSadmix's estimates is the one that NGSadmix
#   reported, once the log(1/3) that NGSadmix counts for each individual
#   without data at a locus is added (we count nothing);
# - gt_radmix() reaches the same maximum, with the same estimates once the
#   clusters are matched;
# - NGSadmix's estimates are a fixed point of our algorithm.
# Our estimates are kept within [1e-5, 1 - 1e-5], as in ADMIXTURE, whereas
# those of NGSadmix go down to 1e-9, so our maximum is slightly lower
# (by 0.02 and 0.08 for K = 2 and 3); within our bounds, it is at least the
# log-likelihood of NGSadmix's estimates.
#
# See data-raw/radmix/README.md for how the fixtures were created and how
# to regenerate them.

fixture_file <- test_path("fixtures", "ngsadmix_gl_reference.rds")
skip_if_not(file.exists(fixture_file), "NGSadmix fixtures not available")
fixture <- readRDS(fixture_file)
# the codes (two rows per individual, u1 then u2, one column per locus) are
# the log-likelihoods of 1 and 2 copies of allele2 relative to 0 copies, in
# units of 1/1024 nat, offset by 32768; both codes 32768 mean no data
u1 <- fixture$codes[c(TRUE, FALSE), ]
u2 <- fixture$codes[c(FALSE, TRUE), ]
gl_ll <- array(0, c(dim(u1), 3))
gl_ll[, , 2] <- (u1 - 32768) / 1024
gl_ll[, , 3] <- (u2 - 32768) / 1024
dat <- radmix_test_gl(gl_ll)
n_no_data <- sum(u1 == 32768 & u2 == 32768)

for (run in fixture$runs) {
  k <- run$k
  # NGSadmix's estimates moved inside our bounds, as gt_radmix() does with
  # starting values, and their log-likelihood (in NGSadmix's convention)
  init_q <- radmix_check_init(run$q, nrow(run$q), k, "init_q", TRUE)
  init_p <- radmix_check_init(run$f, nrow(run$f), k, "init_p", FALSE)
  ll_bounded <- radmix_loglik_cpp(dat, t(init_q), t(init_p), 1) +
    n_no_data * log(1 / 3)

  test_that(paste("our log-likelihood is NGSadmix's, K =", k), {
    expect_gt(n_no_data, 0)
    expect_equal(
      radmix_loglik_cpp(dat, t(run$q), t(run$f), 1) + n_no_data * log(1 / 3),
      run$loglik,
      tolerance = 1e-9
    )
  })

  # the best of 3 random starts, fitted on the prepared data as gt_radmix()
  # does
  fits <- lapply(1:3, function(seed) {
    set.seed(seed)
    radmix_fit(dat, k = k, rtol = 1e-10)
  })
  best <- fits[[which.max(vapply(fits, function(f) f$loglik, numeric(1)))]]

  test_that(paste("gt_radmix() reaches NGSadmix's maximum, K =", k), {
    expect_true(best$converged)
    ll_best <- best$loglik + n_no_data * log(1 / 3)
    expect_lt(abs(ll_best - run$loglik), 0.1)
    # (up to the convergence criterion)
    expect_gt(ll_best, ll_bounded - 1e-4)
    # relabel our clusters to match NGSadmix's, in both Q and P (NGSadmix's
    # F is the frequency of the minor allele, as our P)
    perm <- radmix_match_clusters(best$q, run$q)
    expect_lt(max(abs(best$q[, perm] - run$q)), 0.001)
    expect_lt(max(abs(best$p[, perm] - run$f)), 0.005)
  })

  test_that(paste("NGSadmix's estimates are our fixed point, K =", k), {
    warm <- radmix_fit(dat, k = k, init_q = init_q, init_p = init_p)
    expect_true(warm$converged)
    expect_lte(warm$n_iter, 3)
    expect_gt(warm$loglik + n_no_data * log(1 / 3), ll_bounded - 1e-4)
    expect_lt(max(abs(warm$q - run$q)), 0.001)
  })
}

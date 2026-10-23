# Comparison with ADMIXTURE on the EUR_subset data
#
# ADMIXTURE (Alexander et al. 2009) and OpenADMIXTURE maximise the same
# likelihood with different algorithms and random starts, so we cannot expect
# identical paths, but on well structured data both should reach the same
# maximum. For K = 2 to 4 we compare gt_openadmixture() on the 379
# individuals and 54,051 SNPs of EUR_subset (the data used in the
# OpenADMIXTURE.jl tests) with the estimates of ADMIXTURE, and check that:
# - the log-likelihoods agree;
# - the ancestry proportions (Q) and allele frequencies (P) agree, once the
#   clusters are matched (their labels are arbitrary);
# - the cross-validation errors agree, and rank the values of K in the same
#   order;
# - the estimates of ADMIXTURE are a fixed point of our algorithm: started
#   from them, gt_openadmixture() converges immediately without improving
#   the likelihood, i.e. both programs agree on where the maximum is.
#
# ADMIXTURE is no longer developed (its last version, 1.3.0, dates from 2015)
# and the data do not change, so instead of running ADMIXTURE every time
# (which takes several minutes) we compare with its estimates stored in
# fixtures/openadmixture_admixture_eur.rds, created by
# data-raw/openadmixture/regenerate_admixture_fixtures.R. The data themselves
# are downloaded the first time they are needed (see
# openadmixture_eur_subset_bed() in helper-openadmixture-admixture.R), so the
# tests are skipped without an internet connection.
#
# See data-raw/openadmixture/README.md for how the fixtures were created and how
# to regenerate them.

skip_on_cran()
fixture_file <- test_path("fixtures", "openadmixture_admixture_eur.rds")
skip_if_not(file.exists(fixture_file), "ADMIXTURE fixtures not available")
bed_file <- openadmixture_eur_subset_bed()
skip_if(is.null(bed_file), "EUR_subset could not be downloaded")

fixture <- readRDS(fixture_file)
eur <- tidypopgen::gen_tibble(bed_file,
  backingfile = tempfile("eur_"),
  quiet = TRUE
)
n_cores <- min(4L, parallel::detectCores())

test_that("the data are those used to create the ADMIXTURE fixtures", {
  data_files <- paste0(sub("bed$", "", bed_file), c("bed", "bim", "fam"))
  expect_identical(
    unname(tools::md5sum(data_files)),
    unname(fixture$md5)
  )
  # ADMIXTURE writes P with the loci sorted by chromosome and position; for
  # EUR_subset this is also the order of the input, so the rows of the
  # stored P matrices match the loci of the gen_tibble
  loci <- tidypopgen::show_loci(eur)
  chrom <- as.integer(as.character(loci$chromosome))
  expect_false(is.unsorted(order(chrom, loci$position)))
})

dat <- openadmixture_prepare(eur)
# cross-validation errors of each K, collected in the loop below
oa_cv <- numeric(0)

for (run in fixture$runs) {
  k <- run$k
  q_adm <- run$q_millionths / 1e6
  # ADMIXTURE reports the frequency of the second allele in the bim file
  # (allele_ref), we report the frequency of allele_alt
  p_adm <- 1 - run$p_millionths / 1e6

  test_that(paste("the stored ADMIXTURE estimates are consistent, K =", k), {
    # our log-likelihood of ADMIXTURE's estimates matches the value that
    # ADMIXTURE reported (up to the rounding of the stored estimates); this
    # also confirms how the alleles and loci of the fixtures are oriented
    expect_equal(
      openadmixture_loglik_cpp(dat, t(q_adm), t(p_adm), n_cores),
      run$loglik,
      tolerance = 1e-7
    )
  })

  oa <- gt_openadmixture(eur,
    k = k, seed = 1, n_cores = n_cores,
    crossval = TRUE
  )
  oa_cv <- c(oa_cv, oa$cv)
  q_oa <- unclass(oa$Q[[1]])
  perm <- openadmixture_match_clusters(q_oa, q_adm)
  q_oa <- q_oa[, perm]
  p_oa <- unname(oa$P[[1]])[, perm]

  test_that(paste("log-likelihoods agree with ADMIXTURE, K =", k), {
    expect_true(oa$converged)
    expect_equal(oa$loglik, run$loglik, tolerance = 1e-6)
  })

  test_that(paste("Q and P agree with ADMIXTURE, K =", k), {
    expect_lt(mean(abs(q_oa - q_adm)), 0.005)
    expect_lt(max(abs(q_oa - q_adm)), 0.05)
    expect_gt(stats::cor(as.vector(q_oa), as.vector(q_adm)), 0.999)
    expect_lt(mean(abs(p_oa - p_adm)), 0.005)
    expect_gt(stats::cor(as.vector(p_oa), as.vector(p_adm)), 0.999)
  })

  test_that(paste("cross-validation errors agree with ADMIXTURE, K =", k), {
    # both programs mask 5 random folds of the observed genotypes and report
    # the mean binomial deviance of the masked genotypes, but their folds
    # differ (and ADMIXTURE prints 5 decimal places), so the errors agree
    # closely rather than exactly
    expect_equal(oa$cv, run$cv, tolerance = 1e-3)
  })

  test_that(paste("ADMIXTURE's estimates are our fixed point, K =", k), {
    # the stored estimates are rounded to 6 decimal places, so rows of Q
    # need renormalising before they can be used as starting values
    init_q <- q_adm / rowSums(q_adm)
    warm <- gt_openadmixture(eur,
      k = k, init_q = init_q, init_p = p_adm,
      n_cores = n_cores
    )
    expect_true(warm$converged)
    # converges within a few iterations of starting (the log has one line
    # for the initial log-likelihood, one per iteration, and two at the
    # end)...
    expect_lte(length(warm$log[[1]]) - 3, 7)
    # ...without finding a better maximum than ADMIXTURE
    expect_equal(warm$loglik, run$loglik, tolerance = 1e-7)
    expect_lt(max(abs(unclass(warm$Q[[1]]) - init_q)), 0.01)
  })
}

test_that("cross-validation ranks K as ADMIXTURE does", {
  adm_cv <- vapply(fixture$runs, function(run) run$cv, numeric(1))
  expect_identical(order(oa_cv), order(adm_cv))
})

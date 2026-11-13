# Starts from a fit to a subset of the loci: gt_radmix(init = "subset"),
# with the start itself in R/radmix_subset.R

sub_sim <- radmix_simulate_admixed(
  n = 80, m = 400, k = 3, missing = 0.02, seed = 11
)
sub_gt <- radmix_test_gt(sub_sim$genotypes)

test_that("an EM step from flat frequencies gives ancestry-weighted ones", {
  # radmix_subset_start() extends P to all the loci with this step: from
  # p = 0.5 in every population, the EM step of P gives
  # p_kj = sum_i q_ik g_ij / sum_i n_i q_ik over the observed genotypes
  inp <- radmix_ref_inputs(sub_gt)
  q <- t(sub_sim$q)
  em <- radmix_em_cpp(inp$geno, q, matrix(0.5, 3, 400), 1)
  obs <- !is.na(inp$g)
  g0 <- ifelse(obs, inp$g, 0)
  n_obs <- obs * inp$n
  expect_equal(em$p, (q %*% g0) / (q %*% n_obs), tolerance = 1e-12)
})

test_that("the subset start fits a subset and extends it to all the loci", {
  dat <- radmix_prepare(sub_gt)
  set.seed(3)
  start <- radmix_subset_start(dat, k = 3, frac = 0.25)
  expect_identical(dim(start$q), c(80L, 3L))
  expect_identical(dim(start$p), c(400L, 3L))
  expect_match(start$log[1], "fit of 100 of 400 loci")
  # P is close to the best P given Q: a further Newton step of P changes
  # it little (on average; loci near the bounds converge more slowly)
  p_more <- radmix_update_p_cpp(dat, t(start$q), t(start$p), 1)
  expect_lt(mean(abs(p_more - t(start$p))), 0.01)
  # the start is already close to the fit: Q is near the true ancestry
  perm <- radmix_match_clusters(start$q, sub_sim$q)
  expect_gt(min(diag(stats::cor(start$q[, perm], sub_sim$q))), 0.9)
})

test_that("gt_radmix(init = \"subset\") reaches the fit of a random start", {
  adm_rand <- gt_radmix(sub_gt, k = 3, seed = 1)
  adm_sub <- gt_radmix(sub_gt, k = 3, seed = 1, init = "subset",
    init_subset = 0.25
  )
  expect_true(adm_sub$converged)
  expect_equal(adm_sub$loglik, adm_rand$loglik, tolerance = 1e-6)
  # the log starts with the lines of the subset start, and the fit of all
  # the loci is not primed by EM
  log <- adm_sub$log[[1]]
  expect_match(log[1], "^Subset start: fit of 100 of 400 loci")
  expect_false(any(grepl("EM priming", log)))
  # reproducible with the seed
  expect_identical(
    gt_radmix(sub_gt, k = 3, seed = 1, init = "subset",
      init_subset = 0.25
    )$Q,
    adm_sub$Q
  )
})

test_that("invalid subset starts are errors", {
  for (bad in list(0, 1.5, NA, "0.1", c(0.1, 0.2))) {
    expect_error(
      gt_radmix(sub_gt, k = 2, init = "subset", init_subset = bad),
      "'init_subset'"
    )
  }
  expect_error(
    gt_radmix(sub_gt, k = 3, init = "subset", init_q = sub_sim$q),
    "cannot be combined with 'init_q'"
  )
})

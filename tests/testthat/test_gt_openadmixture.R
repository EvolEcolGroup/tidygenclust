sim <- openadmixture_simulate_admixed(
  n = 80, m = 400, k = 3, missing = 0.02, seed = 42
)
sim_gt <- openadmixture_test_gt(sim$genotypes)

# best match between estimated and true clusters, by correlation of Q columns
# (K = 3: of the 6 permutations, the one that maximises the summed
# correlations between matched columns)
match_clusters <- function(q_est, q_true) {
  perms <- list(
    c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
    c(2, 3, 1), c(3, 1, 2), c(3, 2, 1)
  )
  scores <- vapply(perms, function(x) {
    sum(diag(stats::cor(q_est[, x], q_true)))
  }, numeric(1))
  perms[[which.max(scores)]]
}

test_that("gt_openadmixture returns a valid gt_admix object", {
  adm <- gt_openadmixture(sim_gt, k = 2:3, n_runs = 2, seed = c(1, 2))
  expect_s3_class(adm, "gt_admix")
  expect_identical(adm$k, c(2L, 2L, 3L, 3L))
  expect_length(adm$Q, 4)
  expect_identical(adm$algorithm, "OpenADMIXTURE")
  expect_identical(adm$id, sim_gt$id)
  expect_true(all(adm$converged))
  for (i in seq_along(adm$Q)) {
    q <- adm$Q[[i]]
    expect_s3_class(q, "q_matrix")
    expect_identical(dim(q), c(80L, adm$k[i]))
    expect_equal(unname(rowSums(q)), rep(1, 80), tolerance = 1e-10)
    # estimates are kept within [1e-5, 1 - 1e-5] (with a margin for rounding)
    expect_true(all(q >= 1e-5 - 1e-12))
    p <- adm$P[[i]]
    expect_identical(dim(p), c(400L, adm$k[i]))
    expect_identical(rownames(p), tidypopgen::show_loci(sim_gt)$name)
    expect_true(all(p >= 1e-5 - 1e-12 & p <= 1 - 1e-5 + 1e-12))
  }
  # the log-likelihoods are finite (that they never decrease over the main
  # loop is checked in the next test)
  expect_true(all(is.finite(adm$loglik)))
})

test_that("the log-likelihood is monotone and estimates recover the truth", {
  dat <- openadmixture_prepare(sim_gt)
  set.seed(3)
  fit <- openadmixture_fit(dat, k = 3)
  expect_true(fit$converged)
  # (a small relative tolerance allows for rounding)
  expect_true(all(diff(fit$loglik_history) >= -1e-6 * abs(fit$loglik)))
  perm <- match_clusters(fit$q, sim$q)
  expect_gt(min(diag(stats::cor(fit$q[, perm], sim$q))), 0.9)
  expect_lt(mean(abs(fit$q[, perm] - sim$q)), 0.08)
})

test_that("results are reproducible and the RNG state is restored", {
  # if the state is restored, the first draw after the call is the one that
  # follows set.seed(123) without it
  set.seed(123)
  before <- stats::runif(1)
  set.seed(123)
  adm1 <- gt_openadmixture(sim_gt, k = 3, seed = 99)
  after <- stats::runif(1)
  expect_identical(before, after)
  # the same seed gives the same estimates, with any number of cores
  adm2 <- gt_openadmixture(sim_gt, k = 3, seed = 99, n_cores = 2)
  expect_equal(adm1$Q[[1]], adm2$Q[[1]])
  expect_equal(adm1$loglik, adm2$loglik)
})

test_that("grouped gen_tibbles store the group and work with methods", {
  grouped_gt <- dplyr::group_by(sim_gt, population)
  adm <- gt_openadmixture(grouped_gt, k = 2:3, seed = 1)
  expect_identical(adm$group, sim_gt$population)
  q_mat <- tidypopgen::get_q_matrix(adm, k = 3, run = 1)
  expect_s3_class(q_mat, "q_matrix")
  expect_identical(attr(q_mat, "id"), sim_gt$id)
  tidy_q <- generics::tidy(q_mat, grouped_gt)
  expect_identical(nrow(tidy_q), 80L * 3L)
  expect_s3_class(
    ggplot2::autoplot(adm, k = 3, run = 1, type = "barplot"),
    "ggplot"
  )
  # c() (from tidypopgen) needs the same id, group and algorithm, and sorts
  # the runs by k
  adm_c <- c(adm, gt_openadmixture(grouped_gt, k = 2, seed = 5))
  expect_identical(adm_c$k, c(2L, 2L, 3L))
  expect_no_error(summary(adm_c))
})

test_that("fixing P projects individuals onto known allele frequencies", {
  adm <- gt_openadmixture(sim_gt,
    k = 3, init_p = t(sim$p), fix_p = TRUE,
    seed = 1
  )
  # P is unchanged (up to clamping)
  expect_equal(unname(adm$P[[1]]), pmin(pmax(t(sim$p), 1e-5), 1 - 1e-5))
  q <- unclass(adm$Q[[1]])
  expect_lt(mean(abs(q - sim$q)), 0.05)
})

test_that("fixing Q estimates allele frequencies", {
  adm <- gt_openadmixture(sim_gt,
    k = 3, init_q = sim$q, fix_q = TRUE,
    seed = 1
  )
  expect_equal(unname(unclass(adm$Q[[1]])), sim$q, tolerance = 1e-4)
  expect_gt(stats::cor(as.vector(adm$P[[1]]), as.vector(t(sim$p))), 0.9)
})

test_that("pseudohaploid data can be analysed", {
  hap_genotypes <- sim$genotypes
  # make the first 20 individuals pseudohaploid (homozygous)
  hap_genotypes[1:20, ] <- 2 * (hap_genotypes[1:20, ] >= 1)
  pseudo_gt <- tidypopgen::gt_pseudohaploid(
    openadmixture_test_gt(hap_genotypes)
  )
  expect_identical(tidypopgen::show_ploidy(pseudo_gt), -2)
  adm <- gt_openadmixture(pseudo_gt, k = 3, seed = 1)
  expect_true(adm$converged)
  expect_equal(unname(rowSums(adm$Q[[1]])), rep(1, 80), tolerance = 1e-10)
})

test_that("subsets of individuals and loci are respected", {
  # the same individuals and loci, in the same (unsorted) order, as a subset
  # of a gen_tibble (which keeps the full backing file) and as a new
  # gen_tibble holding only them: with the same seed, the fits must agree
  idx <- c(80:41, 1:20)
  loci_idx <- c(400:201, 1:50)
  sub_gt <- tidypopgen::select_loci(sim_gt[idx, ], tidyselect::all_of(loci_idx))
  sub_mat <- openadmixture_test_gt(sim$genotypes[idx, loci_idx])
  adm_sub <- gt_openadmixture(sub_gt, k = 2, seed = 7)
  adm_mat <- gt_openadmixture(sub_mat, k = 2, seed = 7)
  expect_equal(adm_sub$loglik, adm_mat$loglik)
  expect_equal(unclass(adm_sub$Q[[1]]), unclass(adm_mat$Q[[1]]),
    ignore_attr = TRUE
  )
  expect_identical(adm_sub$id, sim_gt$id[idx])
})

test_that("invalid arguments are caught", {
  expect_error(gt_openadmixture(sim_gt, k = 0), "positive integers")
  expect_error(gt_openadmixture(sim_gt, k = 80), "smaller than")
  expect_error(gt_openadmixture(sim_gt, k = 2, n_runs = 2, seed = 1:3), "seed")
  expect_error(
    gt_openadmixture(sim_gt, k = 2:3, init_p = t(sim$p)),
    "single 'k'"
  )
  expect_error(gt_openadmixture(sim_gt, k = 3, fix_p = TRUE), "requires")
  expect_error(
    gt_openadmixture(sim_gt,
      k = 3, init_q = sim$q, init_p = t(sim$p),
      fix_q = TRUE, fix_p = TRUE
    ),
    "cannot both"
  )
  expect_error(gt_openadmixture(sim_gt, k = 3, init_p = sim$p), "400 rows")
  expect_error(
    gt_openadmixture(sim_gt, k = 3, init_q = sim$q * 2),
    "between 0 and 1|sum to 1"
  )
  expect_error(
    gt_openadmixture(sim_gt, k = 2, use_imputed = TRUE),
    "not been imputed"
  )
})

test_that("verbose progress is reported as messages", {
  small <- openadmixture_test_gt(sim$genotypes[1:30, 1:100])
  expect_message(
    gt_openadmixture(small, k = 2, seed = 1, verbose = TRUE),
    "Initial loglikelihood"
  )
  expect_silent(suppressMessages(
    gt_openadmixture(small, k = 2, seed = 1, verbose = TRUE)
  ))
  # max_iter = 0 returns the starting estimates
  adm0 <- gt_openadmixture(small, k = 2, seed = 1, max_iter = 0)
  expect_false(adm0$converged)
  # the log ends with "Failed to converge after ..." and "Loglikelihood: ..."
  expect_match(utils::tail(adm0$log[[1]], 2)[1], "after 0 iterations")
})

test_that("the log-likelihood is only computed when needed", {
  small <- openadmixture_test_gt(sim$genotypes[1:40, 1:200])
  # without verbose, EM priming is summarised in a single line
  adm <- gt_openadmixture(small, k = 2, seed = 1)
  log_lines <- adm$log[[1]]
  expect_false(any(grepl("^EM iteration", log_lines)))
  expect_true("EM priming: 5 iterations" %in% log_lines)
  # the log-likelihood of the plain sweeps is only computed when the
  # quasi-Newton step is rejected
  iter_lines <- grep("^Iteration", log_lines, value = TRUE)
  expect_true(any(grepl("basic = not computed", iter_lines)))
  # with verbose, the log-likelihood is reported every 5 EM steps, and the
  # estimates are the same
  adm_v <- suppressMessages(
    gt_openadmixture(small, k = 2, seed = 1, verbose = TRUE, n_em_iter = 10)
  )
  em_lines <- grep("^EM iteration", adm_v$log[[1]], value = TRUE)
  expect_identical(
    sub(":.*", "", em_lines),
    c("EM iteration 5", "EM iteration 10")
  )
  # (verbose only adds evaluations of the log-likelihood, which do not
  # change the path of the algorithm)
  adm_10 <- gt_openadmixture(small, k = 2, seed = 1, n_em_iter = 10)
  expect_identical(adm_v$Q, adm_10$Q)
  expect_identical(adm_v$loglik, adm_10$loglik)
})

test_that("a few individuals can be projected onto a fixed P", {
  # with P fixed, each individual is estimated on its own, so k can be as
  # large as the number of individuals (e.g. to project a single new sample
  # onto a reference panel), and the estimates are those obtained by
  # projecting all the individuals together
  p_ref <- t(sim$p)
  k <- nrow(sim$p)
  adm_all <- gt_openadmixture(sim_gt,
    k = k, init_p = p_ref, fix_p = TRUE, seed = 1
  )
  all_q <- unclass(adm_all$Q[[1]])
  for (idx in list(1:2, 5L)) {
    adm <- gt_openadmixture(sim_gt[idx, ],
      k = k, init_p = p_ref, fix_p = TRUE, seed = 1
    )
    expect_equal(unclass(adm$Q[[1]]), all_q[idx, , drop = FALSE],
      tolerance = 1e-8, ignore_attr = TRUE
    )
  }
  # without a fixed P, k must still be smaller than the number of individuals
  expect_error(gt_openadmixture(sim_gt[1:2, ], k = k), "smaller than")
})

test_that("init_p is matched to the loci by its row names", {
  p_ref <- t(sim$p)
  k <- nrow(sim$p)
  loci_names <- tidypopgen::show_loci(sim_gt)$name
  adm <- gt_openadmixture(sim_gt, k = k, init_p = p_ref, fix_p = TRUE, seed = 1)
  # the same P, with the loci in reverse order and named
  p_rev <- p_ref[rev(seq_len(nrow(p_ref))), ]
  rownames(p_rev) <- rev(loci_names)
  adm_rev <- gt_openadmixture(sim_gt,
    k = k, init_p = p_rev, fix_p = TRUE, seed = 1
  )
  expect_identical(adm_rev$Q, adm$Q)
  # a P with other loci as well (e.g. a larger reference panel)
  p_more <- rbind(p_rev, extra = rep(0.5, k))
  adm_more <- gt_openadmixture(sim_gt,
    k = k, init_p = p_more, fix_p = TRUE, seed = 1
  )
  expect_identical(adm_more$Q, adm$Q)
  # a P without some of the loci
  expect_error(
    gt_openadmixture(sim_gt, k = k, init_p = p_rev[-1, ], fix_p = TRUE),
    "do not include all the loci"
  )
})

test_that("cross-validation does not change the estimates, without seed", {
  # cross-validation draws random numbers, which used to shift the random
  # starts of the following runs when no seed was given
  set.seed(10)
  adm <- gt_openadmixture(sim_gt, k = 2:3, n_runs = 2)
  set.seed(10)
  adm_cv <- gt_openadmixture(sim_gt, k = 2:3, n_runs = 2, crossval = TRUE)
  expect_identical(adm_cv$Q, adm$Q)
  expect_identical(adm_cv$loglik, adm$loglik)
  expect_length(adm_cv$cv, 4)
})

test_that("all the arguments are validated", {
  # flags must be TRUE or FALSE (fix_q = NA was taken as TRUE)
  expect_error(
    gt_openadmixture(sim_gt, k = 3, init_q = sim$q, fix_q = NA),
    "'fix_q' should be TRUE or FALSE"
  )
  expect_error(
    gt_openadmixture(sim_gt, k = 2, use_imputed = NA),
    "'use_imputed' should be TRUE or FALSE"
  )
  expect_error(
    gt_openadmixture(sim_gt, k = 2, verbose = "yes"),
    "'verbose' should be TRUE or FALSE"
  )
  # rtol = NA or a negative value disabled convergence
  for (bad in list(NA, -1, Inf, "0.1")) {
    expect_error(gt_openadmixture(sim_gt, k = 2, rtol = bad), "'rtol'")
  }
  # counts are checked in R, with a message naming the argument
  expect_error(
    gt_openadmixture(sim_gt, k = 2, max_iter = -1),
    "'max_iter' should be an integer of at least 0"
  )
  expect_error(gt_openadmixture(sim_gt, k = 2, n_qn = 0), "'n_qn'")
  expect_error(gt_openadmixture(sim_gt, k = 2, n_cores = NA), "'n_cores'")
  expect_error(gt_openadmixture(sim_gt, k = 2, n_em_iter = 1.5), "'n_em_iter'")
  # non-integer values are not truncated
  expect_error(gt_openadmixture(sim_gt, k = 2.5), "'k' should be")
  expect_error(gt_openadmixture(sim_gt, k = 2, n_runs = 1.5), "'n_runs'")
  expect_error(
    gt_openadmixture(sim_gt, k = 2, crossval = TRUE, n_folds = 2.7),
    "'n_folds'"
  )
  # seeds are checked before the first fit
  expect_error(
    gt_openadmixture(sim_gt, k = 2, n_runs = 2, seed = c(1, NA)),
    "'seed' should only contain numbers"
  )
})

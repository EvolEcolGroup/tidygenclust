# Cross-validation, as in ADMIXTURE (see R/radmix_cv.R and the
# cross-validation masks in src/radmix_geno_access.h)

sim <- radmix_simulate_admixed(
  n = 80, m = 400, k = 3, missing = 0.02, seed = 42
)
sim_gt <- radmix_test_gt(sim$genotypes)
dat <- radmix_prepare(sim_gt)

# pure-R reference for the binomial deviance of the genotypes in g (NA for
# missing or hidden), with mu = n * qp; q (K x individuals) and p (K x loci)
# are in the layout of the kernels, and the result is c(sum of the
# deviances, number of genotypes), as returned by radmix_deviance_cpp()
ref_deviance <- function(g, n_trials, q, p) {
  qp <- t(q) %*% p
  n <- matrix(n_trials, nrow(g), ncol(g))
  h <- n - g
  ok <- !is.na(g)
  d <- ifelse(g > 0, g * log(g / (n * qp)), 0) +
    ifelse(h > 0, h * log(h / (n * (1 - qp))), 0)
  c(sum(d[ok]), sum(ok))
}

test_that("folds partition the observed genotypes", {
  n_folds <- 5
  full <- radmix_decode(dat)
  observed <- !is.na(full)
  held_out <- lapply(seq_len(n_folds), function(f) {
    radmix_decode(dat, cv = list(
      mode = "held_out", fold = f, n_folds = n_folds, seed = 7
    ))
  })
  # one row per genotype, one column per fold: is the genotype in that fold?
  in_fold <- sapply(held_out, function(x) as.vector(!is.na(x)))
  # every observed genotype is in exactly one fold, missing ones in none
  # (* 1 turns the logical vector into the doubles returned by rowSums())
  expect_identical(rowSums(in_fold), as.vector(observed) * 1)
  # held-out genotypes keep their values
  for (x in held_out) {
    expect_identical(x[!is.na(x)], full[!is.na(x)])
  }
  # folds are roughly equal (multinomial, 1/5 each)
  sizes <- colSums(in_fold)
  expect_true(all(abs(sizes / sum(observed) - 0.2) < 0.02))
  # the training view of a fold is the complement of its held-out view
  for (f in seq_len(n_folds)) {
    train <- radmix_decode(dat, cv = list(
      mode = "train", fold = f, n_folds = n_folds, seed = 7
    ))
    expect_identical(
      !is.na(train) | !is.na(held_out[[f]]),
      observed
    )
    expect_false(any(!is.na(train) & !is.na(held_out[[f]])))
  }
})

test_that("fold assignment depends on the seed only", {
  mask <- function(seed) {
    !is.na(radmix_decode(dat, cv = list(
      mode = "held_out", fold = 1, n_folds = 5, seed = seed
    )))
  }
  expect_identical(mask(3), mask(3))
  expect_false(identical(mask(3), mask(4)))
})

test_that("invalid masks are rejected", {
  expect_error(
    radmix_decode(dat, cv = list(
      mode = "train", fold = 6, n_folds = 5, seed = 1
    )),
    "cv_fold"
  )
  expect_error(
    radmix_decode(dat, cv = list(
      mode = "train", fold = 1, n_folds = 5, seed = -1
    )),
    "cv_seed"
  )
})

test_that("the deviance kernel matches the R reference", {
  # random q (columns on the simplex) and p, in the layout of the kernels
  set.seed(1)
  k <- 3
  q <- matrix(stats::runif(k * dat$n_indiv), k)
  q <- sweep(q, 2, colSums(q), "/")
  p <- matrix(stats::runif(k * dat$n_loci, 0.01, 0.99), k)
  cv <- list(mode = "held_out", fold = 2, n_folds = 5, seed = 11)
  held_out <- radmix_set_cv(dat, cv = cv)
  g_held_out <- radmix_decode(dat, cv = cv)
  expected <- ref_deviance(g_held_out, dat$n_trials, q, p)
  expect_equal(
    radmix_deviance_cpp(held_out, q, p, 1), expected,
    tolerance = 1e-12
  )
  expect_identical(
    radmix_deviance_cpp(held_out, q, p, 1),
    radmix_deviance_cpp(held_out, q, p, 2)
  )
  # all genotypes (plain view)
  expect_equal(
    radmix_deviance_cpp(dat, q, p, 2),
    ref_deviance(radmix_decode(dat), dat$n_trials, q, p),
    tolerance = 1e-12
  )
})

test_that("the deviance handles pseudohaploids", {
  hap_genotypes <- sim$genotypes
  hap_genotypes[1:20, ] <- 2 * (hap_genotypes[1:20, ] >= 1)
  pseudo_dat <- radmix_prepare(
    tidypopgen::gt_pseudohaploid(radmix_test_gt(hap_genotypes))
  )
  set.seed(2)
  k <- 2
  q <- matrix(stats::runif(k * pseudo_dat$n_indiv), k)
  q <- sweep(q, 2, colSums(q), "/")
  p <- matrix(stats::runif(k * pseudo_dat$n_loci, 0.01, 0.99), k)
  expect_equal(
    radmix_deviance_cpp(pseudo_dat, q, p, 1),
    ref_deviance(radmix_decode(pseudo_dat), pseudo_dat$n_trials, q, p),
    tolerance = 1e-12
  )
})

test_that("gt_radmix computes the cross-validation error", {
  adm <- gt_radmix(sim_gt, k = 2:4, seed = 1, crossval = TRUE)
  expect_length(adm$cv, 3)
  expect_true(all(is.finite(adm$cv) & adm$cv > 0))
  # the data were simulated with K = 3
  expect_identical(adm$k[which.min(adm$cv)], 3L)
  # the log reports each fold and the error, in the format of ADMIXTURE
  expect_length(grep("^CV fold", adm$log[[2]]), 5)
  expect_match(utils::tail(adm$log[[2]], 1), "^CV error \\(K=3\\): ")
  # cross-validation does not change the estimates
  adm_no_cv <- gt_radmix(sim_gt, k = 3, seed = 1)
  expect_equal(adm$Q[[2]], adm_no_cv$Q[[1]])
  expect_equal(adm$loglik[2], adm_no_cv$loglik)
  expect_null(adm_no_cv$cv)
  # reproducible with the same seed, whatever the number of cores
  adm2 <- gt_radmix(sim_gt,
    k = 3, seed = 1, crossval = TRUE,
    n_cores = 2
  )
  expect_identical(adm2$cv, adm$cv[2])
  # the number of folds can be changed
  adm10 <- gt_radmix(sim_gt,
    k = 3, seed = 1, crossval = TRUE,
    n_folds = 10
  )
  expect_length(grep("^CV fold", adm10$log[[1]]), 10)
  # works with tidypopgen's methods (c() sorts the runs by k, so the K = 3
  # run of adm2 comes after those of adm with K = 2 and 3)
  expect_s3_class(ggplot2::autoplot(adm, type = "cv"), "ggplot")
  expect_identical(c(adm, adm2)$cv, c(adm$cv[1:2], adm2$cv, adm$cv[3]))
})

test_that("invalid cross-validation arguments are caught", {
  expect_error(
    gt_radmix(sim_gt, k = 2, crossval = "yes"),
    "'crossval' should be TRUE or FALSE"
  )
  expect_error(
    gt_radmix(sim_gt, k = 2, crossval = TRUE, n_folds = 1),
    "'n_folds'"
  )
  expect_error(
    gt_radmix(sim_gt,
      k = 3, crossval = TRUE, init_p = t(sim$p),
      fix_p = TRUE
    ),
    "cannot be used with"
  )
})

test_that("seeds too large for 64-bit integers are rejected", {
  # through R
  expect_error(
    radmix_decode(dat, cv = list(
      mode = "train", fold = 1, n_folds = 5, seed = 1e20
    )),
    "cv_seed"
  )
  # and in C++, bypassing the check in R (the cast to a 64-bit integer of
  # an infinite or too large value is undefined)
  for (bad_seed in c(Inf, 1e20)) {
    bad <- dat
    bad$cv <- list(mode = 1L, fold = 1L, n_folds = 5L, seed = bad_seed)
    expect_error(radmix_decode_cpp(bad), "cv_seed")
  }
})

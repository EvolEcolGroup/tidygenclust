# Starts from the SVD of the genotypes: gt_radmix(init = "svd"),
# R/radmix_svd.R and src/radmix_als.cpp

# complete genotypes (the SVD needs them, or imputed ones)
svd_sim <- radmix_simulate_admixed(
  n = 60, m = 300, k = 3, missing = 0, seed = 7
)
svd_gt <- radmix_test_gt(svd_sim$genotypes)

test_that("the SVD factors give the allele frequencies at full rank", {
  # with all the components of the centred genotypes (one fewer than the
  # individuals), f 1' + 0.5 Z U' is exactly G' / 2; monomorphic loci,
  # which tidypopgen's default scaling rejects, have zero rows in Z
  g <- svd_sim$genotypes[1:20, 1:50]
  g[, 1] <- 0
  g[, 2] <- 2
  fac <- radmix_svd_factors(radmix_test_gt(g), rank = 19)
  expect_equal(fac$f, colMeans(g) / 2, ignore_attr = TRUE)
  expect_equal(fac$f + 0.5 * fac$z %*% t(fac$u), t(g) / 2,
    ignore_attr = TRUE
  )
  expect_equal(fac$z[1:2, ], matrix(0, 2, 19))
})

test_that("the SVD factors do not depend on the blocks of loci", {
  # bigstatsr scales the loci in blocks of block_size() loci (a single block
  # here, unless the blocks are made smaller): each block must get the
  # centre and scale of its own loci
  fac <- radmix_svd_factors(svd_gt, rank = 5)
  withr::local_options(bigstatsr.block.sizeGB = 50 * 8 * 60 / 1024^3)
  expect_identical(bigstatsr::block_size(60), 50)
  fac_blocks <- radmix_svd_factors(svd_gt, rank = 5)
  expect_equal(fac_blocks$z %*% t(fac_blocks$u), fac$z %*% t(fac$u))
  expect_equal(fac_blocks$f, fac$f)
  # the blocks of the SVD are made smaller during the call only (to save
  # memory), so the user's setting is restored afterwards
  expect_identical(
    getOption("bigstatsr.block.sizeGB"), 50 * 8 * 60 / 1024^3
  )
})

test_that("the ALS matches the R reference, whatever the number of cores", {
  fac <- radmix_svd_factors(svd_gt, rank = 8)
  za <- cbind(0.5 * fac$z, fac$f)
  ua <- cbind(fac$u, 1)
  w <- 1 / pmax(fac$f * (1 - fac$f), 1e-3)
  w <- w / mean(w)
  set.seed(1)
  init_p <- matrix(stats::runif(300 * 3), 300, 3)
  als <- radmix_als_cpp(za, ua, w, init_p,
    max_iter = 1000, tol = 1e-4, n_cores = 1
  )
  ref <- radmix_ref_als(za, ua, w, init_p)
  expect_true(als$converged)
  expect_identical(als$n_iter, ref$n_iter)
  expect_equal(als$q, ref$q, tolerance = 1e-6)
  expect_equal(als$p, ref$p, tolerance = 1e-6)
  # Q on the simplex and P within the bounds of the estimates
  expect_equal(colSums(als$q), rep(1, 60))
  expect_true(all(als$q >= 1e-5 - 1e-12))
  expect_true(all(als$p >= 1e-5 & als$p <= 1 - 1e-5))
  # each individual and locus is solved on its own
  expect_identical(
    radmix_als_cpp(za, ua, w, init_p,
      max_iter = 1000, tol = 1e-4, n_cores = 2
    ),
    als
  )
})

test_that("the ALS of the SVD start runs with one BLAS thread", {
  # BLAS threads compete with the OpenMP threads of the ALS; the number of
  # BLAS threads is restored afterwards
  blas_threads <- RhpcBLASctl::blas_get_num_procs()
  withr::defer(RhpcBLASctl::blas_set_num_threads(blas_threads))
  RhpcBLASctl::blas_set_num_threads(2)
  skip_if(
    RhpcBLASctl::blas_get_num_procs() != 2,
    "the BLAS does not let us set its number of threads"
  )
  als <- radmix_als_cpp
  seen <- integer(0)
  local_mocked_bindings(radmix_als_cpp = function(...) {
    seen <<- c(seen, RhpcBLASctl::blas_get_num_procs())
    als(...)
  })
  fac <- radmix_svd_factors(svd_gt, rank = 8)
  set.seed(1)
  radmix_svd_start(radmix_prepare(svd_gt), factors = fac, k = 3)
  expect_equal(seen, rep(1, 3))
  expect_equal(RhpcBLASctl::blas_get_num_procs(), 2)
})

test_that("fits from the SVD are valid and reach the maximum", {
  adm <- gt_radmix(svd_gt, k = 2:3, n_runs = 2, seed = 1:2, init = "svd")
  expect_true(all(adm$converged))
  for (i in seq_along(adm$Q)) {
    expect_equal(unname(rowSums(adm$Q[[i]])), rep(1, 60), tolerance = 1e-10)
    p <- adm$P[[i]]
    expect_true(all(p >= 1e-5 - 1e-12 & p <= 1 - 1e-5 + 1e-12))
  }
  # the log records the start (k + 5 components), without EM priming
  log <- adm$log[[3]]
  expect_identical(
    log[1], "SVD start: best of 3 ALS solutions, with 8 components"
  )
  expect_length(grep("^ALS solution [1-3]: ", log), 3)
  expect_length(grep("^Starting from ALS solution [1-3]$", log), 1)
  expect_false(any(grepl("EM", unlist(adm$log))))
  # the same maximum as random starts, and the true ancestry
  rnd <- gt_radmix(svd_gt, k = 3, n_runs = 3, seed = 1:3)
  expect_equal(adm$loglik[3], max(rnd$loglik), tolerance = 1e-6)
  q <- unclass(adm$Q[[3]])
  perm <- radmix_match_clusters(q, svd_sim$q)
  expect_gt(min(diag(stats::cor(q[, perm], svd_sim$q))), 0.9)
})

test_that("SVD fits are reproducible and restore the RNG state", {
  set.seed(123)
  adm <- gt_radmix(svd_gt, k = 3, seed = 1, init = "svd")
  after <- stats::runif(1)
  set.seed(123)
  expect_identical(after, stats::runif(1))
  expect_identical(gt_radmix(svd_gt, k = 3, seed = 1, init = "svd"), adm)
})

test_that("the SVD needs imputed genotypes, which the fit ignores", {
  sim_missing <- radmix_simulate_admixed(
    n = 60, m = 300, k = 3, missing = 0.02, seed = 7
  )
  gt_missing <- radmix_test_gt(sim_missing$genotypes)
  expect_error(
    gt_radmix(gt_missing, k = 3, init = "svd"),
    "gt_impute_simple"
  )
  gt_imputed <- tidypopgen::gt_impute_simple(gt_missing, method = "mean0")
  adm <- gt_radmix(gt_imputed, k = 3, seed = 1, init = "svd")
  expect_true(adm$converged)
  # the log-likelihood is that of the observed genotypes
  dat <- radmix_prepare(gt_imputed)
  expect_equal(
    adm$loglik,
    radmix_loglik_cpp(dat, t(unclass(adm$Q[[1]])), t(adm$P[[1]]), 1)
  )
})

test_that("the SVD start works with grouped and pseudohaploid data", {
  grouped <- dplyr::group_by(svd_gt, population)
  adm <- gt_radmix(grouped, k = 3, seed = 1, init = "svd")
  expect_identical(adm$group, svd_gt$population)
  expect_identical(
    adm$Q, gt_radmix(svd_gt, k = 3, seed = 1, init = "svd")$Q
  )
  pseudo_gt <- radmix_pseudohaploid_gt(svd_sim$genotypes, 1:20)
  adm <- gt_radmix(pseudo_gt, k = 3, seed = 1, init = "svd")
  expect_true(adm$converged)
})

test_that("invalid uses of the SVD start are caught", {
  expect_error(
    gt_radmix(svd_gt, k = 3, init = "svd", init_q = svd_sim$q),
    "cannot be combined"
  )
  expect_error(
    gt_radmix(svd_gt[1:8, ], k = 3, init = "svd"),
    "more than max(k) + 5 = 8",
    fixed = TRUE
  )
  expect_error(gt_radmix(svd_gt, k = 3, init = "pca"), "should be one of")
})

test_that("SVD fits report their progress as messages", {
  expect_message(
    gt_radmix(svd_gt, k = 2, seed = 1, init = "svd", verbose = TRUE),
    "Computing the SVD"
  )
})

# Genotype likelihoods (through the example source, see radmix_test_gl()):
# the SVD of the posterior dosages (src/radmix_gl_svd.cpp and
# radmix_gl_svd_factors())
test_that("the SVD start works with genotype likelihoods", {
  # likelihoods with missing data, and their reference in R: the emMAF
  # frequencies and the standardised posterior dosages
  set.seed(3)
  ll <- radmix_simulate_gl(svd_sim$genotypes, depth = 2)
  ll[cbind(1:5, 1:5, 1)] <- 0
  ll[cbind(1:5, 1:5, 2)] <- 0
  ll[cbind(1:5, 1:5, 3)] <- 0
  gl_dat <- radmix_test_gl(ll)
  expect_true(any(!gl_dat$obs))
  dosage <- function(f) {
    fm <- matrix(f, nrow(gl_dat$e1), length(f), byrow = TRUE)
    w1 <- 2 * gl_dat$e1 * fm * (1 - fm)
    w2 <- gl_dat$e2 * fm^2
    (w1 + 2 * w2) / ((1 - fm)^2 + w1 + w2)
  }
  f <- rep(0.25, ncol(gl_dat$e1))
  for (it in 1:100) {
    f <- colSums(dosage(f) * gl_dat$obs) / (2 * colSums(gl_dat$obs))
  }
  freq <- radmix_gl_freq_cpp(gl_dat, 1, tol = 1e-12, max_iter = 1000)
  expect_equal(freq, f, tolerance = 1e-8)
  expect_identical(radmix_gl_freq_cpp(gl_dat, 2, 1e-12, 1000), freq)
  fm <- matrix(freq, nrow(gl_dat$e1), length(freq), byrow = TRUE)
  x <- (dosage(freq) - 2 * fm) / sqrt(2 * fm * (1 - fm)) * gl_dat$obs
  b <- matrix(stats::rnorm(300 * 4), 300, 4)
  expect_equal(radmix_gl_dosage_prod_cpp(gl_dat, freq, b, FALSE, 1), x %*% b)
  b2 <- matrix(stats::rnorm(60 * 4), 60, 4)
  expect_equal(radmix_gl_dosage_prod_cpp(gl_dat, freq, b2, TRUE, 1),
               t(x) %*% b2)
  expect_identical(
    radmix_gl_dosage_prod_cpp(gl_dat, freq, b2, TRUE, 3),
    radmix_gl_dosage_prod_cpp(gl_dat, freq, b2, TRUE, 1)
  )
  # the randomised SVD gives the leading components of the dosages (those
  # of the structure, K - 1 = 2; the others are noise, with small gaps
  # between their singular values), and nearly the best approximation of
  # rank 5
  fac <- radmix_gl_svd_factors(gl_dat, rank = 5)
  sv <- svd(x, nu = 5, nv = 5)
  fac8 <- radmix_gl_svd_factors(gl_dat, rank = 5, n_power = 8)
  expect_equal(abs(fac8$u[, 1:2]), abs(sv$u[, 1:2]), tolerance = 1e-5)
  fm_t <- t(matrix(freq, nrow(x), ncol(x), byrow = TRUE))
  scale <- sqrt(2 * fm_t * (1 - fm_t))
  x5 <- t(fac$z %*% t(fac$u) / scale)
  best5 <- sv$u %*% diag(sv$d[1:5]) %*% t(sv$v)
  expect_lt(sum((x - x5)^2), 1.01 * sum((x - best5)^2))
  # (the frequencies, to the default tolerance of their EM)
  expect_equal(fac$f, freq, tolerance = 1e-5)
  # the likelihoods of hard calls without missing data give the start of
  # the genotypes, and fits from it reach the same maximum
  hard_ll <- array(-Inf, c(dim(svd_sim$genotypes), 3))
  for (gg in 0:2) hard_ll[, , gg + 1][svd_sim$genotypes == gg] <- 0
  test_gl <- svd_gt
  class(test_gl) <- c("test_svd_likelihoods", class(test_gl))
  local_mocked_s3_method(
    "gt_radmix_data", "test_svd_likelihoods",
    function(x, ...) radmix_test_gl(hard_ll)
  )
  fac_gl <- radmix_gl_svd_factors(radmix_test_gl(hard_ll),
    rank = 5, n_power = 8
  )
  fac_gt <- radmix_svd_factors(svd_gt, rank = 5)
  expect_equal(fac_gl$f, fac_gt$f, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(abs(fac_gl$u[, 1:2]), abs(fac_gt$u[, 1:2]), tolerance = 1e-5)
  # (a stringent rtol, as the log-likelihoods differ by a constant, which
  # changes the relative changes that the convergence criterion tests)
  adm_gl <- gt_radmix(test_gl, k = 3, seed = 1, init = "svd", rtol = 1e-12)
  adm_gt <- gt_radmix(svd_gt, k = 3, seed = 1, init = "svd", rtol = 1e-12)
  expect_true(adm_gl$converged)
  expect_match(adm_gl$log[[1]][1], "^SVD start")
  expect_equal(adm_gl$loglik, adm_gt$loglik + sum(svd_sim$genotypes == 1) *
                 log(2), tolerance = 1e-6)
  # more power iterations give the same maximum; svd_power must be a
  # count
  adm_gl8 <- gt_radmix(test_gl,
    k = 3, seed = 1, init = "svd", rtol = 1e-12, svd_power = 8
  )
  expect_equal(adm_gl8$loglik, adm_gl$loglik, tolerance = 1e-6)
  for (bad in list(-1, 1.5, NA_real_, "3")) {
    expect_error(
      gt_radmix(test_gl, k = 3, init = "svd", svd_power = bad),
      "'svd_power' should be an integer of at least 0"
    )
  }
  # reproducible with the same seed, whatever the number of cores
  expect_identical(
    gt_radmix(test_gl,
      k = 3, seed = 1, init = "svd", rtol = 1e-12, n_cores = 2
    )$Q,
    adm_gl$Q
  )
})

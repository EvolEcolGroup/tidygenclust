skip_on_cran()
# skip if the optional tool is not installed
skip_if_not(tgc_tool_installed("adamixture"))

library(tidypopgen)

# set up the gentibble
vcf_path <- system.file(
  "/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
  package = "tidypopgen"
)
anole_gt <- gen_tibble(
  vcf_path,
  quiet = TRUE,
  backingfile = tempfile("anolis_"),
  parser = "cpp"
)
pops_path <- system.file(
  "/extdata/anolis/plot_order_punctatus_n46.csv",
  package = "tidypopgen"
)
pops <- readr::read_csv(pops_path, show_col_types = FALSE)
anole_gt <- anole_gt %>% mutate(id = gsub("punc_", "", .data$id, ))
anole_gt <- anole_gt %>%
  mutate(population = pops$pop[match(pops$ID, .data$id)])

# check a gt_admix object returned by gt_adamixture() for anole_gt: one Q
# matrix (and one P matrix, with `freqs`) per combination of k and run, in
# order of k and then of run, and one cross-validation error per fit (with
# `cv`)
expect_anole_adamixture <- function(res, k, n_runs, freqs, cv = FALSE) {
  n_fits <- length(k) * n_runs
  testthat::expect_s3_class(res, "gt_admix")
  testthat::expect_identical(res$algorithm, "adamixture")
  testthat::expect_identical(res$id, anole_gt$id)
  testthat::expect_identical(res$k, rep(as.integer(k), each = n_runs))
  testthat::expect_length(res$Q, n_fits)
  for (i in seq_len(n_fits)) {
    testthat::expect_s3_class(res$Q[[i]], "q_matrix")
    testthat::expect_identical(
      dim(res$Q[[i]]), c(nrow(anole_gt), res$k[i])
    )
    testthat::expect_equal(
      unname(rowSums(res$Q[[i]])), rep(1, nrow(anole_gt)),
      tolerance = 1e-6
    )
  }
  if (freqs) {
    testthat::expect_length(res$P, n_fits)
    for (i in seq_len(n_fits)) {
      testthat::expect_identical(
        dim(res$P[[i]]), c(count_loci(anole_gt), res$k[i])
      )
      testthat::expect_true(all(res$P[[i]] >= 0 & res$P[[i]] <= 1))
    }
  } else {
    testthat::expect_false("P" %in% names(res))
  }
  if (cv) {
    testthat::expect_length(res$cv, n_fits)
    testthat::expect_true(all(is.finite(res$cv) & res$cv > 0))
  } else {
    testthat::expect_false("cv" %in% names(res))
  }
}

test_that("gt_adamixture", {
  # multiple k and one run, without P matrices
  res <- gt_adamixture(anole_gt, k = 2:3, verbose = FALSE)
  expect_anole_adamixture(res, k = 2:3, n_runs = 1, freqs = FALSE)

  # one k and multiple runs, with P matrices; the runs use different seeds
  res <- gt_adamixture(anole_gt,
    k = 3, n_runs = 2, seed = c(42, 31), no_freqs = FALSE, verbose = FALSE
  )
  expect_anole_adamixture(res, k = 3, n_runs = 2, freqs = TRUE)

  # the same seed gives the same result
  res1 <- gt_adamixture(anole_gt, k = 3, seed = 1, verbose = FALSE)
  res2 <- gt_adamixture(anole_gt, k = 3, seed = 1, verbose = FALSE)
  expect_identical(res1$Q[[1]], res2$Q[[1]])

  # a grouped gen_tibble adds the groups
  res <- gt_adamixture(dplyr::group_by(anole_gt, .data$population),
    k = 2, verbose = FALSE
  )
  expect_identical(res$group, anole_gt$population)
})

test_that("gt_adamixture with cv", {
  res <- gt_adamixture(anole_gt,
    k = 2:4, n_runs = 2, seed = c(1, 2), cv = 5, verbose = FALSE
  )
  expect_anole_adamixture(res, k = 2:4, n_runs = 2, freqs = FALSE, cv = TRUE)
  res <- gt_adamixture(anole_gt,
    k = 3, cv = 5, no_freqs = FALSE, verbose = FALSE
  )
  expect_anole_adamixture(res, k = 3, n_runs = 1, freqs = TRUE, cv = TRUE)
})

test_that("gt_adamixture output works with tidypopgen methods", {
  res <- gt_adamixture(dplyr::group_by(anole_gt, .data$population),
    k = 2:3, cv = 5, verbose = FALSE
  )
  expect_s3_class(autoplot(res, type = "cv"), "ggplot")
  expect_s3_class(autoplot(res, type = "barplot", k = 3, run = 1), "ggplot")
  tidy_q <- tidy(res$Q[[2]], anole_gt)
  expect_equal(nrow(tidy_q), nrow(anole_gt) * 3)
})

test_that("gt_adamixture uses the individuals and loci of a subset", {
  # filtering, reordering and select_loci() only change the indices of a
  # gen_tibble into its backing file: the results must be the same as for
  # PLINK files written from the subset
  idx <- c(40:7, 1:3)
  loci_idx <- rev(seq(1, count_loci(anole_gt), by = 2))
  sub_gt <- select_loci(anole_gt[idx, ], tidyselect::all_of(loci_idx))
  bed_file <- gt_as_plink(sub_gt, file = tempfile("anolis_sub_"))
  plink_prefix <- sub("\\.bed$", "", bed_file)
  res_gt <- gt_adamixture(sub_gt,
    k = 3, seed = 42, no_freqs = FALSE, verbose = FALSE
  )
  res_plink <- gt_adamixture(plink_prefix,
    k = 3, seed = 42, no_freqs = FALSE, verbose = FALSE
  )
  expect_identical(dim(res_gt$Q[[1]]), c(length(idx), 3L))
  expect_identical(dim(res_gt$P[[1]]), c(length(loci_idx), 3L))
  expect_equal(
    unclass(res_gt$Q[[1]]), unclass(res_plink$Q[[1]]),
    tolerance = 1e-6, ignore_attr = TRUE
  )
  expect_equal(res_gt$P[[1]], res_plink$P[[1]], tolerance = 1e-6)
  expect_identical(res_gt$id, anole_gt$id[idx])
})

test_that("P gives the frequencies of the alternate allele", {
  # P should predict the observed frequencies of the alternate allele, also
  # when it is the most common allele overall (ADAMIXTURE then flips the
  # coding, and returns the frequencies of the other allele)
  check_p <- function(gt) {
    res <- gt_adamixture(gt, k = 2, seed = 1, no_freqs = FALSE,
      verbose = FALSE
    )
    obs_freq <- colMeans(show_genotypes(gt), na.rm = TRUE) / 2
    pred_freq <- res$P[[1]] %*% colMeans(res$Q[[1]])
    expect_gt(cor(obs_freq, as.vector(pred_freq)), 0.9)
    # (not exactly equal with missing genotypes; the frequencies of the
    # other allele would give 1 - mean(obs_freq))
    expect_equal(mean(pred_freq), mean(obs_freq), tolerance = 0.05)
  }
  # in the anolis data, the alternate allele is the less common one
  check_p(anole_gt)
  # simulated data in which it is the most common one
  set.seed(1)
  p <- matrix(stats::runif(2 * 300, 0.6, 0.95), 300, 2)
  pop <- rep(1:2, each = 30)
  genotypes <- t(vapply(pop, function(j) {
    stats::rbinom(300, 2, p[, j])
  }, numeric(300)))
  common_alt <- radmix_test_gt(genotypes)
  expect_gt(mean(show_genotypes(common_alt)) / 2, 0.5)
  check_p(common_alt)
})

test_that("gt_adamixture uses the requested number of threads", {
  omp_vars <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS"
  )
  # unset the variables for this test (earlier calls may have set them)
  withr::local_envvar(stats::setNames(rep(NA, length(omp_vars)), omp_vars))
  timing <- system.time(
    gt_adamixture(anole_gt,
      k = 6, n_runs = 2, seed = 1:2, threads = 1, verbose = FALSE
    )
  )
  # the environment of the R session is not changed (ADAMIXTURE sets these
  # variables)
  expect_true(all(is.na(Sys.getenv(omp_vars, unset = NA))))
  # with a single thread, the CPU time cannot exceed the elapsed time by
  # much
  cpu <- timing[["user.self"]] + timing[["sys.self"]]
  expect_lt(cpu / timing[["elapsed"]], 1.5)
})

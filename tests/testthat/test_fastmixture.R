skip_on_cran()
# skip if the optional tool is not installed
skip_if_not(tgc_tool_installed("fastmixture"))

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
pops <- readr::read_csv(pops_path)
anole_gt <- anole_gt %>% mutate(id = gsub("punc_", "", .data$id, ))
anole_gt <- anole_gt %>%
  mutate(population = pops$pop[match(pops$ID, .data$id)])


# check a gt_admix object returned by gt_fastmixture() for anole_gt: one Q
# matrix (and one P matrix, with `freqs`) per combination of k and run, in
# order of k and then of run, and one cross-validation error per fit (with
# `cv`)
expect_anole_fastmixture <- function(res, k, n_runs, freqs, cv = FALSE) {
  n_fits <- length(k) * n_runs
  testthat::expect_s3_class(res, "gt_admix")
  testthat::expect_identical(res$algorithm, "fastmixture")
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

test_that("gt_fastmixture", {
  # multiple k and one run, without P matrices (with the arguments that are
  # passed to fastmixture set explicitly, all to their defaults except for
  # the tolerance, to keep the tests fast)
  res <- gt_fastmixture(
    anole_gt,
    k = 2:3,
    n_runs = 1,
    threads = 1,
    seed = 42,
    iter = 1000,
    tole = 0.5,
    batches = 32,
    supervised = NULL,
    check = 5,
    power = 11,
    chunk = 8192,
    subsample = 0.7,
    min_subsample = 50000,
    max_subsample = 500000,
    als_iter = 1000,
    als_tole = 1e-4,
    no_freqs = TRUE,
    random_init = FALSE
  )
  expect_anole_fastmixture(res, k = 2:3, n_runs = 1, freqs = FALSE)

  # one k and multiple runs, with P matrices; the runs use different seeds
  res <- gt_fastmixture(
    anole_gt, k = 2, n_runs = 2, seed = c(42, 31), tole = 0.5,
    no_freqs = FALSE
  )
  expect_anole_fastmixture(res, k = 2, n_runs = 2, freqs = TRUE)
  expect_false(isTRUE(all.equal(res$Q[[1]], res$Q[[2]])))

  # multiple k and multiple runs, with P matrices
  res <- gt_fastmixture(
    anole_gt, k = 2:3, n_runs = 2, seed = c(42, 31), tole = 0.5,
    no_freqs = FALSE
  )
  expect_anole_fastmixture(res, k = 2:3, n_runs = 2, freqs = TRUE)

  # one k and one run, with P matrices, for a grouped gen_tibble
  res <- gt_fastmixture(
    dplyr::group_by(anole_gt, .data$population), k = 2, tole = 0.5,
    no_freqs = FALSE
  )
  expect_anole_fastmixture(res, k = 2, n_runs = 1, freqs = TRUE)
  expect_identical(res$group, anole_gt$population)
})

test_that("fastmixture with cv", {
  # one k and one run, without and with P matrices
  for (no_freqs in c(TRUE, FALSE)) {
    res <- gt_fastmixture(
      anole_gt, k = 3, seed = 42, tole = 0.5, no_freqs = no_freqs,
      cv = 7, cv_tole = 1e-7
    )
    expect_anole_fastmixture(res, k = 3, n_runs = 1, freqs = !no_freqs,
      cv = TRUE
    )
  }

  # multiple k and multiple runs, without P matrices
  res <- gt_fastmixture(
    anole_gt, k = 2:4, n_runs = 2, seed = c(123, 234), tole = 0.5,
    cv = 7
  )
  expect_anole_fastmixture(res, k = 2:4, n_runs = 2, freqs = FALSE, cv = TRUE)

  # multiple k and multiple runs, with P matrices
  res <- gt_fastmixture(
    anole_gt, k = 2:7, n_runs = 4, seed = c(123, 234, 345, 456), tole = 0.5,
    no_freqs = FALSE, cv = 5
  )
  expect_anole_fastmixture(res, k = 2:7, n_runs = 4, freqs = TRUE, cv = TRUE)
})

test_that("gt_fastmixture gives the same results for gen_tibble and plink", {
  bed_file <- gt_as_plink(anole_gt, file = tempfile("anolis_plink_"))
  plink_prefix <- sub("\\.bed$", "", bed_file)
  res_gt <- gt_fastmixture(
    anole_gt,
    k = 3,
    seed = 42,
    no_freqs = FALSE,
    random_init = FALSE
  )
  res_plink <- gt_fastmixture(
    plink_prefix,
    k = 3,
    seed = 42,
    no_freqs = FALSE,
    random_init = FALSE
  )
  # both inputs are shuffled and normalised in the same way, so the
  # estimates should be identical
  expect_equal(
    unclass(res_gt$Q[[1]]),
    unclass(res_plink$Q[[1]]),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  expect_equal(res_gt$P[[1]], res_plink$P[[1]], tolerance = 1e-6)
  # P should be returned in the original locus order, so it should predict
  # the observed allele frequencies
  obs_freq <- colMeans(show_genotypes(anole_gt), na.rm = TRUE) / 2
  pred_freq <- res_gt$P[[1]] %*% colMeans(res_gt$Q[[1]])
  expect_gt(cor(obs_freq, as.vector(pred_freq)), 0.9)
})

test_that("gt_fastmixture uses the individuals and loci of a subset", {
  # filtering, reordering and select_loci() only change the indices of a
  # gen_tibble into its backing file: the results must be the same as for
  # PLINK files written from the subset
  idx <- c(40:7, 1:3)
  # loci that are polymorphic in the subset (fastmixture requires them), in
  # reverse order
  maf <- loci_maf(anole_gt[idx, ])
  loci_idx <- rev(which(maf > 0.05))
  sub_gt <- select_loci(anole_gt[idx, ], tidyselect::all_of(loci_idx))
  bed_file <- gt_as_plink(sub_gt, file = tempfile("anolis_sub_"))
  plink_prefix <- sub("\\.bed$", "", bed_file)
  res_gt <- gt_fastmixture(
    sub_gt,
    k = 3,
    seed = 42,
    no_freqs = FALSE,
    random_init = FALSE
  )
  res_plink <- gt_fastmixture(
    plink_prefix,
    k = 3,
    seed = 42,
    no_freqs = FALSE,
    random_init = FALSE
  )
  expect_identical(dim(res_gt$Q[[1]]), c(length(idx), 3L))
  expect_identical(dim(res_gt$P[[1]]), c(length(loci_idx), 3L))
  expect_equal(
    unclass(res_gt$Q[[1]]),
    unclass(res_plink$Q[[1]]),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
  expect_equal(res_gt$P[[1]], res_plink$P[[1]], tolerance = 1e-6)
  expect_identical(res_gt$id, anole_gt$id[idx])
})

test_that("monomorphic loci need random_init = TRUE", {
  # some loci are monomorphic in the first 5 individuals
  sub_gt <- anole_gt[1:5, ]
  expect_true(any(loci_maf(sub_gt) == 0))
  # the SVD-based initialisation (the default) needs polymorphic loci
  expect_error(
    gt_fastmixture(sub_gt, k = 2, seed = 42, tole = 0.5),
    "random_init = TRUE"
  )
  # random starting values do not
  res <- gt_fastmixture(sub_gt, k = 2, seed = 42, tole = 0.5,
    random_init = TRUE
  )
  expect_identical(dim(res$Q[[1]]), c(5L, 2L))
})

test_that("gt_fastmixture uses the requested number of threads", {
  omp_vars <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS"
  )
  # unset the variables for this test (earlier calls may have set them)
  withr::local_envvar(stats::setNames(rep(NA, length(omp_vars)), omp_vars))
  timing <- system.time(
    gt_fastmixture(anole_gt, k = 6, n_runs = 2, seed = 1:2, threads = 1)
  )
  # the environment of the R session is not changed
  expect_true(all(is.na(Sys.getenv(omp_vars, unset = NA))))
  # with a single thread, the CPU time cannot exceed the elapsed time by
  # much (it was several times larger when threads was ignored)
  cpu <- timing[["user.self"]] + timing[["sys.self"]]
  expect_lt(cpu / timing[["elapsed"]], 1.5)
})

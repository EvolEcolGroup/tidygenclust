# Streaming of data larger than the memory (src/radmix_geno_access.h,
# "Streaming larger-than-memory data", and src/radmix_kernels.cpp,
# "Streaming"): the passes walk the loci in chunks, prefetch the next one
# and alternate their direction, with results bit-identical to those of the
# usual layout. Tiny chunks (chunk_bytes) force many of them.

# the data with streaming on, in chunks of about `chunk_bytes` bytes
stream_data <- function(dat, chunk_bytes = 1) {
  dat$stream <- TRUE
  dat$chunk_bytes <- chunk_bytes
  dat
}

# every pass of the kernels, as a list of results, on data `dat` (with
# individuals `fixed` fixed, or none)
all_passes <- function(dat, q, p, fixed = NULL) {
  list(
    loglik = radmix_loglik_cpp(dat, q, p, 2),
    deviance = radmix_deviance_cpp(dat, q, p, 2),
    em = radmix_em_cpp(dat, q, p, 2),
    grad_q = radmix_grad_q_cpp(dat, q, p, 2),
    grad_p = radmix_grad_p_cpp(dat, q, p, 2),
    update_q = radmix_update_q_cpp(dat, q, p, 2, fixed),
    update_q_loglik = radmix_update_q_loglik_cpp(dat, q, p, 2, fixed),
    update_p = radmix_update_p_cpp(dat, q, p, 2)
  )
}

random_qp <- function(k, n, m, seed) {
  set.seed(seed)
  q <- matrix(stats::runif(k * n), k, n)
  q <- sweep(q, 2, colSums(q), "/")
  p <- matrix(stats::runif(k * m, 0.05, 0.95), k, m)
  list(q = q, p = p)
}

sim <- radmix_simulate_admixed(n = 150, m = 1100, k = 3, seed = 3)
sim_gt <- radmix_test_gt(sim$genotypes)
dat <- radmix_prepare(sim_gt)

test_that("streaming passes are identical (Q side in segments)", {
  # 150 individuals: 3 blocks, so the Q side splits the loci into segments
  # (4 of 275 loci)
  expect_identical(radmix_q_segments_cpp(150, 1100), 4L)
  qp <- random_qp(3, 150, 1100, 1)
  ref <- all_passes(dat, qp$q, qp$p)
  # chunks of one block of loci (and one segment), of two segments, and of
  # all the loci
  for (cb in c(1, 150 * 640, 1e9)) {
    expect_identical(all_passes(stream_data(dat, cb), qp$q, qp$p), ref)
  }
  # a supervised analysis (fixed individuals) and a training view of
  # cross-validation
  fixed <- rep(c(TRUE, FALSE, FALSE), 50)
  expect_identical(
    all_passes(stream_data(dat), qp$q, qp$p, fixed),
    all_passes(dat, qp$q, qp$p, fixed)
  )
  cv <- list(mode = "train", fold = 2, n_folds = 5, seed = 11)
  expect_identical(
    all_passes(stream_data(radmix_set_cv(dat, cv)), qp$q, qp$p),
    all_passes(radmix_set_cv(dat, cv), qp$q, qp$p)
  )
})

test_that("streaming passes are identical (Q side in one segment)", {
  skip_on_cran()
  # 8200 individuals (more than MIN_Q_TASKS blocks): one segment, so a
  # streaming Q side carries the sums of each block from chunk to chunk;
  # 201 loci give chunks of 64 loci and an odd last one
  set.seed(5)
  n <- 8200
  m <- 201
  g <- matrix(stats::rbinom(n * m, 2, 0.3), n, m)
  g[sample(length(g), 2000)] <- NA
  big <- radmix_prepare(radmix_test_gt(g))
  expect_identical(radmix_q_segments_cpp(n, m), 1L)
  qp <- random_qp(2, n, m, 2)
  ref <- all_passes(big, qp$q, qp$p)
  expect_identical(all_passes(stream_data(big), qp$q, qp$p), ref)
  fixed <- seq_len(n) %% 7 == 0
  expect_identical(
    all_passes(stream_data(big), qp$q, qp$p, fixed),
    all_passes(big, qp$q, qp$p, fixed)
  )
})

test_that("streaming passes are identical for genotype likelihoods", {
  gl <- radmix_test_gl(radmix_simulate_gl(sim$genotypes))
  qp <- random_qp(3, 150, 1100, 4)
  # without locus_bytes, a chunk holds chunk_bytes loci (here 70)
  expect_identical(
    all_passes(stream_data(gl, 70), qp$q, qp$p),
    all_passes(gl, qp$q, qp$p)
  )
})

test_that("fits that stream are identical to those that do not", {
  fit <- function(stream) {
    withr::local_options(
      tidygenclust.radmix_stream = stream,
      tidygenclust.radmix_chunk_mb = 1e-4
    )
    gt_radmix(sim_gt, k = 3, seed = 1, n_cores = 2)
  }
  expect_identical(fit(TRUE), fit(FALSE))
})

test_that("streaming is chosen from the size of the data", {
  withr::local_options(tidygenclust.radmix_stream = NULL)
  # these data are small
  s <- radmix_stream_mode(dat)
  expect_false(s$stream)
  expect_identical(s$data_bytes, 150 * 1100)
  expect_gte(s$chunk_bytes, 16 * 2^20)
  expect_lte(s$chunk_bytes, 256 * 2^20)
  # the options force it, and set the size of the chunks
  withr::local_options(
    tidygenclust.radmix_stream = TRUE,
    tidygenclust.radmix_chunk_mb = 2
  )
  s <- radmix_stream_mode(dat)
  expect_true(s$stream)
  expect_identical(s$chunk_bytes, 2 * 2^20)
  withr::local_options(tidygenclust.radmix_chunk_mb = -1)
  expect_error(radmix_stream_mode(dat), "radmix_chunk_mb")
  # the memory available is a positive number of bytes
  expect_gt(radmix_mem_available(), 0)
  # a source of genotype likelihoods without locus_bytes has no known size
  gl <- radmix_test_gl(radmix_simulate_gl(sim$genotypes[1:5, 1:5]))
  expect_identical(radmix_data_bytes(gl), NA_real_)
  gl$locus_bytes <- 20
  expect_identical(radmix_data_bytes(gl), 100)
})

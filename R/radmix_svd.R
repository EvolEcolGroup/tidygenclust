#' The factors of a truncated SVD of the genotypes, for the starts of
#' gt_radmix(init = "svd")
#'
#' Computes a partial SVD of the genotypes of `x` with
#' [tidypopgen::gt_pca_partialSVD()], each locus being centred by twice the
#' frequency f of its alternate allele and scaled by sqrt(2 f (1 - f)), as
#' with the default scaling of tidypopgen's PCAs (Patterson et al. 2006), and
#' returns it as factors of the expected allele frequencies of the
#' individuals. With G the genotypes (individuals x loci), the SVD gives
#' (G - 2 1 f') S^-1 ~ U D V', with S the diagonal matrix of the scales, so
#' the expected frequencies G' / 2 (loci x individuals) are approximately
#' f 1' + 0.5 Z U', with Z = S V D.
#'
#' Missing genotypes are not allowed in the SVD: they must have been imputed
#' (e.g. with [tidypopgen::gt_impute_simple()]), and the SVD then uses the
#' imputed values (gt_pca_partialSVD() switches them on). f is the frequency
#' among the observed genotypes, as the allele frequencies of ADAMIXTURE's
#' start. Monomorphic loci, which tidypopgen's default scaling rejects, are
#' left unscaled: their centred genotypes are all zero, so they do not
#' contribute to the SVD, and their rows of Z are zero.
#' @param x a [`tidypopgen::gen_tibble`]
#' @param rank the number of components
#' @param n_cores the number of cores used for the allele frequencies and the
#'   missingness
#' @returns a list with elements `z` (loci x rank), `u` (individuals x
#'   rank) and `f` (the frequency of the alternate allele at each locus), in
#'   the order of the loci and individuals of `x`
#' @keywords internal
#' @noRd
radmix_svd_factors <- function(x, rank, n_cores = 1) {
  # the allele frequencies of a grouped gen_tibble would be computed by
  # group
  x <- dplyr::ungroup(x)
  # (gt_has_imputed() tells whether the backing file holds imputed values,
  # which then replace all the missing genotypes)
  has_missing <- !tidypopgen::gt_has_imputed(x) &&
    any(tidypopgen::loci_missingness(x, n_cores = n_cores) > 0)
  if (has_missing) {
    stop(
      "init = \"svd\" needs a gen_tibble without missing genotypes, or ",
      "with imputed ones: impute them first with ",
      "tidypopgen::gt_impute_simple() (the imputed genotypes are only used ",
      "for the starting values; the fit ignores them, unless ",
      "'use_imputed = TRUE')",
      call. = FALSE
    )
  }
  f <- tidypopgen::loci_alt_freq(x, n_cores = n_cores)
  # (a locus without any observed genotype, whose imputed values are
  # arbitrary)
  f[!is.finite(f)] <- 0.5
  scale <- sqrt(2 * f * (1 - f))
  scale[scale == 0] <- 1
  # bigstatsr calls the scaling function with the FBM and the indices of
  # the selected individuals and loci, for one block of loci at a time
  # (blocks of bigstatsr::block_size() loci: a single block for small data,
  # several when n_indiv x n_loci exceeds ~1.3e8), and expects the values of
  # the loci in `ind.col` only; f and scale are in the order of the loci of
  # x, i.e. of their columns in the FBM
  fbm_cols <- tidypopgen::show_loci(x)$big_index
  scaling <- function(X, ind.row, ind.col) { # nolint: object_name_linter
    pos <- match(ind.col, fbm_cols)
    data.frame(center = 2 * f[pos], scale = scale[pos])
  }
  pca <- tidypopgen::gt_pca_partialSVD(x,
    k = rank, fun_scaling = scaling, total_var = FALSE
  )
  list(
    z = unname(sweep(pca$v, 2, pca$d, "*") * scale),
    u = unname(pca$u),
    f = f
  )
}

#' The factors of a truncated SVD of genotype likelihoods, for the starts of
#' gt_radmix(init = "svd")
#'
#' The counterpart of `radmix_svd_factors()` for genotype likelihoods (a
#' source given by a method of [gt_radmix_data()]): the genotypes are not
#' observed, so the SVD is that of their posterior dosages (expected
#' numbers of copies of `allele2`) under Hardy-Weinberg proportions with the
#' frequency f of `allele2` at each locus, estimated from the likelihoods by
#' EM (emMAF, Kim et al. 2011), centred by 2 f and scaled by
#' sqrt(2 f (1 - f)), as in the first step of PCAngsd (Meisner and
#' Albrechtsen 2018); an individual without data at a locus gets the mean
#' dosage, 2 f (src/radmix_gl_svd.cpp). For the likelihoods of hard calls
#' without missing genotypes, the dosages are the genotypes, and the SVD is
#' that of `radmix_svd_factors()`.
#'
#' The dosages are never stored: the truncated SVD is a randomised SVD
#' (Halko, Martinsson and Tropp 2011, algorithm 4.4) with `rank + 10`
#' vectors and `n_power` power iterations, whose products with the dosages
#' are passes over the likelihoods (2 n_power + 2 passes, plus one for the
#' frequencies). Its random test matrix is drawn from R's generator.
#' @param data the list of kind `"gl_source"` given by [gt_radmix_data()]
#' @param rank the number of components
#' @param n_cores the number of threads of the passes over the likelihoods
#' @param n_power the number of power iterations
#' @returns a list with elements `z` (loci x rank), `u` (individuals x
#'   rank) and `f` (the frequency of `allele2` at each locus), as
#'   `radmix_svd_factors()`
#' @references Halko N, Martinsson PG, Tropp JA (2011) Finding structure
#'   with randomness: probabilistic algorithms for constructing approximate
#'   matrix decompositions. SIAM Review 53: 217-288.
#'
#'   Kim SY, Lohmueller KE, Albrechtsen A, et al. (2011) Estimation of allele
#'   frequency and association mapping using next-generation sequencing
#'   data. BMC Bioinformatics 12: 231.
#'
#'   Meisner J, Albrechtsen A (2018) Inferring population structure and
#'   admixture proportions in low-depth NGS data. Genetics 210: 719-731.
#' @keywords internal
#' @noRd
radmix_gl_svd_factors <- function(data, rank, n_cores = 1, n_power = 3) {
  f <- radmix_gl_freq_cpp(data, n_cores)
  n_vec <- rank + 10
  # the range of the dosages: Y = X Omega, refined by power iterations
  # (each product orthonormalised first, for numerical stability)
  omega <- matrix(stats::rnorm(data$n_loci * n_vec), data$n_loci, n_vec)
  y <- radmix_gl_dosage_prod_cpp(data, f, omega, FALSE, n_cores)
  for (it in seq_len(n_power)) {
    z <- radmix_gl_dosage_prod_cpp(data, f, qr.Q(qr(y)), TRUE, n_cores)
    y <- radmix_gl_dosage_prod_cpp(data, f, qr.Q(qr(z)), FALSE, n_cores)
  }
  qy <- qr.Q(qr(y))
  # X ~ Qy Qy' X = Qy B', with B = X' Qy (loci x n_vec), whose SVD
  # B = W D V' gives X ~ (Qy V) D W'
  b <- radmix_gl_dosage_prod_cpp(data, f, qy, TRUE, n_cores)
  sv <- svd(b, nu = rank, nv = rank)
  scale <- sqrt(2 * f * (1 - f))
  scale[scale == 0] <- 1
  list(
    z = unname(sweep(sv$u, 2, sv$d[seq_len(rank)], "*") * scale),
    u = unname(qy %*% sv$v),
    f = f
  )
}

#' A start for gt_radmix(init = "svd")
#'
#' Fits `n_starts` ALS solutions (`radmix_als_cpp()`, src/radmix_als.cpp)
#' to the factors of the SVD with K + 5 components, each from random allele
#' frequencies drawn from R's random number generator, and returns the one
#' with the highest log-likelihood (one pass over the genotypes each). R's
#' BLAS is set to one thread while the ALS runs, and restored afterwards.
#'
#' The inputs of the ALS (the factors scaled by 0.5 and augmented with the
#' allele frequencies and a column of ones, and the weights of the loci)
#' are prepared as in ALS() of ADAMIXTURE 1.8.9 (adamixture/model/als.py,
#' Copyright (c) 2024, Ioannidis Lab, BSD 3-clause licence), from which the
#' ALS itself is adapted (see src/radmix_als.cpp and inst/COPYRIGHTS).
#' @param data the output of `radmix_prepare()`
#' @param factors the output of `radmix_svd_factors()`, with at least K + 5
#'   components
#' @param k the number of ancestral populations
#' @param n_starts the number of ALS solutions
#' @param n_cores the number of cores (OpenMP threads)
#' @param verbose boolean; should the log be printed as it comes?
#' @returns a list with elements `q` (individuals x K), `p` (loci x K) and
#'   `log` (lines to add to the log of the fit)
#' @keywords internal
#' @noRd
radmix_svd_start <- function(data, factors, k, n_starts = 3, n_cores = 1,
                             verbose = FALSE) {
  # ADAMIXTURE uses K components; K + 5 found better starts in our tests
  # (data-raw/radmix/explore_svd_start_refined.R), for little extra cost
  rank <- k + 5
  # Pi ~ za ua' (see radmix_svd_factors())
  za <- cbind(0.5 * factors$z[, seq_len(rank), drop = FALSE], factors$f)
  ua <- cbind(factors$u[, seq_len(rank), drop = FALSE], 1)
  # the weights of the loci in the Q step: the inverse of their binomial
  # variance, bounded for loci with rare alleles, and normalised to mean one
  # (as in ADAMIXTURE)
  w <- 1 / pmax(factors$f * (1 - factors$f), 1e-3)
  w <- w / mean(w)
  add_log <- function(log, line) {
    if (verbose) message(line)
    c(log, line)
  }
  log <- add_log(character(0), paste0(
    "SVD start: best of ", n_starts, " ALS solutions, with ", rank,
    " components"
  ))
  # the ALS with one BLAS thread: its products are tall and skinny, so BLAS
  # threads (OpenBLAS uses all the cores by default) gain nothing and compete
  # with the OpenMP threads of its bounded least squares (HGDP, 1043 x
  # 130,943, K = 10, n_cores = 8: 4.7 s with 24 OpenBLAS threads, 2.4 s with
  # one)
  blas_threads <- RhpcBLASctl::blas_get_num_procs()
  if (isTRUE(blas_threads > 1)) {
    RhpcBLASctl::blas_set_num_threads(1)
    on.exit(RhpcBLASctl::blas_set_num_threads(blas_threads), add = TRUE)
  }
  best <- NULL
  for (s in seq_len(n_starts)) {
    init_p <- matrix(stats::runif(data$n_loci * k), data$n_loci, k)
    als <- radmix_als_cpp(za, ua, w, init_p,
      max_iter = 1000, tol = 1e-4, n_cores = n_cores
    )
    als$loglik <- radmix_loglik_cpp(data, als$q, als$p, n_cores)
    als$index <- s
    log <- add_log(log, paste0(
      "ALS solution ", s, ": ", als$n_iter, " iterations",
      if (!als$converged) " (not converged)",
      ", loglikelihood = ", format(als$loglik, digits = 12)
    ))
    # (isTRUE() keeps the first solution if the log-likelihoods are NaN)
    if (is.null(best) || isTRUE(als$loglik > best$loglik)) {
      best <- als
    }
  }
  log <- add_log(log, paste0("Starting from ALS solution ", best$index))
  # the user layout (individuals x K and loci x K), as radmix_fit() expects
  list(q = t(best$q), p = t(best$p), log = log)
}

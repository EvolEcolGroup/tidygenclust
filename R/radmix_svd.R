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

#' A start for gt_radmix(init = "svd")
#'
#' Fits `n_starts` ALS solutions (`radmix_als_cpp()`, src/radmix_als.cpp)
#' to the factors of the SVD with K + 5 components, each from random allele
#' frequencies drawn from R's random number generator, and returns the one
#' with the highest log-likelihood (one pass over the genotypes each).
#' @param data the output of `radmix_prepare()`
#' @param factors the output of `radmix_svd_factors()`, with at least K + 5
#'   components
#' @param k the number of ancestral populations
#' @param n_starts the number of ALS solutions
#' @param n_cores the number of cores
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

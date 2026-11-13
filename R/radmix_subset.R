#' A start for gt_radmix(init = "subset")
#'
#' Fits the model to a random subset of the loci, from a random start, and
#' extends its estimates to all the loci: the Q of the subset fit, and, for
#' every locus, the allele frequencies that best fit the genotypes given
#' that Q. The fit of all the loci then starts close to a mode of the
#' likelihood, and needs fewer iterations than from a random start, while
#' the subset fit costs only a fraction of a pass over the genotypes per
#' pass (measurements: data-raw/radmix/newspeed_algo/).
#'
#' The steps:
#' 1. a random subset of `ceiling(frac * n_loci)` loci (at least K + 1, and
#'    at most all of them) is drawn from R's random number generator;
#' 2. the model is fitted to the genotypes of these loci by radmix_fit(),
#'    from a random start (also drawn from R's generator) primed by
#'    `n_em_iter` EM steps, to the relative tolerance `rtol_subset` (the
#'    fit of all the loci refines it, so it need not be tight);
#' 3. given its Q, P is estimated for every locus: one EM step from allele
#'    frequencies of 0.5 in every population, which gives the frequencies
#'    weighted by ancestry, p_kj = sum_i q_ik g_ij / sum_i n_i q_ik over the
#'    observed genotypes (the EM step of Q is discarded), then two Newton
#'    steps of P (the loci are independent given Q, so the Newton steps
#'    converge quickly). This costs three passes over the genotypes.
#'
#' @param data the output of `radmix_prepare()`
#' @param k the number of ancestral populations
#' @param frac the fraction of the loci in the subset
#' @param rtol_subset the convergence threshold of the subset fit
#' @param max_iter,n_em_iter,n_qn,n_cores passed to radmix_fit()
#' @param verbose boolean; should the log be printed as it comes?
#' @returns a list with elements `q` (individuals x K), `p` (loci x K) and
#'   `log` (lines to add to the log of the fit)
#' @keywords internal
#' @noRd
radmix_subset_start <- function(data, k, frac, rtol_subset = 1e-6,
                                max_iter = 1000, n_em_iter = 5, n_qn = 3,
                                n_cores = 1, verbose = FALSE) {
  add_log <- function(log, line) {
    if (verbose) message(line)
    c(log, line)
  }
  # 1. The subset, sorted so that the loci are read in the order of the
  # backing file. The kernels read the loci listed in `cols` (see
  # radmix_prepare()), so restricting `cols` gives a view of the subset
  # without copying any genotypes.
  n_sub <- min(data$n_loci, max(k + 1, ceiling(frac * data$n_loci)))
  idx <- sort(sample.int(data$n_loci, n_sub))
  sub <- data
  sub$cols <- data$cols[idx]
  sub$n_loci <- n_sub
  # 2. The fit of the subset, from a random start primed by EM
  fit_sub <- radmix_fit(sub,
    k = k, max_iter = max_iter, rtol = rtol_subset, n_em_iter = n_em_iter,
    n_qn = n_qn, n_cores = n_cores
  )
  log <- add_log(character(0), paste0(
    "Subset start: fit of ", n_sub, " of ", data$n_loci, " loci, ",
    fit_sub$n_iter, " iterations, loglikelihood = ",
    format(fit_sub$loglik, digits = 12)
  ))
  # 3. P for all the loci given the Q of the subset. The kernels take Q and
  # P in the internal layout (K x individuals and K x loci).
  q <- t(fit_sub$q)
  flat <- matrix(0.5, k, data$n_loci)
  p <- radmix_em_cpp(data, q, flat, n_cores)$p
  for (step in 1:2) {
    p <- radmix_update_p_cpp(data, q, p, n_cores)
  }
  log <- add_log(log, paste0(
    "Subset start: P estimated for all the loci given the Q of the subset ",
    "(1 EM and 2 Newton steps)"
  ))
  # the user layout (individuals x K and loci x K), as radmix_fit() expects
  list(q = fit_sub$q, p = t(p), log = log)
}

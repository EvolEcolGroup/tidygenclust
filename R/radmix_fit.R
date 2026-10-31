#' Fit the ADMIXTURE model with gt_radmix() for a single value of K
#'
#' This is the engine behind `gt_radmix()`. It follows
#' `admixture_qn!()` in OpenADMIXTURE.jl: a random start (or user supplied
#' starting values), a few EM iterations to prime the estimates, and then
#' block-wise constrained Newton updates of Q and P accelerated by the
#' quasi-Newton scheme of Zhou, Alexander and Lange (2011).
#'
#' Each iteration of the main loop computes two rounds of block updates
#' (Q, then P with the new Q), giving three successive estimates x0, x1 and
#' x2. The secant pairs (x1 - x0, x2 - x1) of the last `n_qn` iterations
#' define a quasi-Newton step from x1, which is projected back onto the
#' feasible set and accepted only if it improves on the previous
#' log-likelihood; otherwise x2 is used.
#'
#' @param data the output of `radmix_prepare()` (possibly with a
#'   cross-validation mask, see `radmix_set_cv()`)
#' @param k the number of ancestral populations
#' @param init_q optional matrix (individuals x K) of starting ancestry
#'   proportions
#' @param init_p optional matrix (loci x K) of starting allele frequencies
#' @param fix_q,fix_p booleans; if `TRUE`, the corresponding parameter block
#'   is kept at its starting value
#' @param max_iter maximum number of iterations of the main loop
#' @param rtol convergence threshold on the relative change in log-likelihood
#' @param n_em_iter number of EM iterations used to prime the estimates
#'   (skipped if both `init_q` and `init_p` are given)
#' @param n_qn number of secant pairs used by the quasi-Newton step
#' @param n_cores number of threads
#' @param verbose boolean; print progress
#' @returns a list with elements `q` (individuals x K), `p` (loci x K),
#'   `loglik`, `loglik_history`, `n_iter`, `converged` and `log`
#' @keywords internal
#' @noRd
radmix_fit <- function(data, k, init_q = NULL, init_p = NULL,
                       fix_q = FALSE, fix_p = FALSE, max_iter = 1000,
                       rtol = 1e-7, n_em_iter = 5, n_qn = 3, n_cores = 1,
                       verbose = FALSE) {
  # The whole fit runs in C++ (src/radmix_fit.cpp), which documents the
  # model and the algorithm step by step. Running the loop in C++ means that the
  # estimates live in buffers allocated once and reused at every iteration:
  # when the loop was written in R, each pass over the genotypes returned new
  # R matrices, and the discarded ones (each a copy of P, with many loci)
  # accumulated until R's garbage collector ran, more than doubling the peak
  # memory of a fit.
  #
  # The random start is drawn from R's random number generator inside
  # radmix_fit_cpp(), so set.seed() makes fits reproducible
  # (gt_radmix() handles the seeds). Starting values are passed in the
  # user layout (individuals x K and loci x K), as are the estimates returned.
  # The C++ function takes them as numeric matrices or NULL. (Rcpp would also
  # convert an integer matrix, so this conversion only makes the storage
  # mode explicit: the callers in the package already pass doubles.)
  as_double_matrix <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    storage.mode(x) <- "double"
    x
  }
  radmix_fit_cpp(
    data,
    k = k,
    init_q = as_double_matrix(init_q),
    init_p = as_double_matrix(init_p),
    fix_q = fix_q, fix_p = fix_p, max_iter = max_iter, rtol = rtol,
    n_em_iter = n_em_iter, n_qn = n_qn, n_cores = n_cores,
    verbose = verbose
  )
}

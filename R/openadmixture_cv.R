#' Cross-validation error of a fitted model, as in ADMIXTURE
#'
#' Implements the cross-validation of ADMIXTURE (Alexander and Lange 2011):
#' the observed genotypes are split into `n_folds` folds; for each fold, the
#' model is refitted with the genotypes of that fold masked (treated as
#' missing), and the refitted model is used to predict them. The prediction
#' error is the binomial deviance of the masked genotypes, averaged over all
#' masked genotypes of all folds.
#'
#' @param data the output of `openadmixture_prepare()`
#' @param fit the fit to the full data, as returned by `openadmixture_fit()`,
#' used as
#'   the starting point of each fold
#' @param k the number of ancestral populations
#' @param n_folds the number of folds
#' @param cv_seed seed of the assignment of genotypes to folds (a
#'   non-negative integer)
#' @param max_iter,rtol,n_qn,n_cores passed to `openadmixture_fit()`
#' @returns a list with elements `cv` (the cross-validation error),
#'   `fold_cv` (the error of each fold), `fold_n` (the number of genotypes in
#'   each fold) and `log` (log lines, one per fold plus a summary line in the
#'   format used by ADMIXTURE)
#' @keywords internal
#' @noRd
openadmixture_cv <- function(data, fit, k, n_folds, cv_seed, max_iter = 1000,
                             rtol = 1e-7, n_qn = 3, n_cores = 1) {
  fold_dev <- numeric(n_folds)
  fold_n <- numeric(n_folds)
  log_lines <- character(0)
  for (fold in seq_len(n_folds)) {
    # Fit the model without the genotypes of this fold. The training view
    # hides them from every kernel, exactly as if they were missing (see
    # src/openadmixture_geno_access.h). As ADMIXTURE does not document how its
    # fold fits are started, we start from the estimates for the full data (a
    # "warm start", as in fastmixture): they are close to the optimum without
    # the fold, so each fold converges in a few iterations, and we skip the EM
    # priming (both init_q and init_p are given). The convergence criteria
    # are the same as for the full data.
    train <- data
    train$geno_ptr <- openadmixture_geno(data, cv = list(
      mode = "train", fold = fold, n_folds = n_folds, seed = cv_seed
    ))
    fold_fit <- openadmixture_fit(train,
      k = k, init_q = fit$q, init_p = fit$p,
      max_iter = max_iter, rtol = rtol, n_qn = n_qn, n_cores = n_cores
    )
    # Score the predictions of the fold estimates on the masked genotypes:
    # the held-out view shows only the genotypes of this fold, and the
    # deviance kernel sums d(g, mu) over them (src/openadmixture_kernels.cpp)
    held_out_ptr <- openadmixture_geno(data, cv = list(
      mode = "held_out", fold = fold, n_folds = n_folds, seed = cv_seed
    ))
    dev <- openadmixture_deviance_cpp(
      held_out_ptr, t(fold_fit$q), t(fold_fit$p), n_cores
    )
    fold_dev[fold] <- dev[1]
    fold_n[fold] <- dev[2]
    log_lines <- c(log_lines, paste0(
      "CV fold ", fold, ": ", fold_n[fold], " genotypes, ",
      fold_fit$n_iter, " iterations, error = ",
      format(fold_dev[fold] / fold_n[fold], digits = 6)
    ))
  }
  # ADMIXTURE averages the deviance "across all masked entries over all
  # folds": pool the folds rather than averaging the fold means (the two are
  # almost identical, as folds are of roughly equal size)
  cv <- sum(fold_dev) / sum(fold_n)
  log_lines <- c(
    log_lines,
    paste0("CV error (K=", k, "): ", formatC(cv, format = "f", digits = 5))
  )
  list(
    cv = cv,
    fold_cv = fold_dev / fold_n,
    fold_n = fold_n,
    log = log_lines
  )
}

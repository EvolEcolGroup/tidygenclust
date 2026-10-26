#' Prepare a gen_tibble for the OpenADMIXTURE kernels
#'
#' Extracts the FBM, the row (individual) and column (locus) indices of the
#' current selection, the lookup table used to decode the raw bytes, and the
#' number of trials (allele copies) per individual.
#'
#' @param x a [`tidypopgen::gen_tibble`]
#' @param use_imputed boolean; if `TRUE`, imputed genotypes (bytes 4-6 in
#'   bigsnpr's `CODE_IMPUTE_PRED`) are used as observed genotypes; otherwise
#'   they are treated as missing.
#' @returns a list with elements `fbm`, `rows`, `cols`, `lut`, `n_trials`,
#'   `scale`, `n_indiv` and `n_loci`
#' @keywords internal
#' @noRd
openadmixture_prepare <- function(x, use_imputed = FALSE) {
  # Everything the C++ code needs is gathered here, on the R side, so that the
  # kernels only deal with plain vectors and the FBM (the C++ functions build
  # their view of the genotypes from this list, see
  # src/openadmixture_geno_access.h).
  if (!inherits(x, "gen_tbl")) {
    stop("'x' should be a gen_tibble", call. = FALSE)
  }
  # an empty gen_tibble would reach the C++ code, which needs at least one
  # individual and one locus
  if (nrow(x) == 0 || nrow(tidypopgen::show_loci(x)) == 0) {
    stop("'x' should have at least one individual and one locus",
      call. = FALSE
    )
  }
  # The likelihood is binomial with n_i = 2 (diploids) or n_i = 1
  # (pseudohaploids). tidypopgen codes a gen_tibble's ploidy as 2 (diploid),
  # -2 (pseudohaploid, i.e. a mix of diploids and pseudohaploids), 0 (mixed
  # ploidy) or k (uniform polyploid); only the first two are supported.
  ploidy <- tidypopgen::show_ploidy(x)
  if (!ploidy %in% c(2, -2)) {
    stop(
      "gt_openadmixture() only works with diploid or pseudohaploid data ",
      "(ploidy 2 or -2); this gen_tibble has ploidy ", ploidy,
      call. = FALSE
    )
  }
  if (use_imputed && !tidypopgen::gt_has_imputed(x)) {
    stop("'use_imputed = TRUE' but this gen_tibble has not been imputed",
      call. = FALSE
    )
  }

  # The genotypes column of a gen_tibble is a vctrs vector whose data are the
  # (1-based) rows of the FBM for each individual, and which carries the FBM
  # as an attribute. The loci table records the FBM column of each locus in
  # `big_index`. Filtering, reordering and select_loci() only change these
  # indices, never the backing file, so we pass them on to C++; they can be
  # in any order.
  # `rows` follows the rows of x (hence the ids and the rows of Q), and
  # `cols` the loci as listed by show_loci() (hence the rows of P).
  fbm <- attr(x$genotypes, "fbm")
  rows <- vctrs::vec_data(x$genotypes)
  cols <- tidypopgen::show_loci(x)$big_index

  # decode raw bytes ourselves rather than relying on the live code256
  # Bytes 0, 1 and 2 are the counts of the alternate allele (R indices 1:3
  # correspond to bytes 0:2). Imputed genotypes are stored as bytes 4, 5 and
  # 6 (bigsnpr's CODE_IMPUTE_PRED), and are only used if requested. Every
  # other byte (including 3, the code for missing data) maps to NA, which the
  # C++ code treats as missing.
  lut <- rep(NA_real_, 256)
  lut[1:3] <- 0:2
  if (use_imputed) {
    lut[5:7] <- 0:2
  }

  # pseudohaploid individuals are stored as 0/2 but carry a single allele
  # indiv_ploidy() gives one value per row of x: 2 for every individual of a
  # diploid gen_tibble, and 1 or 2 (looked up by FBM row) for a
  # pseudohaploid one. The check below is a safeguard: given the check of
  # show_ploidy() above, other values should not occur.
  indiv_ploidy <- tidypopgen::indiv_ploidy(x)
  if (!all(indiv_ploidy %in% c(1, 2))) {
    stop("individual ploidies should be 1 (pseudohaploid) or 2",
      call. = FALSE
    )
  }
  # n_trials is the number of allele copies of each individual (n_i in the
  # binomial likelihood); scale converts the stored dosage into a count out
  # of n_i (a stored 2 becomes 1 for a pseudohaploid individual)
  n_trials <- as.numeric(indiv_ploidy)
  scale <- n_trials / 2

  # The C++ code reads the elements it needs by name; n_indiv and n_loci are
  # only used in R (the C++ code takes the dimensions from rows and cols).
  list(
    fbm = fbm,
    rows = as.integer(rows),
    cols = as.integer(cols),
    lut = lut,
    n_trials = n_trials,
    scale = scale,
    n_indiv = length(rows),
    n_loci = length(cols)
  )
}

#' Decode the genotypes selected in a prepared dataset into a dense matrix
#'
#' Mostly useful for testing.
#' @param data the output of `openadmixture_prepare()`
#' @param cv an optional cross-validation mask, see `openadmixture_set_cv()`
#' @returns a numeric matrix (individuals x loci), `NA` for missing (or
#'   masked) values
#' @keywords internal
#' @noRd
openadmixture_decode <- function(data, cv = NULL) {
  openadmixture_decode_cpp(openadmixture_set_cv(data, cv))
}

#' Add a cross-validation mask to a prepared dataset
#'
#' The C++ functions read the genotypes of a prepared dataset through a view
#' (src/openadmixture_geno_access.h), which they build from the dataset at
#' the start of each call. The view can hide one fold of the observed
#' genotypes (to fit the model without them, `mode = "train"`) or show only
#' that fold (to score the predictions of that model, `mode = "held_out"`).
#' Which fold each genotype belongs to is computed in C++ from a hash of
#' `seed` and the position of the genotype, so the same seed always gives the
#' same folds.
#' @param data the output of `openadmixture_prepare()`
#' @param cv `NULL` (no mask) or a list with elements `mode` ("train" or
#'   "held_out"), `fold` (1 to `n_folds`), `n_folds` and `seed` (a
#'   non-negative integer, stored as a double)
#' @returns `data`, with an element `cv` holding the integer `mode`
#'   (0 = none, 1 = train, 2 = held out) and the other elements, as expected
#'   by the C++ functions
#' @keywords internal
#' @noRd
openadmixture_set_cv <- function(data, cv = NULL) {
  # `data` is a plain list, so the mask is added to a copy (the vectors and
  # the FBM are shared, not duplicated): the caller's dataset is unchanged.
  if (is.null(cv)) {
    # mode 0: no mask (this also removes a mask that `data` already had);
    # the C++ code then ignores fold, n_folds and seed
    data$cv <- list(mode = 0L, fold = 1L, n_folds = 1L, seed = 0)
    return(data)
  }
  # The C++ code checks the seed too, but is.finite() also rejects Inf,
  # which it would let through before casting the seed to a 64-bit integer.
  # (fold and n_folds are checked in C++.)
  # (below 2^53, the largest integer that a double holds exactly; the C++
  # code converts it to a 64-bit integer)
  seed_ok <- length(cv$seed) == 1 && is.finite(cv$seed) && cv$seed >= 0 &&
    cv$seed < 2^53 && cv$seed == floor(cv$seed)
  if (!seed_ok) {
    stop("'cv_seed' must be a non-negative integer", call. = FALSE)
  }
  data$cv <- list(
    # 1 = train, 2 = held out; any other mode gives NA, rejected in C++
    mode = match(cv$mode, c("train", "held_out")),
    fold = as.integer(cv$fold),
    n_folds = as.integer(cv$n_folds),
    seed = as.numeric(cv$seed)
  )
  data
}

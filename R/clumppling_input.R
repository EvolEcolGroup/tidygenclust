# Input handling for the native clumppling port: validation and rescaling of
# the Q matrices of a gt_admix object, and grouping of replicates by K,
# following clumppling/parseInput/utils.py and _prepare_q_list() in
# data-raw/clumppling/py_rclumppling.py (the python wrapper used by previous
# versions of tidygenclust). Q files are read by tidypopgen.

#' Prepare the Q matrices of a gt_admix object
#'
#' Each matrix is validated and rescaled, matrices with K = 1 are dropped, and
#' the remaining ones are grouped by K and named as the processed files of
#' clumppling (e.g. "3_K2R1" is the first replicate with K = 2, and the third
#' matrix overall once sorted by K).
#' @param q_list a list of Q matrices (individuals x K), e.g. the element `Q`
#'   of a `gt_admix` object
#' @param tolerance tolerance used to check that rows sum to 1
#' @returns a list with elements `q_names`, `q_mats` (both sorted by K),
#'   `k_range` (the sorted unique values of K) and `k` (the K of each matrix)
#' @keywords internal
#' @noRd
clumppling_prepare_q_list <- function(q_list, tolerance = 1e-6) {
  q_list <- lapply(seq_along(q_list), function(i) {
    q <- clumppling_as_numeric_matrix(q_list[[i]])
    tryCatch(
      clumppling_validate_q(q, tolerance),
      error = function(e) {
        stop("Q matrix ", i, ": ", conditionMessage(e), call. = FALSE)
      }
    )
    clumppling_rescale_q(q)
  })
  k <- vapply(q_list, ncol, integer(1))
  keep <- k > 1
  for (i in which(!keep)) {
    message("Ignored Q matrix ", i, " with K=1.")
  }
  n_rows <- vapply(q_list, nrow, integer(1))
  expected_nrows <- n_rows[keep][1]
  bad_rows <- which(keep & n_rows != expected_nrows)
  if (length(bad_rows) > 0) {
    stop(
      "Q matrix ", bad_rows[1], " has ", n_rows[bad_rows[1]],
      " individuals, expected ", expected_nrows, "."
    )
  }
  # group by K; order() is stable, so replicates keep their order within
  # each K
  k_order <- order(k[keep], method = "radix")
  q_mats <- q_list[keep][k_order]
  k <- k[keep][k_order]
  replicate <- stats::ave(k, k, FUN = seq_along)
  q_names <- paste0(seq_along(k), "_K", k, "R", replicate)
  list(
    q_names = q_names,
    q_mats = q_mats,
    k_range = unique(k),
    k = k
  )
}

#' Check that a matrix is a valid membership matrix
#'
#' As `validate_membership_matrix()` in clumppling: all values must be in
#' \[0, 1\], and rows must sum to 1. numpy's `allclose()`, used by
#' clumppling, adds a relative tolerance of 1e-5 to `tolerance`.
#' @param q a numeric matrix
#' @param tolerance absolute tolerance for the row sums
#' @returns TRUE (invisibly), or an error
#' @keywords internal
#' @noRd
clumppling_validate_q <- function(q, tolerance = 1e-6) {
  if (anyNA(q) || any(q < 0 | q > 1)) {
    stop("Matrix contains values outside [0, 1].", call. = FALSE)
  }
  row_sums <- rowSums(q)
  bad_rows <- which(!(abs(row_sums - 1) <= tolerance + 1e-5))
  if (length(bad_rows) > 0) {
    stop(
      "Rows ", paste(utils::head(bad_rows, 10), collapse = ", "),
      if (length(bad_rows) > 10) ", ...",
      " do not sum to 1 (+/-", tolerance, ").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Rescale a membership matrix so that rows sum to 1
#'
#' As `rescale_membership_matrix()` in clumppling, a small value is added to
#' all entries to avoid zeros.
#' @param q a numeric matrix
#' @param tolerance the value added to all entries
#' @returns the rescaled matrix
#' @keywords internal
#' @noRd
clumppling_rescale_q <- function(q, tolerance = 1e-9) {
  q <- q + tolerance
  q / rowSums(q)
}

#' Coerce a Q matrix (e.g. a `q_matrix` or a data.frame) to a plain numeric
#' matrix
#' @param q a matrix-like object
#' @returns a numeric matrix without attributes other than dimensions
#' @keywords internal
#' @noRd
clumppling_as_numeric_matrix <- function(q) {
  q <- as.matrix(q)
  matrix(as.numeric(q), nrow = nrow(q), ncol = ncol(q))
}

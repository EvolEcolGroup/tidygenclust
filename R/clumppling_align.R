# Helpers to align and compare membership matrices for the native clumppling
# port (clumppling/core.py and clumppling/utils.py). The optimal alignments
# themselves are computed in C++ (src/clumppling_assignment.cpp).

#' Align all replicates with the same K
#'
#' As `align_within_k()` in clumppling: every pair of replicates is aligned,
#' and the cost of the optimal alignment is recorded.
#' @param q_list a list of Q matrices with the same dimensions
#' @returns a list with the symmetric `cost` matrix and a K x R x R array
#'   `pattern`, where `pattern[, i, j]` maps the columns of replicate j onto
#'   those of replicate i
#' @keywords internal
#' @noRd
clumppling_align_within_k <- function(q_list) {
  clumppling_align_within_k_cpp(q_list)
}

#' Reorder (and possibly merge) the columns of Q to match those of P
#'
#' As `alignQ_wrtP()` in clumppling. With `merge = TRUE`, columns of Q mapped
#' to the same column of P are summed, giving a matrix with the dimensions of
#' P. With `merge = FALSE`, among the columns of Q mapped to the same column
#' of P only the closest one is placed there, and the others are moved to new
#' columns after those of P.
#' @param p the reference matrix (N x K_P)
#' @param q the matrix to align (N x K_Q, with K_Q >= K_P)
#' @param idx for each column of Q, the column of P it maps to (1-based)
#' @param merge boolean, whether to merge columns of Q
#' @returns a list with the aligned matrix `q`, and the `pattern` giving the
#'   column of the aligned matrix for each column of Q
#' @keywords internal
#' @noRd
clumppling_align_q_wrt_p <- function(p, q, idx, merge = TRUE) {
  kp <- ncol(p)
  if (merge) {
    return(list(q = clumppling_move_columns(q, idx, kp), pattern = idx))
  }
  # the Euclidean distance of each column of Q to the column of P it maps to
  dist <- sqrt(colSums((q - p[, idx, drop = FALSE])^2))
  # among the columns of Q mapped to the same column of P, the closest one
  # (the first one on ties, as order() is stable) is placed there; the others
  # go to new columns kp + 1, kp + 2, ... in their order in Q. When every
  # column of P receives a column of Q (as in optimal alignments), the
  # pattern is a permutation of 1:kq.
  by_dist <- order(idx, dist)
  closest <- by_dist[!duplicated(idx[by_dist])]
  others <- setdiff(seq_len(ncol(q)), closest)
  pattern <- integer(ncol(q))
  pattern[closest] <- idx[closest]
  pattern[others] <- kp + seq_along(others)
  list(q = clumppling_move_columns(q, pattern), pattern = pattern)
}

#' Move the columns of Q according to a pattern
#'
#' Column j of Q is added to column `pattern[j]` of the result (as done in
#' `reorderQ_across_k()` in clumppling), so that columns with the same
#' destination are summed.
#' @param q a matrix
#' @param pattern the destination column of each column of Q
#' @param k the number of columns of the result
#' @returns a matrix with the same number of rows as Q, and `k` columns
#' @keywords internal
#' @noRd
clumppling_move_columns <- function(q, pattern, k = ncol(q)) {
  aligned <- matrix(0, nrow(q), k)
  for (j in seq_len(ncol(q))) {
    aligned[, pattern[j]] <- aligned[, pattern[j]] + q[, j]
  }
  aligned
}

#' Cost between two membership matrices of the same size
#'
#' As `cost_membership()` in clumppling: the squared Frobenius norm of the
#' difference, divided by twice the number of individuals.
#' @param p,q matrices with the same dimensions
#' @returns the cost
#' @keywords internal
#' @noRd
clumppling_cost_membership <- function(p, q) {
  if (!identical(dim(p), dim(q))) {
    stop("P and Q must have the same size.")
  }
  # numpy computes the Frobenius norm, which is then squared
  norm <- sqrt(sum((p - q)^2))
  norm^2 / (2 * nrow(p))
}

#' Cost of the clusters of P matched to a single cluster of Q
#'
#' As `cost_membership_sep()` in clumppling: the sum of the squared distances
#' between the columns of P matched to exactly one column of Q and that
#' column, divided by twice the number of individuals.
#' @param p the reference matrix
#' @param q the aligned matrix
#' @param idx for each column of Q, the column of P it maps to (1-based)
#' @returns the cost
#' @keywords internal
#' @noRd
clumppling_cost_membership_sep <- function(p, q, idx) {
  cost <- 0
  for (pj in seq_len(ncol(p))) {
    q_idx <- which(idx == pj)
    if (length(q_idx) == 1) {
      cost <- cost + sum((p[, pj] - q[, q_idx])^2)
    }
  }
  cost / (2 * nrow(p))
}

#' Convert a cost into a similarity
#'
#' As `C2Gprime()` in clumppling.
#' @param cost a numeric vector or matrix of costs
#' @returns 1 - sqrt(cost)
#' @keywords internal
#' @noRd
clumppling_cost_to_similarity <- function(cost) {
  1 - sqrt(cost)
}

#' Format an alignment pattern as a string
#'
#' As `pattern_to_str()` in clumppling (e.g. "2 1 3").
#' @param pattern an integer vector (1-based)
#' @returns a string
#' @keywords internal
#' @noRd
clumppling_pattern_to_str <- function(pattern) {
  paste(pattern, collapse = " ")
}

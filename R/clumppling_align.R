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
  kq <- ncol(q)
  if (merge) {
    aligned <- matrix(0, nrow(q), kp)
    for (j in seq_len(kq)) {
      aligned[, idx[j]] <- aligned[, idx[j]] + q[, j]
    }
    return(list(q = aligned, pattern = idx))
  }
  # columns of P that receive several columns of Q
  dups <- which(tabulate(idx, nbins = kp) > 1)
  aligned <- matrix(0, nrow(q), kq)
  new_pattern <- integer(kq)
  # for each such column of P, the closest column of Q (the first one on
  # ties) and its Euclidean distance
  dups_min_diff <- rep(Inf, kp)
  dups_min_q <- rep(NA_integer_, kp)
  # first pass: columns of Q that are the only ones matched to their column
  # of P go there; find the closest column for the other columns of P
  for (j in seq_len(kq)) {
    pj <- idx[j]
    if (!pj %in% dups) {
      new_pattern[j] <- pj
      aligned[, pj] <- q[, j]
    } else {
      diff <- sqrt(sum((q[, j] - p[, pj])^2))
      if (dups_min_diff[pj] > diff) {
        dups_min_diff[pj] <- diff
        dups_min_q[pj] <- j
      }
    }
  }
  # second pass: columns that are not the closest to their match go to new
  # columns kp + 1, kp + 2, ... in their order in Q. When every column of P
  # receives a column of Q (as in optimal alignments), the pattern is a
  # permutation of 1:kq.
  extras <- integer(0)
  for (j in seq_len(kq)) {
    pj <- idx[j]
    if (pj %in% dups) {
      if (j == dups_min_q[pj]) {
        new_pattern[j] <- pj
        aligned[, pj] <- q[, j]
      } else {
        extras <- c(extras, j)
        new_pattern[j] <- kp + length(extras)
      }
    }
  }
  for (i in seq_along(extras)) {
    aligned[, kp + i] <- q[, extras[i]]
  }
  list(q = aligned, pattern = new_pattern)
}

#' Move the columns of Q according to a pattern
#'
#' Column j of Q is added to column `pattern[j]` of the result (as done in
#' `reorderQ_across_k()` in clumppling).
#' @param q a matrix
#' @param pattern the destination column of each column of Q
#' @returns a matrix with the same dimensions as Q
#' @keywords internal
#' @noRd
clumppling_move_columns <- function(q, pattern) {
  aligned <- matrix(0, nrow(q), ncol(q))
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

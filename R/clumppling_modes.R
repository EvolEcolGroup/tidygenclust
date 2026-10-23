# Extraction of modes (groups of similar replicates with the same K) for the
# native clumppling port (community_labels_to_modes() in
# clumppling/detectMode/utils.py, find_repr_modes() and
# extract_modes_and_stats() in clumppling/core.py).

#' Convert community labels into modes
#'
#' @param communities 0-based community labels, consecutive from 0
#' @returns a list with the (1-based) indices of the replicates in each mode
#' @keywords internal
#' @noRd
clumppling_labels_to_modes <- function(communities) {
  if (max(communities) != length(unique(communities)) - 1) {
    stop("Community labels should be zero-indexed and consecutive.")
  }
  lapply(seq_len(max(communities) + 1) - 1, function(i) {
    which(communities == i)
  })
}

#' Find the representative replicate of each mode
#'
#' The representative is the replicate with the smallest total cost to the
#' other replicates in the mode (the first one in case of ties).
#' @param modes a list of replicate indices
#' @param cost_mat the cost matrix of the replicates
#' @returns the index of the representative replicate of each mode
#' @keywords internal
#' @noRd
clumppling_find_repr_modes <- function(modes, cost_mat) {
  vapply(modes, function(idx) {
    if (length(idx) == 1) {
      return(idx)
    }
    sub_cost <- cost_mat[idx, idx, drop = FALSE]
    # column sums accumulated in double precision, row by row, as numpy does
    col_sums <- sub_cost[1, ]
    for (i in seq_len(nrow(sub_cost))[-1]) {
      col_sums <- col_sums + sub_cost[i, ]
    }
    idx[which.min(col_sums)]
  }, integer(1))
}

#' Extract modes and their statistics
#'
#' As `extract_modes_and_stats()` in clumppling: for each mode, the Q matrix
#' of its representative replicate, the average of the Q matrices of its
#' replicates (once aligned to the representative), and summary statistics
#' (size, average cost and average similarity between pairs of replicates).
#' @param modes a list of replicate indices (see
#'   `clumppling_labels_to_modes()`)
#' @param cost_mat the cost matrix of the replicates
#' @param q_names the names of the replicates
#' @param q_list the Q matrices of the replicates
#' @param pattern the alignment patterns between replicates (see
#'   `clumppling_align_within_k()`)
#' @param label_prefix prefix for the mode labels (e.g. "K3M")
#' @returns a list with elements `modes`, `repr_modes`, `mode_stats` (a
#'   data.frame), `align_info` (a data.frame with the alignment of each
#'   replicate to the representative of its mode), `rep_q_modes` and
#'   `avg_q_modes`
#' @keywords internal
#' @noRd
clumppling_extract_modes <- function(modes, cost_mat, q_names, q_list,
                                     pattern, label_prefix = "M") {
  repr_modes <- clumppling_find_repr_modes(modes, cost_mat)
  mode_labels <- paste0(label_prefix, seq_along(modes))
  rep_q_modes <- list()
  avg_q_modes <- list()
  mode_stats <- list()
  align_info <- list()
  for (mode_idx in seq_along(modes)) {
    mode_label <- mode_labels[mode_idx]
    all_indices <- modes[[mode_idx]]
    i_ref <- repr_modes[mode_idx]
    repr_name <- q_names[i_ref]
    q_ref <- q_list[[i_ref]]
    identity_pattern <- clumppling_pattern_to_str(seq_len(ncol(q_ref)))
    rep_q_modes[[mode_idx]] <- q_ref
    if (length(all_indices) == 1) {
      mode_stats[[mode_idx]] <- data.frame(
        Mode = mode_label, Representative = repr_name, Size = 1, Cost = 0,
        Performance = 1
      )
      avg_q_modes[[mode_idx]] <- q_ref
      align_info[[mode_idx]] <- data.frame(
        Mode = mode_label, Representative = repr_name,
        Replicate = repr_name, Alignment = identity_pattern
      )
      next
    }
    sub_cost <- cost_mat[all_indices, all_indices, drop = FALSE]
    mode_stats[[mode_idx]] <- data.frame(
      Mode = mode_label, Representative = repr_name,
      Size = as.numeric(length(all_indices)),
      Cost = clumppling_avg_tril(sub_cost),
      Performance = clumppling_avg_tril(clumppling_cost_to_similarity(sub_cost))
    )
    # align the other replicates to the representative, and average them
    q_sum <- matrix(0, nrow(q_ref), ncol(q_ref))
    alignments <- character(length(all_indices))
    for (i in seq_along(all_indices)) {
      r <- all_indices[i]
      if (r == i_ref) {
        q_sum <- q_sum + q_ref
        alignments[i] <- identity_pattern
      } else {
        idx <- pattern[, i_ref, r]
        q_sum <- q_sum + clumppling_align_q_wrt_p(q_ref, q_list[[r]], idx)$q
        alignments[i] <- clumppling_pattern_to_str(idx)
      }
    }
    avg_q_modes[[mode_idx]] <- q_sum / length(all_indices)
    align_info[[mode_idx]] <- data.frame(
      Mode = mode_label, Representative = repr_name,
      Replicate = q_names[all_indices], Alignment = alignments
    )
  }
  list(
    modes = modes,
    repr_modes = repr_modes,
    mode_stats = do.call(rbind, mode_stats),
    align_info = do.call(rbind, align_info),
    rep_q_modes = rep_q_modes,
    avg_q_modes = avg_q_modes
  )
}

#' Average of the lower triangle of a square matrix (excluding the diagonal)
#'
#' @param a a square matrix
#' @returns the average
#' @keywords internal
#' @noRd
clumppling_avg_tril <- function(a) {
  n <- nrow(a)
  sum(a[lower.tri(a)]) / (n * (n - 1) / 2)
}

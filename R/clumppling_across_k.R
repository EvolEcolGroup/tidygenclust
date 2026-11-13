# Alignment of modes across different values of K for the native clumppling
# port (align_across_k() and reorderQ_across_k() in clumppling/core.py).
#
# Conventions: modes are named "K<K>M<m>" (e.g. "K3M2"), and a pair of modes
# is labelled "<mode P>-<mode Q>", where P is the mode with fewer (or as
# many) clusters; the alignment of the pair is an integer vector `idx` that
# gives, for each column of Q, the (1-based) column of P it is matched to.
# Every column of P receives at least one column of Q: with the same K the
# alignment is a permutation, otherwise some columns of P receive several
# columns of Q (e.g. a cluster of K = 2 that splits into two clusters of
# K = 3).

#' Align modes across K
#'
#' As `align_across_k()` in clumppling: every mode is aligned to every mode
#' of the next smaller K (merging two clusters if `merge` is TRUE and the two
#' K differ by one), and to every mode of the same K.
#' @param k_range the sorted values of K
#' @param q_modes_list a list (one element per K) of lists of mode matrices
#' @param mode_names_list a list (one element per K) of mode names
#' @param merge boolean, whether to merge two clusters when aligning K + 1 to K
#' @returns a list with elements `alignment_acrossK` and `cost_acrossK`
#'   (named lists, with names "smaller K mode-larger K mode"), and
#'   `best_acrossK` and `major_acrossK` (data.frames summarising, for each
#'   pair of consecutive K, the best aligned pair of modes and the pair of
#'   major modes)
#' @keywords internal
#' @noRd
clumppling_align_across_k <- function(k_range, q_modes_list, mode_names_list,
                                      merge = TRUE) {
  stopifnot(
    length(k_range) == length(q_modes_list),
    length(k_range) == length(mode_names_list)
  )
  # from large K to small K
  k_sorted <- sort(k_range, decreasing = TRUE)
  n_k <- length(k_sorted)
  # the pairs of K to align: each K with the next smaller one, then each K
  # with itself (in this order, which is that of the names of the results)
  k_comb <- c(
    lapply(seq_len(n_k - 1), function(i) k_sorted[c(i, i + 1)]),
    lapply(seq_len(n_k), function(i) k_sorted[c(i, i)])
  )

  alignment_across_k <- list()
  cost_across_k <- list()
  # for each pair of K, the indices (ri, rj) of the modes with the smallest
  # alignment objective (the first such pair on ties)
  best_pairs <- list()
  for (k_pair in k_comb) {
    k1 <- k_pair[1]
    k2 <- k_pair[2]
    i_k1 <- match(k1, k_range)
    i_k2 <- match(k2, k_range)
    best_obj <- Inf
    best_pair <- c(1L, 1L)
    for (ri in seq_along(q_modes_list[[i_k1]])) {
      for (rj in seq_along(q_modes_list[[i_k2]])) {
        q1 <- q_modes_list[[i_k1]][[ri]] # the one with more clusters
        q2 <- q_modes_list[[i_k2]][[rj]]
        # q1 (larger K) is aligned to q2: res$idx maps the columns of q1
        # onto those of q2. With `merge`, when K differs by one, two
        # clusters of q1 are merged and the others matched one to one; a
        # mode is aligned to itself with the identity.
        if (k1 == k2 && ri == rj) {
          res <- list(obj = 0, idx = seq_len(k1))
        } else if (merge && k1 - k2 == 1) {
          res <- clumppling_align_ilp_diff1_cpp(q2, q1)
        } else {
          res <- clumppling_align_ilp_cpp(q2, q1)
        }
        if (best_obj > res$obj) {
          best_obj <- res$obj
          best_pair <- c(ri, rj)
        }
        # the cost of the pair is that of q2 against q1 with its matched
        # columns summed (whatever `merge`), so it measures how well the
        # clusters of the larger K nest within those of the smaller K
        aligned <- clumppling_align_q_wrt_p(q2, q1, res$idx, merge = TRUE)$q
        pair_label <- paste0(
          mode_names_list[[i_k2]][rj], "-", mode_names_list[[i_k1]][ri]
        )
        cost_across_k[[pair_label]] <- clumppling_cost_membership(q2, aligned)
        alignment_across_k[[pair_label]] <- res$idx
      }
    }
    best_pairs[[paste(k1, k2)]] <- best_pair
  }

  # summarise the best and the major (first mode) pairs between consecutive K
  # (from the largest K); the separate-cluster cost only counts the clusters
  # of the smaller K that are matched to a single cluster of the larger K
  summarise_pairs <- function(get_pair, pair_col) {
    rows <- lapply(seq_len(n_k - 1), function(i) {
      i_k1 <- match(k_sorted[i], k_range)
      i_k2 <- match(k_sorted[i + 1], k_range)
      pair <- get_pair(k_sorted[i], k_sorted[i + 1])
      pair_label <- paste0(
        mode_names_list[[i_k2]][pair[2]], "-", mode_names_list[[i_k1]][pair[1]]
      )
      cost_pair <- cost_across_k[[pair_label]]
      cost_sep <- clumppling_cost_membership_sep(
        q_modes_list[[i_k2]][[pair[2]]], q_modes_list[[i_k1]][[pair[1]]],
        alignment_across_k[[pair_label]]
      )
      data.frame(
        pair = pair_label,
        Cost = cost_pair,
        Performance = clumppling_cost_to_similarity(cost_pair),
        sep_cost = cost_sep,
        sep_perf = clumppling_cost_to_similarity(cost_sep)
      )
    })
    out <- do.call(rbind, rows)
    names(out) <- c(
      pair_col, "Cost", "Performance", "Separate-Cluster Cost",
      "Separate-Cluster Performance"
    )
    out
  }
  best_across_k <- summarise_pairs(
    function(k1, k2) best_pairs[[paste(k1, k2)]], "Best Pair"
  )
  major_across_k <- summarise_pairs(function(k1, k2) c(1L, 1L), "Major Pair")

  list(
    alignment_acrossK = alignment_across_k,
    cost_acrossK = cost_across_k,
    best_acrossK = best_across_k,
    major_acrossK = major_across_k
  )
}

#' Reorder the clusters of all modes so that they are aligned across K
#'
#' As `reorderQ_across_k()` in clumppling: starting from the anchor mode of
#' the smallest K, the anchor pairs (one per pair of consecutive K) are used
#' to propagate the order of the clusters to larger K, and the other modes of
#' each K are aligned to the anchor mode of that K. New clusters appearing at
#' larger K are placed after the existing ones.
#' @param k_range the sorted values of K (at least two)
#' @param q_modes_list a list (one element per K) of lists of mode matrices
#' @param mode_names_list a list (one element per K) of mode names
#' @param alignment_across_k the alignments returned by
#'   `clumppling_align_across_k()`
#' @param anchor_pairs the labels of the anchor pairs, from the largest to the
#'   smallest K (e.g. the best pairs)
#' @returns a named list of aligned mode matrices
#' @keywords internal
#' @noRd
clumppling_reorder_q_across_k <- function(k_range, q_modes_list,
                                          mode_names_list, alignment_across_k,
                                          anchor_pairs) {
  stopifnot(length(anchor_pairs) == length(k_range) - 1)
  # the two modes of each anchor pair, from the smallest K: the i-th pair
  # links a mode of k_range[i] to a mode of k_range[i + 1]
  anchor_modes <- strsplit(rev(anchor_pairs), "-", fixed = TRUE)
  get_mode <- function(i, mode) {
    q_modes_list[[i]][[match(mode, mode_names_list[[i]])]]
  }
  # the pattern of a mode q aligned (without merging) to a mode p whose
  # clusters go to the positions given by p_pattern
  follow <- function(p_pattern, p, q, pair) {
    aligned <- clumppling_align_q_wrt_p(p, q, alignment_across_k[[pair]],
      merge = FALSE
    )
    p_pattern[aligned$pattern]
  }
  # for each mode, the final (1-based) position of each of its columns, so
  # that column j of the mode becomes column patterns[[mode]][j] of its
  # aligned matrix. Positions are consistent across K: a cluster matched to
  # a cluster of the smaller K takes its position, and the new clusters of a
  # larger K go after those of the smaller K.
  patterns <- list()
  aligned_qs <- list()
  for (i in seq_along(k_range)) {
    if (i == 1) {
      # the anchor mode of the smallest K keeps its order
      anchor <- anchor_modes[[1]][1]
      patterns[[anchor]] <- seq_len(k_range[1])
    } else {
      # the anchor mode of this K follows the mode of the previous K in its
      # anchor pair (not necessarily the anchor of that K, but all the modes
      # of that K already have a pattern): the closest match of each cluster
      # of that mode takes its position, and the other clusters of this K
      # (the new ones) go after the clusters of the previous K
      prev <- anchor_modes[[i - 1]][1]
      anchor <- anchor_modes[[i - 1]][2]
      new_clusters <- seq.int(k_range[i - 1] + 1, k_range[i])
      patterns[[anchor]] <- follow(
        c(patterns[[prev]], new_clusters), get_mode(i - 1, prev),
        get_mode(i, anchor), paste0(prev, "-", anchor)
      )
    }
    # the other modes of this K are aligned to its anchor mode (with the same
    # K the alignments are permutations), and then follow its order
    q_anchor <- get_mode(i, anchor)
    for (mode in c(anchor, setdiff(mode_names_list[[i]], anchor))) {
      q <- get_mode(i, mode)
      if (mode != anchor) {
        patterns[[mode]] <- follow(
          patterns[[anchor]], q_anchor, q, paste0(anchor, "-", mode)
        )
      }
      aligned_qs[[mode]] <- clumppling_move_columns(q, patterns[[mode]])
    }
  }
  aligned_qs
}

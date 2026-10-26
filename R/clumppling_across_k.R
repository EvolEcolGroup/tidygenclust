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
#' @param k_range the sorted values of K
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
  # anchor pairs from the smallest K: the i-th pair links a mode of
  # k_range[i] to a mode of k_range[i + 1]
  anchor_pairs_rev <- rev(anchor_pairs)
  aligned_qs <- list()
  # for each mode, the final (1-based) position of each of its columns, so
  # that column j of the mode becomes column base_patterns[[mode]][j] of its
  # aligned matrix. Positions are consistent across K: a cluster matched to
  # a cluster of the smaller K takes its position, and the new clusters of a
  # larger K go after those of the smaller K.
  base_patterns <- list()
  split_pair <- function(pair) strsplit(pair, "-", fixed = TRUE)[[1]]

  # the smallest K: the anchor mode keeps its order, the other modes are
  # aligned to it (as their alignments to it are permutations, they are also
  # their base patterns)
  mode_names <- mode_names_list[[1]]
  # there is always an anchor pair, as there are at least two K; the
  # alternative is kept from clumppling
  if (length(anchor_pairs) > 0) {
    m1 <- split_pair(anchor_pairs_rev[1])[1]
    m_m1 <- match(m1, mode_names)
    base_patterns[[m1]] <- seq_len(k_range[1])
  } else {
    m_m1 <- 1L
    m1 <- mode_names[1]
  }
  aligned_qs[[m1]] <- q_modes_list[[1]][[m_m1]]
  for (m in seq_along(q_modes_list[[1]])) {
    if (m == m_m1) next
    m2 <- mode_names[m]
    ali_pat <- alignment_across_k[[paste0(m1, "-", m2)]]
    aligned_qs[[m2]] <- clumppling_move_columns(q_modes_list[[1]][[m]], ali_pat)
    base_patterns[[m2]] <- ali_pat
  }

  # propagate to larger K through the anchor pairs
  for (i_pair in seq_along(anchor_pairs_rev)) {
    # m1: the mode of k_range[i_pair] in this pair (not necessarily the
    # anchor of the previous step, but all the modes of k_range[i_pair]
    # already have a base pattern); m2: the mode of the next K (i_k), which
    # becomes the anchor of that K
    pair <- anchor_pairs_rev[i_pair]
    m1 <- split_pair(pair)[1]
    m2 <- split_pair(pair)[2]
    i_k <- i_pair + 1
    m_m2 <- match(m2, mode_names_list[[i_k]])
    q <- q_modes_list[[i_k]][[m_m2]]
    m1_k <- k_range[i_pair]
    p <- q_modes_list[[i_pair]][[match(m1, mode_names_list[[i_pair]])]]
    # without merging, each column of P keeps the closest of the columns of Q
    # matched to it (positions 1 to m1_k), and the other columns of Q get
    # positions m1_k + 1, m1_k + 2, ...
    new_pattern <- clumppling_align_q_wrt_p(p, q, alignment_across_k[[pair]],
      merge = FALSE
    )$pattern
    # clusters matched to the smaller K follow its order, new ones are added
    # at the end
    pat <- vapply(new_pattern, function(i) {
      if (i <= m1_k) base_patterns[[m1]][i] else i
    }, integer(1))
    aligned_qs[[m2]] <- clumppling_move_columns(q, pat)
    base_patterns[[m2]] <- pat

    # the other modes of this K are aligned to its anchor mode (with the same
    # K the alignments are permutations), and then follow its order
    for (m in seq_along(q_modes_list[[i_k]])) {
      if (m == m_m2) next
      m3 <- mode_names_list[[i_k]][m]
      q3 <- q_modes_list[[i_k]][[m]]
      new_pattern <- clumppling_align_q_wrt_p(q, q3,
        alignment_across_k[[paste0(m2, "-", m3)]],
        merge = FALSE
      )$pattern
      pat <- base_patterns[[m2]][new_pattern]
      aligned_qs[[m3]] <- clumppling_move_columns(q3, pat)
      base_patterns[[m3]] <- pat
    }
  }
  aligned_qs
}

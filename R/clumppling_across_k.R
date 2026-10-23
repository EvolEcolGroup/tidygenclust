# Alignment of modes across different values of K for the native clumppling
# port (align_across_k() and reorderQ_across_k() in clumppling/core.py).

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
  k_comb <- c(
    lapply(seq_len(n_k - 1), function(i) k_sorted[c(i, i + 1)]),
    lapply(seq_len(n_k), function(i) k_sorted[c(i, i)])
  )

  alignment_across_k <- list()
  cost_across_k <- list()
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
  anchor_pairs_rev <- rev(anchor_pairs)
  aligned_qs <- list()
  base_patterns <- list()
  split_pair <- function(pair) strsplit(pair, "-", fixed = TRUE)[[1]]

  # the smallest K: the anchor mode keeps its order, the other modes are
  # aligned to it
  mode_names <- mode_names_list[[1]]
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
    pair <- anchor_pairs_rev[i_pair]
    m1 <- split_pair(pair)[1]
    m2 <- split_pair(pair)[2]
    i_k <- i_pair + 1
    m_m2 <- match(m2, mode_names_list[[i_k]])
    q <- q_modes_list[[i_k]][[m_m2]]
    m1_k <- k_range[i_pair]
    p <- q_modes_list[[i_pair]][[match(m1, mode_names_list[[i_pair]])]]
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

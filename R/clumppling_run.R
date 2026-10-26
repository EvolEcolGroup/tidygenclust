#' Run the clumppling pipeline
#'
#' The native equivalent of `clumppling_run()` in
#' data-raw/clumppling/py_rclumppling.py (itself derived from clumppling's
#' `__main__.py`): align replicates within each K, detect modes, align modes
#' across K, and reorder the clusters of the modes so that they are consistent
#' across K.
#' @param q_data a list of Q matrices and their metadata, as returned by
#'   `clumppling_prepare_q_list()`
#' @param cd_method,cd_res,test_comm,comm_min,comm_max see [gt_clumppling()]
#' @param merge,use_rep,use_best_pair see [gt_clumppling()]
#' @returns a list with the results, including intermediate ones: `N`,
#'   `K_range`, `mode_names`, `aligned_modes` (a list of matrices),
#'   `cost_acrossK`, `alignment_acrossK`, `best_acrossK`, `major_acrossK`,
#'   `mode_stats`, `mode_alignment`, and, for each K, the `cost_withinK`
#'   matrices and the `communities` of replicates
#' @keywords internal
#' @noRd
clumppling_run <- function(q_data, cd_method = "louvain", cd_res = 1,
                           test_comm = FALSE, comm_min = 1e-6,
                           comm_max = 1e-2, merge = TRUE, use_rep = TRUE,
                           use_best_pair = TRUE) {
  k_range <- q_data$k_range
  if (length(k_range) < 2) {
    stop(
      "At least two K values are required for alignment across K: K=",
      paste(k_range, collapse = ", ")
    )
  }

  # align within K, and detect and extract modes; the results are lists with
  # one element per K (in increasing order), named "K2", "K3", ...
  cd_res_list <- list()
  cost_within_k <- list()
  communities_list <- list()
  for (k in k_range) {
    # the replicates with this K, in their order in the input (replicate r
    # of this K is named "<i>_K<k>R<r>", see clumppling_prepare_q_list())
    idx_k <- which(q_data$k == k)
    q_names_k <- q_data$q_names[idx_k]
    q_list_k <- q_data$q_mats[idx_k]
    # 1. optimal alignment and cost of every pair of replicates
    aligned <- clumppling_align_within_k(q_list_k)
    # 2. modes: communities in the network of replicates (0-based labels,
    # the largest community first), turned into lists of the (1-based)
    # indices of their replicates within this K
    communities <- clumppling_detect_communities(aligned$cost,
      test_comm = test_comm, method = cd_method, res = cd_res,
      min_threshold = comm_min, max_threshold = comm_max
    )
    modes <- clumppling_labels_to_modes(communities)
    # 3. representative and average Q matrix, and statistics, of each mode
    # (named "K<k>M1", "K<k>M2", ...)
    k_label <- paste0("K", k)
    cd_res_list[[k_label]] <- clumppling_extract_modes(modes,
      aligned$cost, q_names_k, q_list_k, aligned$pattern,
      label_prefix = paste0(k_label, "M")
    )
    cost_within_k[[k_label]] <- aligned$cost
    communities_list[[k_label]] <- communities
  }

  # align modes across K: each mode is summarised by the Q matrix of its
  # representative replicate (use_rep) or by the average of its replicates
  # (aligned to the representative); both lists have one element per K
  mode_names_list <- lapply(cd_res_list, function(x) x$mode_stats$Mode)
  q_modes_list <- lapply(cd_res_list, function(x) {
    if (use_rep) x$rep_q_modes else x$avg_q_modes
  })
  across_k <- clumppling_align_across_k(k_range, q_modes_list,
    mode_names_list,
    merge = merge
  )

  # reorder the Q matrices of the modes according to the alignment. The
  # anchor pairs (one per pair of consecutive K, from the largest K, e.g.
  # "K4M1-K5M2") link the modes through which the order of the clusters is
  # propagated from the smallest to the largest K: the pairs of modes with
  # the best alignment (use_best_pair), or the pairs of major modes (the
  # first, i.e. largest, mode of each K)
  anchor_pairs <- if (use_best_pair) {
    across_k$best_acrossK[["Best Pair"]]
  } else {
    across_k$major_acrossK[["Major Pair"]]
  }
  aligned_qs <- clumppling_reorder_q_across_k(k_range, q_modes_list,
    mode_names_list, across_k$alignment_acrossK, anchor_pairs
  )

  # all modes, by increasing K and then by mode number; the data.frames of
  # the different K are stacked, as pd.concat() in python
  mode_names <- unlist(mode_names_list, use.names = FALSE)
  rbind_cd <- function(element) {
    out <- do.call(rbind, unname(lapply(cd_res_list, `[[`, element)))
    rownames(out) <- NULL
    out
  }
  list(
    N = nrow(q_data$q_mats[[1]]),
    K_range = as.integer(k_range),
    mode_names = mode_names,
    aligned_modes = unname(aligned_qs[mode_names]),
    cost_acrossK = across_k$cost_acrossK,
    alignment_acrossK = across_k$alignment_acrossK,
    best_acrossK = across_k$best_acrossK,
    major_acrossK = across_k$major_acrossK,
    mode_stats = rbind_cd("mode_stats"),
    mode_alignment = rbind_cd("align_info"),
    cost_withinK = cost_within_k,
    communities = communities_list
  )
}

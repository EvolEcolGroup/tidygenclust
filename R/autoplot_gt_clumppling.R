#' autoplot for clumppling objects
#'
#' An autoplot method to generate quick visualisations for [`gt_clumppling`]
#' objects. Available types are:
#' - 'modes': all aligned modes in structure plots over a multipartite graph,
#' where better alignment between the modes is indicated by the darker color of
#' the edges connecting their structure plots (the colour scale of the costs of
#' optimal alignment is shown in an empty cell of the grid, if there is one),
#' and the number of replicates in each mode is given in brackets.
#' - 'modes_within_k': all the modes with `k` clusters, in a series of
#' structure plots.
#' - 'major_modes': the major modes of each K aligned in a series of structure
#' plots.
#' - 'all_modes': all aligned modes in a series of structure plots.
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are not
#' customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' If you would like to generate an annotated autoplot, ensure that all
#' individuals from the same population are adjacent to one another in the
#' Q-matrix or gt_admix object supplied to gt_clumppling. Autoplot 'group'
#' argument requires that all individuals from the same group are adjacent.
#'
#' @param object a [`gt_clumppling`] object
#' @param type the type of plot, one of 'modes', 'modes_within_k', 'major_modes'
#'   or 'all_modes'.
#' @param group a vector of membership to a-priori groups (e.g. populations).
#'   Note that individuals from the same group need to be adjacent to each other
#' @param k the k value to be plotted if 'type' is 'modes_within_k'
#' @param ... not used at the moment
#' @returns a plot
#' @rdname autoplot_gt_clumppling
#' @export

autoplot.gt_clumppling <- function(
  object,
  type = c("modes", "modes_within_k", "major_modes", "all_modes"),
  group = NULL,
  k = NULL,
  ...
) {
  type <- match.arg(type)
  rlang::check_dots_empty()
  # check that group, if given, is coherent with the rest of the data
  if (!is.null(group)) {
    # groups are defined by their values (as text): rle() does not accept
    # factors, which are often used for populations
    group <- as.character(group)
    # first check that labels are in blocks: if each group forms a single
    # block of consecutive individuals, the number of runs of identical
    # values equals the number of distinct values
    if (length(rle(group)$values) != length(unique(group))) {
      stop(
        paste0(
          "values in 'group' are not ordered (they should be ",
          "in consecutive blocks, one per group)"
        )
      )
    }
    if (length(group) != object$N) {
      stop(
        paste0(
          "'group' should be of the same length ",
          "as the original data (as found in object$N)"
        )
      )
    }
    # the position of the last individual of each group, with the groups in
    # their order of appearance (a named vector, with the names of the
    # groups): plot_q() draws the boundaries between groups after these
    # positions, and labels each group in the middle of its block
    group_x <- cumsum(table(fct_inorder_base(group)))
  } else {
    group_x <- NULL
  }

  # the 'modes' plot is given the original group vector, as it calls
  # autoplot() again to create the structure plots of all the modes
  if (type == "modes") {
    plot_modes(object, group)
  } else if (type == "modes_within_k") {
    plot_modes_within_k(object, group_x, k = k)
  } else if (type == "major_modes") {
    plot_major_modes(object, group_x)
  } else if (type == "all_modes") {
    plot_all_modes(object, group_x)
  }
}

# Plot major modes
plot_major_modes <- function(object, group_x) {
  if (!is.null(attr(object, "subset_indivs"))) {
    warning(
      "You are plotting a gt_clumppling object that has been subsetted by ",
      "individuals. ",
      "This function can be used to create plot insets or explore data, but ",
      "modes will not be representative of this subset. ",
      "To visualise modes within a subset of individuals, please subset your ",
      "gt_admix object and re-run gt_clumppling."
    )
  }

  # one row per major mode (the mode of rank 1 of each K)
  major_modes <- tidy(object, matrix = "major_modes")
  # label each plot with its value of K (all major modes have m == 1)
  y_labels <- rep("", nrow(major_modes))
  y_labels[major_modes$m == 1] <- paste0(
    "K = ",
    major_modes$k[major_modes$m == 1]
  )
  # one structure plot per major mode: the i-th plot uses the i-th Q matrix
  # and the i-th label
  plot_list <- lapply(
    seq_len(nrow(major_modes)),
    plot_q_from_list,
    object$aligned_modes[major_modes$label],
    y_labels,
    group_x
  )
  # stack the plots in a single column; identical x axes (the group labels,
  # if any) are collected, so that they are only shown below the last plot
  patchwork::wrap_plots(plot_list, ncol = 1) +
    patchwork::plot_layout(axes = "collect")
}

# Plot all modes
plot_all_modes <- function(object, group_x) {
  if (!is.null(attr(object, "subset_indivs"))) {
    warning(
      "You are plotting a gt_clumppling object that has been subsetted by ",
      "individuals. ",
      "This function can be used to create plot insets or explore data, but ",
      "modes will not be representative of this subset. ",
      "To visualise modes within a subset of individuals, please subset your ",
      "gt_admix object and re-run gt_clumppling."
    )
  }

  # one row per mode, in the same order as object$aligned_modes
  all_modes <- tidy(object, matrix = "modes")
  # only the first mode of each K is labelled with its value of K, so that the
  # labels mark where each block of modes with the same K starts
  y_labels <- rep("", nrow(all_modes))
  y_labels[all_modes$m == 1] <- paste0("K = ", all_modes$k[all_modes$m == 1])
  # one structure plot per mode, stacked in a single column (see
  # plot_major_modes())
  plot_list <- lapply(
    seq_len(nrow(all_modes)),
    plot_q_from_list,
    object$aligned_modes,
    y_labels,
    group_x
  )
  patchwork::wrap_plots(plot_list, ncol = 1) +
    patchwork::plot_layout(axes = "collect")
}

# Plot all the modes with a given value of K
plot_modes_within_k <- function(object, group_x, k) {
  if (!is.null(attr(object, "subset_indivs"))) {
    warning(
      "You are plotting a gt_clumppling object that has been subsetted by ",
      "individuals. ",
      "This function can be used to create plot insets or explore data, but ",
      "modes will not be representative of this subset. ",
      "To visualise modes within a subset of individuals, please subset your ",
      "gt_admix object and re-run gt_clumppling."
    )
  }

  if (is.null(k)) {
    stop("k must be provided for modes_within_k")
  }
  # a single value of K with modes (otherwise the plot would be empty, and
  # only fail when printed, or a vector would be recycled)
  if (length(k) != 1 || !k %in% object$K_range) {
    stop(
      "k must be one of the values of K in K_range (",
      paste(object$K_range, collapse = ", "), ")"
    )
  }
  # select the modes with the requested K, and label each plot with the name
  # of its mode (e.g. "K3M2")
  all_modes <- tidy(object, matrix = "modes")
  k_modes <- all_modes[all_modes$k == k, ]
  y_labels <- k_modes$label
  # one structure plot per mode, stacked in a single column (see
  # plot_major_modes())
  plot_list <- lapply(
    seq_len(nrow(k_modes)),
    plot_q_from_list,
    object$aligned_modes[k_modes$label],
    y_labels,
    group_x
  )
  patchwork::wrap_plots(plot_list, ncol = 1) +
    patchwork::plot_layout(axes = "collect")
}

# a compact q plot to use as a panel
# @param q_tidied a q values tidied with tidy
# @param y_lab the y label if present
# @param group_x the x values where groups switch, a black line will be
#   plotted here
plot_q <- function(q_tidied, y_lab = "", group_x) {
  # a structure plot: one bar per individual (at x = its index), split into
  # its membership to each cluster. Bars have width 1, so that there is no
  # space between individuals, and the clusters are stacked in reverse order
  # (the first cluster at the bottom); the distruct theme and colours are the
  # ones used by tidypopgen for gt_admix objects
  # the bars are drawn at positions 1 to N: after
  # subset_gt_clumppling(indivs = ), the ids are the original indices of the
  # individuals (e.g. 301 to 399), whereas the boundaries of the groups
  # (group_x) and the annotations of plot_modes() count from the first
  # individual of the subset (the individuals of a subset are consecutive)
  q_tidied$id <- q_tidied$id - min(q_tidied$id) + 1L
  # clusters are mapped to colours (and stacked) by their number: as text,
  # "10" would sort before "2", and the same aligned cluster would get a
  # different colour for K < 10 and K >= 10
  plt <- ggplot2::ggplot(
    q_tidied,
    ggplot2::aes(.data$id, .data$percentage,
      fill = factor(as.integer(.data$q))
    )
  ) +
    ggplot2::geom_col(
      width = 1,
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    tidypopgen::theme_distruct() +
    ggplot2::labs(y = y_lab) +
    # memberships sum to 1; drawing is not clipped to the panel, so that
    # annotations above the bars (e.g. the number of replicates added by
    # plot_modes()) are visible
    ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off") +
    tidypopgen::scale_fill_distruct()
  # add vertical lines if we have groups
  if (!is.null(group_x)) {
    # a vertical segment, spanning the whole bar height, after the last
    # individual of each group (i.e. between two bars)
    segment_data <- data.frame(
      x = group_x + 0.5,
      xend = group_x + 0.5,
      y = rep(0, length(group_x)),
      yend = rep(1, length(group_x))
    )

    plt <- plt +
      ggplot2::geom_segment(
        data = segment_data,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          xend = .data$xend,
          yend = .data$yend
        ),
        inherit.aes = FALSE
      )
    # the midpoints between consecutive values of a vector
    get_mid_points <- function(vec) {
      (vec[-length(vec)] + vec[-1L]) / 2.
    }
    # label each group on the x axis, in the middle of its block of bars
    # (between the end of the previous group, or 0, and its own end)
    plt <- plt +
      ggplot2::scale_x_continuous(
        breaks = get_mid_points(c(0, group_x)),
        labels = names(group_x)
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
      )
  } else {
    # without groups, the individual indices on the x axis are not useful
    plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_blank(), )
  }

  plt
}

# plot the i element from a q list (a helper to use plot_q() with lapply()
# over the indices of the modes, so that each plot gets the corresponding
# Q matrix and y label)
plot_q_from_list <- function(i, q_list, y_labs, group_x) {
  plot_q(tidy_q(q_list[[i]]), y_labs[i], group_x)
}

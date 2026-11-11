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
    group_x <- cumsum(table(factor(group, levels = unique(group))))
  } else {
    group_x <- NULL
  }
  if (type == "modes_within_k") {
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
  }
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

  if (type == "modes") {
    return(plot_modes(object, group_x))
  }
  # the other types are a column of structure plots, one per mode (in the
  # order of object$aligned_modes), which differ only in the modes they show
  # and in the labels of the plots: the name of each mode (e.g. "K3M2") for
  # the modes of a single K, and otherwise the value of K for the first mode
  # of each K, marking where each block of modes with the same K starts
  modes <- tidy(object, matrix = "modes")
  modes <- switch(type,
    modes_within_k = modes[modes$k == k, ],
    major_modes = modes[modes$m == 1, ],
    all_modes = modes
  )
  y_labels <- if (type == "modes_within_k") {
    modes$label
  } else {
    ifelse(modes$m == 1, paste0("K = ", modes$k), "")
  }
  plot_list <- plot_mode_list(object, modes$label, y_labels, group_x)
  # stack the plots in a single column; identical x axes (the group labels,
  # if any) are collected, so that they are only shown below the last plot
  patchwork::wrap_plots(plot_list, ncol = 1) +
    patchwork::plot_layout(axes = "collect")
}

# a structure plot (see plot_q()) for each of the modes with the given
# labels, with the corresponding y labels (recycled)
plot_mode_list <- function(object, labels, y_labels, group_x) {
  Map(function(label, y_lab) {
    plot_q(tidy_q(object$aligned_modes[[label]]), y_lab, group_x)
  }, labels, y_labels, USE.NAMES = FALSE)
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

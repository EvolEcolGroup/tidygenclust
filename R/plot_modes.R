# Plot the modes of a gt_clumppling object over a multipartite graph (the
# 'modes' type of autoplot()).
#
# The modes are arranged on a grid, with one row per value of K (in the order
# of object$K_range) and one column per mode rank (M1, M2, ...), so that the
# major modes are in the first column. The figure is a single empty ggplot
# whose x and y scales both span [0, 1], so that its coordinates match the
# "npc" units used by patchwork::inset_element(): the edges of the graph are
# drawn as lines in this ggplot, and every other element (the structure plot
# of each mode, the legend of the edge colours and the row and column labels)
# is then placed on top of it as an inset (see mode_cells() for the layout).
# @param object a gt_clumppling object
# @param group_x the position of the last individual of each group (see
#   autoplot.gt_clumppling()), or NULL
plot_modes <- function(object, group_x) {
  k_range <- object$K_range
  modes <- tidy(object, matrix = "modes")
  # the number of columns of the grid: the largest number of modes of a K
  n_col <- max(modes$m)
  cells <- mode_cells(modes$k, modes$m, k_range, n_col)

  # the edges of the graph: two rows (the two ends of the edge) for each pair
  # of modes of consecutive values of K, coloured by their alignment cost
  # (with the default colour gradient, lower costs, i.e. better alignments,
  # are darker)
  multi_plot <- ggplot2::ggplot(
    data = mode_edges(object, n_col),
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$pair)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(colour = .data$cost),
      linewidth = 2,
      lineend = "round"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "native")) +
    # fix both scales to [0, 1] with no expansion, so that the data coordinates
    # are the same as the npc coordinates of the insets
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0))

  # extract the colour legend as a grob, so that it can be placed in an empty
  # cell of the grid, and then remove it from the main plot
  legend <- .get_legend(multi_plot, position = "bottom")
  multi_plot <- multi_plot + ggplot2::guides(color = "none")

  # the structure plot of each mode, without axes and labels (the group
  # boundaries drawn by plot_q() are kept), with the number of replicates in
  # the mode above the right end of the plot (outside the panel, which is
  # drawn with clip = "off")
  mode_plots <- plot_mode_list(object, modes$label, "", group_x)
  for (i in seq_along(mode_plots)) {
    n <- length(object$mode_replicates[[modes$label[i]]])
    plt <- mode_plots[[i]] +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "native")) +
      ggplot2::annotate(
        "text",
        x = object$N * 0.95,
        y = 1.15,
        label = paste0("(", n, ")"),
        size = 3
      )
    multi_plot <- add_inset(multi_plot, plt, cells[i, ])
  }

  # if a cell of the last column is free (i.e. a value of K with fewer modes
  # than the maximum), use the first free one from the top for the legend;
  # otherwise, the legend is not shown
  k_free <- setdiff(k_range, modes$k[modes$m == n_col])
  if (length(k_free) > 0) {
    multi_plot <- add_inset(
      multi_plot, legend, mode_cells(k_free[1], n_col, k_range, n_col)
    )
  }

  # the labels of the columns ("Mode 1", ...), between the top of the first
  # row of plots and the top of the figure, and of the rows ("K=2", ...), in
  # the left margin (the 10% of the width reserved for them)
  top_labels <- mode_cells(k_range[1], seq_len(n_col), k_range, n_col)
  top_labels$bottom <- top_labels$top
  top_labels$top <- 1
  top_labels$text <- paste("Mode", seq_len(n_col))
  side_labels <- mode_cells(k_range, 1, k_range, n_col)
  side_labels$left <- 0
  side_labels$right <- 0.1
  side_labels$text <- paste0("K=", k_range)
  text_labels <- rbind(top_labels, side_labels)
  for (i in seq_len(nrow(text_labels))) {
    # a plot with a single text annotation (its coordinates are arbitrary, as
    # there is only one point, which is centred in the inset)
    plt <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 4, y = 25, label = text_labels$text[i]) +
      ggplot2::theme_void()
    multi_plot <- add_inset(multi_plot, plt, text_labels[i, ])
  }
  multi_plot
}

# add a plot (or a grob) as an inset filling the given cell (a list or a row
# of a data.frame with left, right, bottom and top, in npc units)
add_inset <- function(multi_plot, plt, cell) {
  multi_plot +
    patchwork::inset_element(
      plt,
      left = cell$left,
      right = cell$right,
      bottom = cell$bottom,
      top = cell$top,
      align_to = "full"
    )
}

# The cells of the grid of the 'modes' plot for modes with the given K and
# rank m (recycled), with n_col columns: a data.frame with the coordinates of
# their left, right, top and bottom sides, and of their horizontal centre, in
# npc units (i.e. in [0, 1], with y pointing up). The positions are computed
# in "grid units", measured from the top-left corner, and then rescaled:
# - each value of K takes 2 vertical units: the structure plot is 1 unit tall,
#   with half a unit of empty space above and below it, so that consecutive
#   rows are separated by 1 unit, in which the edges are drawn (the half unit
#   above the first row holds the "Mode x" labels);
# - each mode rank takes 1 horizontal unit (the insets of consecutive modes
#   touch; the white space between their bars comes from the default
#   expansion of the x scale of the structure plots), and the leftmost 10% of
#   the width is reserved for the "K=x" labels.
mode_cells <- function(k, m, k_range, n_col) {
  row_size <- 1 / (length(k_range) * 2)
  col_size <- (1 - 0.1) / n_col
  # the row is the position of K in k_range (which allows gaps between the
  # values of K); y is flipped so that the first value of K is at the top
  top <- (match(k, k_range) - 1) * 2 + 0.5
  data.frame(
    left = (m - 1) * col_size + 0.1,
    right = m * col_size + 0.1,
    centre = (m - 0.5) * col_size + 0.1,
    top = 1 - top * row_size,
    bottom = 1 - (top + 1) * row_size
  )
}

# The edges of the 'modes' plot, in npc units: for each pair of modes of
# different values of K in cost_acrossK (named "<mode1>-<mode2>", e.g.
# "K2M1-K3M2", where mode1 has the smaller K), an edge from the top centre of
# the cell of mode2 to the bottom centre of the cell of mode1, which is in the
# previous row. The data are in long format, with one row per end of each
# edge, and the two ends are joined by geom_line() as they share the same
# 'pair'.
mode_edges <- function(object, n_col) {
  pairs <- names(object$cost_acrossK)
  ends <- matrix(
    as.character(unlist(strsplit(pairs, "-", fixed = TRUE))),
    ncol = 2, byrow = TRUE
  )
  # remove links at the same level (including self links)
  keep <- mode_k(ends[, 1]) != mode_k(ends[, 2])
  ends <- ends[keep, , drop = FALSE]
  cost <- as.numeric(unlist(object$cost_acrossK))[keep]
  cell_of <- function(label) {
    mode_cells(mode_k(label), mode_rank(label), object$K_range, n_col)
  }
  tops <- cell_of(ends[, 2])
  bottoms <- cell_of(ends[, 1])
  data.frame(
    pair = rep(pairs[keep], 2),
    cost = rep(cost, 2),
    x = c(tops$centre, bottoms$centre),
    y = c(tops$top, bottoms$bottom)
  )
}

# Return the legend of a plot as a grob (or NULL if the plot has no legend),
# adapted from the unexported function of ggpubr.
# @param p a ggplot
# @param position an optional position (e.g. "bottom") to reposition the
#   legend before extracting it
.get_legend <- function(p, position = NULL) {
  if (is.null(p)) {
    return(NULL)
  }
  if (!is.null(position)) {
    p <- p + ggplot2::theme(legend.position = position)
  }
  # build the plot into a gtable, and find the grob that holds the legend
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) leg <- tmp$grobs[[leg]] else leg <- NULL
  leg
}

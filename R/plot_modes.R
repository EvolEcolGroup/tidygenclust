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
# is then placed on top of it as an inset at the same coordinates.
#
# The positions are first computed in "grid units" (see get_top_xy() and
# get_inset_corners()), measured from the top-left corner, and then rescaled
# to [0, 1] (with y pointing up):
# - each value of K takes 2 vertical units: the structure plot is 1 unit tall,
#   with half a unit of empty space above and below it, so that consecutive
#   rows are separated by 1 unit, in which the edges are drawn (the half unit
#   above the first row holds the "Mode x" labels);
# - each mode rank takes 1 horizontal unit, and the leftmost 10% of the width
#   is reserved for the "K=x" labels.
# @param object a gt_clumppling object
# @param group the (already checked) vector of group membership, or NULL
plot_modes <- function(object, group) {
  # the edges of the graph: two rows (the two ends of the edge) for each pair
  # of modes of consecutive values of K, coloured by their alignment cost
  # (with the default colour gradient, lower costs, i.e. better alignments,
  # are darker)
  connection_df <- make_connection_df(object)
  multi_plot <- ggplot2::ggplot(
    data = connection_df,
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

  # add the admixture plots as inset_elements, with the appropriate positions
  # (the 'all_modes' patchwork holds one structure plot per mode, in the same
  # order as names(object$aligned_modes))
  all_modes <- autoplot(object, type = "all_modes", group = group)
  all_modes_labels <- names(object$aligned_modes)
  inset_corners <- get_inset_corners(all_modes_labels, object$K_range)

  for (i in seq_along(all_modes_labels)) {
    # the number of replicates in this mode, shown above the right end of its
    # structure plot (outside the panel, which is drawn with clip = "off")
    n <- length(object$mode_replicates[[all_modes_labels[i]]])
    # remove the axes and labels from each structure plot (the group
    # boundaries drawn by plot_q() are kept)
    plt <- all_modes[[i]] +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "native")) +
      ggplot2::annotate(
        "text",
        x = object$N * 0.95,
        y = 1.15,
        label = paste0("(", n, ")"),
        size = 3
      )
    multi_plot <- multi_plot +
      patchwork::inset_element(
        plt,
        left = inset_corners$left[i],
        right = inset_corners$right[i],
        bottom = inset_corners$bottom[i],
        top = inset_corners$top[i],
        align_to = "full"
      )
  }

  # if a cell of the last column is free (i.e. a value of K with fewer modes
  # than the maximum), use the first free one from the top for the legend;
  # otherwise, the legend is not shown
  n_col <- max(as.integer(gsub(".*M(\\d+).*", "\\1", all_modes_labels)))
  possible_legend_slots <- paste0("K", object$K_range, "M", n_col)
  possible_legend_slots <- possible_legend_slots[
    !possible_legend_slots %in% all_modes_labels
  ]
  if (length(possible_legend_slots) > 0) {
    legend_slot <- possible_legend_slots[1]
    # the labels of the modes are passed as well, so that the corners are
    # computed with the same number of columns as for the modes
    legend_corners <- get_inset_corners(
      c(legend_slot, all_modes_labels),
      object$K_range
    )
    legend_corners <- legend_corners[1, ]
    multi_plot <- multi_plot +
      patchwork::inset_element(
        legend,
        left = legend_corners["left"],
        right = legend_corners["right"],
        bottom = legend_corners["bottom"],
        top = legend_corners["top"],
        align_to = "full"
      )
  }
  ## add top labels (Mode X)
  # we create corners for all the top row of plots
  top_labels_corners <- get_inset_corners(
    paste0("K", min(object$K_range), "M", seq_len(n_col)),
    object$K_range
  )
  # switch bottom for top, so that each label fills the space between the
  # top of its column of plots and the top of the figure
  top_labels_corners$bottom <- top_labels_corners$top
  top_labels_corners$top <- 1

  for (i in seq_len(n_col)) {
    text <- paste("Mode", i)
    # a plot with a single text annotation (its coordinates are arbitrary, as
    # there is only one point, which is centred in the inset)
    plt <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 4, y = 25, label = text) +
      ggplot2::theme_void()
    multi_plot <- multi_plot +
      patchwork::inset_element(
        plt,
        left = top_labels_corners$left[i],
        right = top_labels_corners$right[i],
        bottom = top_labels_corners$bottom[i],
        top = top_labels_corners$top[i],
        align_to = "full"
      )
  }

  ## add side labels (K=X)
  # we create corners for the first column of plots (one per value of K)
  side_labels_corners <- get_inset_corners(
    paste0("K", object$K_range, "M1"),
    object$K_range
  )
  # move the labels to the left margin (the 10% of the width reserved for
  # them), keeping the vertical position of each row
  side_labels_corners$left <- 0
  side_labels_corners$right <- 0.1

  for (i in seq_len(nrow(side_labels_corners))) {
    text <- paste0("K=", object$K_range[i])
    plt <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 4, y = 25, label = text) +
      ggplot2::theme_void()
    multi_plot <- multi_plot +
      patchwork::inset_element(
        plt,
        left = side_labels_corners$left[i],
        right = side_labels_corners$right[i],
        bottom = side_labels_corners$bottom[i],
        top = side_labels_corners$top[i],
        align_to = "full"
      )
  }
  multi_plot
}

#############################################################
# create a df of connections, with x and y in row and col units,
# as if we had all rows and all columns
# to be scaled accordingly later on
make_connection_df <- function(x) {
  # parse names of connections: each name is "<mode1>-<mode2>" (e.g.
  # "K2M1-K3M2"), where mode1 has the smaller K
  pairs_links <- strsplit(names(x$cost_acrossK), "-")
  # split them into a matrix, with the two modes of each pair in V1 and V2
  pairs_links <- matrix(
    unlist(pairs_links),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("V1", "V2"))
  )
  pairs_links <- tibble::as_tibble(pairs_links)
  pairs_links$cost <- unlist(x$cost_acrossK)
  pairs_links$pair <- names(x$cost_acrossK)
  # remove links at the same level (including self links)
  pairs_links <- pairs_links[
    gsub(".*K(\\d+)M.*", "\\1", pairs_links$V1) !=
      gsub(".*K(\\d+)M.*", "\\1", pairs_links$V2),
  ]
  # convert labels to x and y coords: each edge goes from the bottom centre
  # of the plot of the mode with the smaller K (V1) to the top centre of the
  # plot of the mode with the larger K (V2), which is in the next row. The
  # data are in long format, with one row per end of each edge, and the two
  # ends are joined by geom_line() as they share the same 'pair'
  tops <- pairs_links %>%
    dplyr::select(dplyr::all_of(c("pair", "V2", "cost"))) %>%
    dplyr::rename(label = "V2") %>%
    dplyr::mutate(position = "top")
  tops <- tops %>% dplyr::bind_cols(get_top_xy(tops$label, x$K_range))
  bottoms <- pairs_links %>%
    dplyr::select(dplyr::all_of(c("pair", "V1", "cost"))) %>%
    dplyr::rename(label = "V1") %>%
    dplyr::mutate(position = "bottom")
  bottoms <- bottoms %>%
    dplyr::bind_cols(get_bottom_xy(bottoms$label, x$K_range))
  connection_df <- rbind(tops, bottoms)
  # rescale from grid units to [0, 1], as in get_inset_corners(): a row unit
  # is 1 / (2 * number of K values), a column unit is 0.9 / number of
  # columns (after the 10% left margin), and y is flipped so that the first
  # value of K is at the top. The number of columns is the largest number of
  # modes for a K, taken from the names of all the modes (with a single K
  # there are no edges, so it cannot be taken from connection_df)
  n_row <- length(x$K_range)
  n_col <- max(as.integer(gsub(".*M(\\d+).*", "\\1", names(x$aligned_modes))))
  row_size <- 1 / (n_row * 2)
  col_size <- (1 - 0.1) / n_col
  connection_df$x <- connection_df$x * col_size + 0.1
  connection_df$y <- 1 - connection_df$y * row_size
  return(connection_df)
}

# get the coordinates for the top of plot of a given label (e.g. K2M3)
get_top_xy <- function(label, k_range) {
  # the row is the position of K in k_range (which allows gaps between the
  # values of K), and the column is the rank of the mode
  k <- as.integer(gsub(".*K(\\d+)M.*", "\\1", label))
  row <- match(k, k_range)
  col <- as.integer(gsub(".*M(\\d+).*", "\\1", label))
  # y is in row units, from the top: the plot of row r starts half a unit
  # below the top of its 2-unit band
  y <- (row - 1) * 2 + 0.5
  # x is in column units: the centre of column col
  x <- col - 0.5
  return(cbind(x, y)) # in row and col units
}

# get the coordinates for the bottom of plot of a given label (e.g. K2M3)
get_bottom_xy <- function(label, k_range) {
  xy <- get_top_xy(label, k_range)
  # each plot is one row unit tall (y increases downwards in grid units)
  xy[, 2] <- xy[, 2] + 1
  return(xy)
}

# get the corners of the inset of each label (e.g. K2M3), in npc units (i.e.
# in [0, 1], with y pointing up), as a data.frame with columns left, right,
# top and bottom and one row per label. Note that the number of columns of
# the grid is the largest mode rank among the labels, so all the labels of
# the grid should be included
get_inset_corners <- function(label, k_range) {
  # vertical position in grid units, from the top (see get_top_xy())
  k <- as.integer(gsub(".*K(\\d+)M.*", "\\1", label))
  top <- (match(k, k_range) - 1) * 2 + 0.5
  bottom <- top + 1
  # horizontal position: each mode takes a whole column unit (the insets of
  # consecutive modes touch; the white space between their bars comes from
  # the default expansion of the x scale of the structure plots)
  col <- as.integer(gsub(".*M(\\d+).*", "\\1", label))
  left <- col - 1
  right <- col
  coords <- as.data.frame(cbind(left, right, top, bottom))
  rownames(coords) <- label
  # rescale: x to [0.1, 1] (leaving the left margin for the "K=x" labels),
  # and y to [0, 1], flipping it so that it points up
  n_row <- length(k_range)
  n_col <- max(as.integer(gsub(".*M(\\d+).*", "\\1", label)))
  row_size <- 1 / (n_row * 2)
  col_size <- (1 - 0.1) / n_col
  coords$left <- coords$left * col_size + 0.1
  coords$right <- coords$right * col_size + 0.1
  coords$top <- 1 - coords$top * row_size
  coords$bottom <- 1 - coords$bottom * row_size
  return(coords)
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

plot_modes <- function(object, group){
  connection_df <- make_connection_df(object)
  multi_plot <- ggplot2::ggplot(data=connection_df,
                                ggplot2::aes(x=.data$x, y=.data$y,group=.data$pair)) +
    ggplot2::geom_line(ggplot2::aes(colour=.data$cost),linewidth =2, lineend = "round")+
    ggplot2::theme_void()+
    ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "native"))+
    ggplot2::scale_y_continuous(limits=c(0,1), expand = c(0,0)) +
    ggplot2::scale_x_continuous(limits = c(0,1), expand = c(0,0))

  legend <- .get_legend(multi_plot, position = "bottom")

  multi_plot <- multi_plot + ggplot2::guides(color = "none")

  # add the admixture plots as inset_elements, with the appropriate positions
  all_modes <- autoplot(object, type="all_modes", group = group)
  all_modes_labels <- names(object$aligned_modes)
  inset_corners <- get_inset_corners(all_modes_labels, object$K_range)

  for (i in seq_len(length(all_modes_labels))){
    k <- gsub(".*K(\\d+)M.*", "\\1", all_modes_labels[i])
    m <- as.integer(gsub(".*M(\\d+).*", "\\1", all_modes_labels[i]))
    n <- length(object$mode_replicates[[all_modes_labels[i]]])
    #n <- length(object$modes_allK[[k]][[m-1]])
    plt <- all_modes[[i]]+ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "native"))+
      ggplot2::annotate("text",x=object$N*0.95,y=1.15,
               label=paste0("(",n,")"),size=3)
    multi_plot <- multi_plot+patchwork::inset_element(plt,left = inset_corners$left[i],
                                                      right = inset_corners$right[i],
                                                      bottom = inset_corners$bottom[i],
                                                      top = inset_corners$top[i],
                                                      align_to = "full")

  }

  # if the top right panel is free, use it for the legend
  n_col <- max(as.integer(gsub(".*M(\\d+).*", "\\1", all_modes_labels)))
  possible_legend_slots <- paste0("K", object$K_range,"M",n_col)
  possible_legend_slots <- possible_legend_slots[!possible_legend_slots %in% all_modes_labels]
  if (length(possible_legend_slots)>0){
    legend_slot <- possible_legend_slots[1]
    legend_corners <- get_inset_corners(c(legend_slot,all_modes_labels),
                                         object$K_range)
    legend_corners <- legend_corners[1,]
    multi_plot <- multi_plot+patchwork::inset_element(legend,
                                                      left = legend_corners["left"],
                                                      right = legend_corners["right"],
                                                      bottom = legend_corners["bottom"],
                                                      top = legend_corners["top"],
                                                      align_to = "full")
  }
  ## add top labels (Mode X)
  # we create corners for all the top row of plots
  top_labels_corners <- get_inset_corners(
    paste0("K",min(object$K_range),"M",seq_len(n_col)), object$K_range)
  # switch bottom for top
  top_labels_corners$bottom <- top_labels_corners$top
  top_labels_corners$top <- 1

  for (i in seq_len(n_col)){
    text <- paste("Mode",i)
    plt <-   ggplot2::ggplot() +
      ggplot2::annotate("text", x = 4, y = 25, label = text) +
      ggplot2::theme_void()
    multi_plot <- multi_plot+patchwork::inset_element(plt,
                                                      left = top_labels_corners$left[i],
                                                      right = top_labels_corners$right[i],
                                                      bottom = top_labels_corners$bottom[i],
                                                      top = top_labels_corners$top[i],
                                                      align_to = "full")
  }

  ## add side labels (K=X)
  # we create corners for all the top row of plots
  side_labels_corners <- get_inset_corners(
    paste0("K",object$K_range,"M1"), object$K_range)
  # switch bottom for top
  side_labels_corners$left <- 0
  side_labels_corners$right <- 0.1

  for (i in seq_len(nrow(side_labels_corners))){
    text <- paste0("K=",object$K_range[i])
    plt <-   ggplot2::ggplot() +
      ggplot2::annotate("text", x = 4, y = 25, label = text) +
      ggplot2::theme_void()
    multi_plot <- multi_plot+patchwork::inset_element(plt,
                                                      left = side_labels_corners$left[i],
                                                      right = side_labels_corners$right[i],
                                                      bottom = side_labels_corners$bottom[i],
                                                      top = side_labels_corners$top[i],
                                                      align_to = "full")
  }
  multi_plot
}

#############################################################
# create a df of connections, with x and y in row and col units,
# as if we had all rows and all columns
# to be scaled accordingly later on
make_connection_df <- function(x){
  # parse names of connections
  pairs_links <- strsplit(names(x$cost_acrossK),"-")
  pairs_links <- matrix(unlist(pairs_links),ncol=2, byrow=TRUE,
                        dimnames=list(NULL, c("V1","V2")))
  # split them into a matrix
  pairs_links <- tibble::as_tibble(pairs_links)
  pairs_links$cost <- unlist(x$cost_acrossK)
  pairs_links$pair <- names(x$cost_acrossK)
  # remove links at the same level (including self links)
  pairs_links <- pairs_links[gsub(".*K(\\d+)M.*", "\\1", pairs_links$V1) !=
                               gsub(".*K(\\d+)M.*", "\\1", pairs_links$V2),]
  # convert labels to x and y coords
  tops <- pairs_links %>% dplyr::select(dplyr::all_of(c("pair", "V2", "cost"))) %>%
    dplyr::rename(label="V2") %>% dplyr::mutate(position="top")
  tops <- tops %>% dplyr::bind_cols(get_top_xy(tops$label, x$K_range))
  bottoms <- pairs_links %>% dplyr::select(dplyr::all_of(c("pair", "V1", "cost"))) %>%
    dplyr::rename(label="V1") %>% dplyr::mutate(position="bottom")
  bottoms <- bottoms %>% dplyr::bind_cols(get_bottom_xy(bottoms$label, x$K_range))
  connection_df <- rbind(tops,bottoms)
  # rescale
  n_row <- length(x$K_range)
  n_col <- max(as.integer(gsub(".*M(\\d+).*", "\\1", connection_df$label)))
  row_size <- 1/(n_row*2)
  col_size <- (1-0.1)/n_col
  connection_df$x <- connection_df$x * col_size + 0.1
  connection_df$y <- 1- connection_df$y * row_size
  return(connection_df)

}

# get the coordinates for the top of plot of a given label (e.g. K2M3)
get_top_xy <- function(label, k_range){
  k <- as.integer(gsub(".*K(\\d+)M.*", "\\1", label))
  row <- match(k,k_range)
  col <- as.integer(gsub(".*M(\\d+).*", "\\1", label))
  # x is in row units
  y <- (row-1) * 2+0.5
  x <- col - 0.5
  return(cbind(x,y)) # in row and col units
}

# get the coordinates for the bottom of plot of a given label (e.g. K2M3)
get_bottom_xy <- function(label, k_range){
  xy <- get_top_xy(label, k_range)
  xy[,2]<- xy[,2]+1
  return(xy)
}

get_inset_corners <- function(label, k_range){
  k <- as.integer(gsub(".*K(\\d+)M.*", "\\1", label))
  top <- (match(k,k_range)-1)*2+0.5
  bottom <- top+1
  col <- as.integer(gsub(".*M(\\d+).*", "\\1", label))
  left <- col-1
  right <- col
  coords<-as.data.frame(cbind(left, right, top, bottom))
  rownames(coords) <- label
  # rescale
  n_row <- length(k_range)
  n_col <- max(as.integer(gsub(".*M(\\d+).*", "\\1", label)))
  row_size <- 1/(n_row*2)
  col_size <- (1-0.1)/n_col
  coords$left <-coords$left * col_size+0.1
  coords$right <-coords$right * col_size+0.1
  coords$top <- 1- coords$top * row_size
  coords$bottom <-1 - coords$bottom * row_size
  return(coords)
}

# Return legend for one plot
# unexported function in ggpubr
# p is a plot
# position is an optional parameter to reposition the legend before extracting it
.get_legend <- function(p, position = NULL){

  if(is.null(p)) return(NULL)
  if(!is.null(position)){
    p <- p + ggplot2::theme(legend.position = position)
  }
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg) > 0) leg <- tmp$grobs[[leg]]
  else leg <- NULL
  leg
}


#' autoplot for clumppling objects
#'
#' An autoplot method to generate quick visualisations for [`clumppling`] objects.
#' Available types are:
#' - 'modes': all aligned modes in structure plots over a multipartite graph,
#' where better alignment between the modes is indicated by the darker color
#' of the edges connecting their structure plots, and the cost of optimal
#' alignment is labelled on each edge.
#' - 'modes_within_K': A set of figures, one for each number of clusters, with
#' all modes with the same number of clusters in structure plots in one figure.
#' - 'major_modes': the major modes of each K aligned in a series of structure
#' plots.
#' - 'all_modes': all aligned modes in a series of structure plots.
#'
#' Currently `autoplot` wraps the python plotting functions, which use
#' `matlibplot`. In the future, these will be replaced with `ggplot2` native
#' plots.
#'
#' @param object a [`clumppling`] object
#' @param type the type of plot, one of 'modes', 'modes_within_K',
#' 'major_modes' or 'all_modes'.
#' @param group a vector of membership to a-priori groups (e.g. populations)
#' @param k the k value to be plotted if 'type' is 'modes_within_k'
#' @param ... not used at the moment
#' @returns a plot
#' @rdname autoplot_gt_clumppling
#' @export

autoplot.gt_clumppling <- function(object,
                                type = c("modes", "modes_within_k",
                                         "major_modes", "all_modes"),
                                group = NULL,
                                k = NULL,
                                ...) {
  type <- match.arg(type)
  rlang::check_dots_empty()
  # check that group, if given, is cohere with the rest of the data
  if (!is.null(group)){
    # first check that labels are in blocks
    if (length(rle(group)$values)!=length(unique(group))) {
      stop("values in 'group' are not ordered (they should be in consecutive blocks, one per group")
    }
    if (length(group)!=object$N){
      stop("'groups' should be of the same lenght as the original data (as found in object$N)")
    }
    group_x <- cumsum(table(forcats::fct_inorder(group)))
  } else {
    group_x <- NULL
  }

 if (type=="modes"){
   plot_modes(object, group)
 } else if (type == "modes_within_k") {
   plot_modes_within_k(object, group_x, k=k)
 } else if (type == "major_modes") {
   plot_major_modes(object, group_x)
 } else if (type == "all_modes") {
   plot_all_modes(object, group_x)
 }
}

# Plot major modes
plot_major_modes <- function(object, group_x){
  major_modes <- tidy(object, matrix="major_modes")
  y_labels <- rep("",nrow(major_modes))
  y_labels[major_modes$m==1]<-paste0("K = ",major_modes$k[major_modes$m==1])
  plot_list<-lapply(1:nrow(major_modes),plot_q_from_list,
                    object$aligned_modes[major_modes$label],y_labels, group_x)
  patchwork::wrap_plots(plot_list, ncol = 1) +
    patchwork::plot_layout(axes = "collect")
}

# Plot all modes
plot_all_modes <- function(object, group_x){
  all_modes <- tidy(object, matrix="modes")
  y_labels <- rep("",nrow(all_modes))
  y_labels[all_modes$m==1]<-paste0("K = ",all_modes$k[all_modes$m==1])
  plot_list<-lapply(1:nrow(all_modes),plot_q_from_list,object$aligned_modes,y_labels, group_x)
  patchwork::wrap_plots(plot_list, ncol = 1) +
    patchwork::plot_layout(axes = "collect")
}

# Plot major modes
plot_modes_within_k <- function(object, group_x, k){
  if(is.null(k)){
    stop("k must be provided for modes_within_k")
  }
  #browser()
  all_modes <- tidy(object, matrix="modes")
  k_modes <- all_modes [all_modes$k == k,]
  y_labels <- k_modes$label
  plot_list<-lapply(seq_len(nrow(k_modes)),plot_q_from_list,
                    object$aligned_modes[k_modes$label],y_labels, group_x)
  patchwork::wrap_plots(plot_list, ncol = 1) +
    patchwork::plot_layout(axes = "collect")
}

# a compact q plot to use as a panel
# @param q_tidied a q values tidied with tidy
# @ the y label if present
# @ group_x the x values where groups switch, a black line will be plotted here
plot_q <- function(q_tidied, y_lab = "", group_x){
  plt <- ggplot2::ggplot(q_tidied,
                ggplot2::aes(.data$id,
                    .data$percentage,
                    fill = .data$q)) +
    ggplot2::geom_col(width = 1,
                      position = ggplot2::position_stack(reverse = TRUE))+
    # ggplot2::theme_minimal() + # remove most thick marks etc
    #  ggplot2::theme( # adjust title position and remove panel grid
    #    panel.grid = ggplot2::element_blank(),
    #    axis.text.y = ggplot2::element_blank(),
    #    axis.title.x = ggplot2::element_blank()
    #  ) +
    tidypopgen::theme_distruct()+
     ggplot2::labs(y = y_lab)+
     ggplot2::coord_cartesian(ylim=c(0,1),clip="off")+
     tidypopgen::scale_fill_distruct()
  # add vertical lines if we have groups
  if (!is.null(group_x)){
    segment_data = data.frame(
      x = group_x+0.5,
      xend = group_x+0.5,
      y = rep(0,length(group_x)),
      yend = rep(1, length(group_x))
    )

    plt <- plt + ggplot2::geom_segment(data = segment_data,
                             ggplot2::aes(x = .data$x,
                                          y = .data$y,
                                          xend = .data$xend,
                                          yend = .data$yend),
                             inherit.aes = FALSE)
    get_mid_points <- function (vec){
      (vec[-length(vec)] + vec[-1L])/2.
    }
    plt <- plt + ggplot2::scale_x_continuous(breaks = get_mid_points(c(0,group_x)),
                                         labels = names(group_x))+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1))
  } else {
    plt <- plt +ggplot2::theme(axis.text.x = ggplot2::element_blank(),
)
  }

  plt
}

# plot the i element from a q list
plot_q_from_list <- function(i, q_list, y_labs, group_x){
  plot_q(tidy_q(q_list[[i]]), y_labs[i], group_x)
}

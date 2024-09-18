#' run clumppling
#'
#' This function runs the clumppling algorithm.
#'
#' @param input_path the path where the Q files are stored, either a directory
#' or a zip archive
#' @param input_format a string defining the format of the input files, one of
#' 'admixture' (default)
#' @param cd_param the cd_param
#' use_rep boolean on whether a representative repeat should be used as a consesus
#' for a mode. Defaults to FALSE, which leads to the computatoin of an average
#' @param merge_cls boolean,
#' @param cd_default boolean
#' @param use_rep boolean
#' @param custom_cmap a custom colour map (this should be moved
#' to the plotting functions!) a string that is either empty (default “”)
#' or contains a list of colors to be used for a customized colormap. If
#' the string is empty then the default colormap will be used.
#' The customized colormap, if provided, should be a list of colors in hex
#' code in a comma-delimited string. An example colormap with five
#' colors looks like "#FF0000,#FFFF00,#00EAFF #AA00FF,#FF7F00". If
#' the provided colormap does not have enough colors for all clusters,
#' colors will be cycled.
#' @param output_path (optional) the clumppling functions in python save
#' everything to file. By default, R stores the information in objects in the
#' environment, and sends those files to a temporary directory
#' that will be cleared at the end of a session. `output_path` allows to
#' change the location
#' of those files. This is only useful to those interested in recovering the
#' same files as created by the python clumppling module, or for debugging.
#' @returns a list of class `clumpling` TODO describe the elements of the list
#' @export

clumppling <- function (input_path,
                        input_format = "admixture",
                        cd_param =1.0,
                        use_rep = 0,
                        merge_cls=0,
                        cd_default=1,
                        custom_cmap = "",
                        output_path = tempfile("clumpling")){
 if (tools::file_ext(input_path)=="zip"){
   temp_zip_path <- tempfile()
   utils::unzip(input_path,exdir = temp_zip_path)
   input_path <- temp_zip_path
 }
  # create a namespace object with all the inputs
  argparse <- reticulate::import("argparse")
  rclump_args <- argparse$Namespace(input_path = input_path,
                                    output_path = output_path,
                                    input_format = input_format,
                                    vis = 1,
                                    cd_param =cd_param,
                                    use_rep = use_rep,
                                    merge_cls= merge_cls,
                                    cd_default=cd_default,
                                    plot_modes = 0,
                                    plot_modes_withinK = 0,
                                    plot_major_modes = 0,
                                    plot_all_modes=0,
                                    custom_cmap = custom_cmap)


  clump_res<-.py_rclumppling$clumppling_run(args = rclump_args)
  names(clump_res) <- c("args", "cmap", "Q_list", "K_list", "Q_files", "R",
                        "N", "K_range", "K_max", "K2IDs", "alignment_withinK",
                        "cost_withinK", "modes_allK", "cost_matrices", "msg",
                        "mode_labels", "rep_modes", "repQ_modes", "avgQ_modes",
                        "alignment_to_modes", "stats", "costs",
                        "alignment_acrossK_cons", "cost_acrossK_cons",
                        "best_acrossK_cons")
  # cast back to integers
  # (reticulate casts int64 to double, which is incompatible with certain
  # functions)
  clump_res$K_range <- as.integer(clump_res$K_range)
  clump_res$alignment_acrossK_cons <- lapply(clump_res$alignment_acrossK_cons,
                                             as.integer)
  # read modes_aligned (used for ggplots)
  modes_aligned_path <- file.path(clump_res$args$output_path,
                                             "modes_aligned")
  q_files <- list.files(modes_aligned_path, pattern="Q")
  clump_res$aligned_modes<-lapply(file.path(modes_aligned_path,q_files),
                                  utils::read.table)
  get_code <- function(x){strsplit(x,"_")[[1]][1]}
  names(clump_res$aligned_modes) <- unlist(lapply(q_files,get_code))

  class(clump_res) <- "clumppling"
  return(clump_res)
}

#' run clumppling
#'
#' This function runs the clumppling algorithm.
#'
#' If you would like to generate an annotated autoplot from your gt_clumppling object,
#' ensure that all individuals from the same population are adjacent to one another
#' in the Q-matrix or gt_admix object supplied to gt_clumppling. Autoplot 'group'
#' argument requires that all individuals from the same group are adjacent.
#'
#' @param input_path the path where the Q files are stored, either a directory
#' or a zip archive, or a `q_matrix_list` object
#' @param input_format a string defining the format of the input files, one of
#' 'admixture' (default)
#' @param cd_param the cd_param
#' use_rep boolean on whether a representative repeat should be used as a consesus
#' for a mode. Defaults to FALSE, which leads to the computatoin of an average
#' @param merge_cls boolean,
#' @param cd_default boolean
#' @param use_rep boolean
#' @param output_path (optional) the clumppling functions in python save
#' everything to file. By default, R stores the information in objects in the
#' environment, and sends those files to a temporary directory
#' that will be cleared at the end of a session. `output_path` allows to
#' change the location
#' of those files. This is only useful to those interested in recovering the
#' same files as created by the python clumppling module, or for debugging.
#' @returns a list of class `clumpling` TODO describe the elements of the list
#' @export

gt_clumppling <- function (input_path,
                        input_format = "admixture",
                        cd_param = 1.0,
                        use_rep = 0,
                        merge_cls= 0,
                        cd_default= 1,
                        output_path = tempfile("clump_out")){

  # Check if input_path is a zip file and unzip it
  if (is.character(input_path)) {
    if (tools::file_ext(input_path) == "zip") {
      temp_zip_path <- tempfile()
      utils::unzip(input_path, exdir = temp_zip_path)
      input_path <- temp_zip_path
    }
  }

  # check if input_path is a "gt_admix" object
  if (inherits(input_path, "gt_admix")) {
    # extract flat q-matrix list
    q_mat_list <- input_path$Q
    #check q_mat_list has more than 1 q_matrix
    if (length(q_mat_list) < 2) {
      stop("Input list must contain more than one Q-matrix")
    }

    # create different temp dir each run
    temp_q_dir <- tempfile(pattern = "clumppling_q_dir")
    dir.create(temp_q_dir)
    invisible(lapply(seq_along(q_mat_list), function(i) {
      # create a unique filename for each list item
      file_name <- paste0("gt_admix_", i, ".Q")
      q_filename <- file.path(temp_q_dir, file_name)
      utils::write.table(q_mat_list[[i]], file = q_filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
    }))
    input_path <- temp_q_dir
  }

  # create command line for clumppling
  clump_args <- paste0("-m clumppling -i ", input_path, " -o ", output_path,
                  " -f ", input_format, " -v=1 --cd_param=", cd_param,
                  " --use_rep=", use_rep, " --merge_cls=", merge_cls, " --cd_default=", cd_default,
                  " --plot_modes=0 --plot_modes_withinK=0 --plot_major_modes=0 --plot_all_modes=0")
  reticulate::conda_run2(cmd = "python",
                         args = clump_args,
                         envname = "cclumppling")
  # now we read the output files and create an output object



  # elements that we use:
  # modes_allK (only used to find n when plotting modes, can we do it in another way???)
  # cost_acrossK_cons

  # clump_res<-.py_rclumppling$clumppling_run(args = rclump_args)
  # names(clump_res) <- c("args", "cmap", "Q_list", "K_list", "Q_files", "R",
  #                       "N", "K_range", "K_max", "K2IDs", "alignment_withinK",
  #                       "cost_withinK", "modes_allK", "cost_matrices", "msg",
  #                       "mode_labels", "rep_modes", "repQ_modes", "avgQ_modes",
  #                       "alignment_to_modes", "stats", "costs",
  #                       "alignment_acrossK_cons", "cost_acrossK_cons",
  #                       "best_acrossK_cons")
  # # cast back to integers
  # # (reticulate casts int64 to double, which is incompatible with certain
  # # functions)
  # clump_res$K_range <- as.integer(clump_res$K_range)
  # clump_res$alignment_acrossK_cons <- lapply(clump_res$alignment_acrossK_cons,
  #                                            as.integer)

  clump_res <- list()
  # read modes_aligned (used for ggplots)
  modes_aligned_path <- file.path(output_path,
                                             "modes_aligned")
  q_files <- list.files(modes_aligned_path, pattern="Q")
  clump_res$aligned_modes<-lapply(file.path(modes_aligned_path,q_files),
                                  utils::read.table)
  get_code <- function(x){strsplit(x,"_")[[1]][1]}
  names(clump_res$aligned_modes) <- unlist(lapply(q_files,get_code))
  # read cost_acrossK
  cost_acrossK_path <- file.path(output_path,
                                 "alignment_acrossK")
  cost_file <- list.files(cost_acrossK_path, pattern="alignment_acrossK")
  costs <- utils::read.csv(file.path(cost_acrossK_path,cost_file[1]))
  clump_res$cost_acrossK <- as.list(costs$Cost)
  names(clump_res$cost_acrossK) <- costs$Mode1.Mode2
  # read modes
  modes_path <- file.path(output_path,
                          "modes")
  mode_alignments <- utils::read.csv(file.path(modes_path,"mode_alignments.txt"))

  # for each mode, get the replicates
  # add replicate id column by parsing string and getting number after R
  mode_alignments$ReplicateID <- as.integer(sub(".*R(\\d+).*", "\\1", mode_alignments$Replicate))
  split(mode_alignments$ReplicateID, mode_alignments$Mode) -> mode_replicates
  clump_res$mode_replicates <- mode_replicates

  # add K_range
  mode_alignments$K <- as.integer(sub(".*K(\\d+).*", "\\1", mode_alignments$Mode))
  clump_res$K_range <- unique(mode_alignments$K)

  # add N
  clump_res$N <- nrow(clump_res$aligned_modes[[1]])

  class(clump_res) <- "gt_clumppling"
  return(clump_res)
}

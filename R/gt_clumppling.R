#' run clumppling
#'
#' This function runs the clumppling algorithm.
#'
#' If you would like to generate an annotated autoplot from your gt_clumppling
#' object, ensure that all individuals from the same population are adjacent to
#' one another in the Q-matrix or gt_admix object supplied to gt_clumppling.
#' Autoplot 'group' argument requires that all individuals from the same group
#' are adjacent.
#'
#' @param input_path the path where the Q files are stored, either a directory
#'   or a zip archive, or a `q_matrix_list` object
#' @param input_format a string defining the format of the input files, one of
#'   'admixture' (default),'structure','fastStructure' or 'generalQ'
#' @param cd_method the community detection method to use, one of 'louvain'
#'  (default), 'leiden', 'infomap', 'markov_clustering', 'label_propagation',
#'   'walktrap', 'custom'
#' @param use_rep boolean, whether to use representative modes (alternative:
#'  average), defaults to TRUE
#' @param merge boolean, whether to merge two clusters when aligning K+1 to K,
#'  defaults to TRUE
#' @param use_best_pair boolean, whether to use best pair as anchor for across-K
#'  alignment (alternative: major), defaults to TRUE
#' @param extension (optional) if loading from files rather than a
#'   `q_matrix_list` object specify the extension e.g. ".Q" or ".indivq"
#' @param output_path (optional) the clumppling functions in python save
#'   everything to file. By default, R stores the information in objects in the
#'   environment, and sends those files to a temporary directory that will be
#'   cleared at the end of a session. `output_path` allows to change the
#'   location of those files. This is only useful to those interested in
#'   recovering the same files as created by the python clumppling module, or
#'   for debugging.
#' @returns a list of class `gt_clumppling` containing:
#' - N: number of individuals
#' - K_range: vector of K values analyzed
#' - mode_replicates: a list of replicate indices for each mode
#' - cost_acrossK: a named list of costs for each pairwise K alignment
#' - aligned_modes: a list of data.frames, each data.frame is a Q-matrix
#' @export

gt_clumppling <- function(
  input_path,
  output_path = tempfile("clump_out"),
  input_format = "admixture",
  use_rep = TRUE,
  merge = TRUE,
  cd_method = "louvain",
  use_best_pair = TRUE,
  extension = ".Q"
) {
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
    # check q_mat_list has more than 1 q_matrix
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
      utils::write.table(
        q_mat_list[[i]],
        file = q_filename,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
      )
    }))
    input_path <- temp_q_dir
  }

  if (merge == TRUE) {
    merge <- "True"
  } else if (merge == FALSE) {
    merge <- "False"
  }

  if (use_rep == TRUE) {
    use_rep <- "True"
  } else if (use_rep == FALSE) {
    use_rep <- "False"
  }

  if (use_best_pair == TRUE) {
    use_best_pair <- "True"
  } else if (use_best_pair == FALSE) {
    use_best_pair <- "False"
  }


  clump_args <- paste0(
    "-m clumppling ",
    "-i ", input_path,
    " -o ", output_path,
    " -f ", input_format,
    " -v=False",
    " --cd_method=", cd_method,
    " --test_comm=False",
    " --comm_min=1e-6", # set to default
    " --comm_max=1e-2", # set to default
    " --merge=", merge,
    " --use_rep=", use_rep,
    " --use_best_pair=", use_best_pair,
    " --extension=", extension
  )


  reticulate::conda_run2(
    cmd = "python",
    args = clump_args,
    envname = "cclumppling"
  )
  # now we read the output files and create an output object
  clump_res <- list()
  # read modes_aligned (used for ggplots)
  modes_aligned_path <- file.path(
    output_path,
    "modes_aligned"
  )
  q_files <- list.files(modes_aligned_path, pattern = "Q")
  clump_res$aligned_modes <- lapply(
    file.path(modes_aligned_path, q_files),
    utils::read.table
  )
  get_code <- function(x) {
    strsplit(x, "_")[[1]][1]
  }
  names(clump_res$aligned_modes) <- unlist(lapply(q_files, get_code))

  # Reorder the list, if needed
  k <- as.numeric(gsub(".*K(\\d+)M.*", "\\1", names(clump_res$aligned_modes)))
  indices <- order(k)
  clump_res$aligned_modes <- clump_res$aligned_modes[indices]

  # read cost_acrossK
  cost_across_k_path <- file.path(
    output_path,
    "alignment_acrossK"
  )
  cost_file <- list.files(cost_across_k_path, pattern = "alignment_acrossK")
  costs <- utils::read.csv(file.path(cost_across_k_path, cost_file[1]))
  clump_res$cost_acrossK <- as.list(costs$Cost)
  names(clump_res$cost_acrossK) <- costs$Mode1.Mode2
  # read modes
  modes_path <- file.path(
    output_path,
    "modes"
  )
  mode_alignment <- utils::read.csv(file.path(
    modes_path,
    "mode_alignment.txt"
  ))

  # for each mode, get the replicates
  # add replicate id column by parsing string and getting number after R
  mode_alignment$ReplicateID <- as.integer(sub(
    ".*R(\\d+).*",
    "\\1",
    mode_alignment$Replicate
  ))
  mode_replicates <- split(mode_alignment$ReplicateID, mode_alignment$Mode)
  clump_res$mode_replicates <- mode_replicates

  # Reorder the list, if needed
  clump_res$mode_replicates <- clump_res$mode_replicates[indices]

  # add K_range
  mode_alignment$K <- as.integer(sub(
    ".*K(\\d+).*",
    "\\1",
    mode_alignment$Mode
  ))
  clump_res$K_range <- unique(mode_alignment$K)

  # add N
  clump_res$N <- nrow(clump_res$aligned_modes[[1]])

  class(clump_res) <- "gt_clumppling"
  return(clump_res)
}


#' Subset a `gt_clumppling` object
#'
#' This function subsets `gt_clumppling` objects to a set of individuals or a
#' set of values of K. This is intended to create plot insets, or to visualise a
#' subset of individuals during data analysis. To understand the modes within a
#' subset of individuals in your data, you should subset your `gt_admix` object
#' and re-run `gt_clumppling`.
#'
#' @param x a gt_clumppling object
#' @param indivs a vector of individual indices to keep
#' @param k a vector of k values to subset to
#' @returns a gt_clumppling object subsetted to the individuals specified
#' @export
subset_gt_clumppling <- function(x, k = NULL, indivs = NULL) {
  if (!is.null(indivs)) {
    indivs <- as.integer(indivs)
    # check indivs is an integer vector of consecutive values from min to max
    indivs_sorted <- sort(indivs)

    if (!identical(indivs, indivs_sorted) ||
      !all(diff(indivs) == 1)) { # nolint
      stop( # nolint
        "indivs values must be strictly increasing and consecutive from min ",
        "to max. If your individuals are not consecutive, please reorder your ",
        "individuals using gt_admix_reorder_q() and re-run gt_clumppling."
      )
    }

    x$N <- length(indivs)
    x$aligned_modes <- lapply(
      x$aligned_modes,
      function(mat) {
        mat[indivs, , drop = FALSE]
      }
    )

    # add an attribute to indicate this gt_clumppling has been subset by indiv
    attr(x, "subset_indivs") <- TRUE
  }

  if (!is.null(k)) {
    k <- as.integer(k)
    k_sorted <- sort(k)
    expected <- seq(min(k_sorted), max(k_sorted))

    # check k is an integer vector of consecutive values from min to max
    if (!identical(k, expected)) {
      stop("k values must be consecutive from min to max")
    }

    # find any x$K_range not in k
    k_remove <- which(!x$K_range %in% k)
    k_rm <- x$K_range[k_remove]

    # k values must be consecutive from min to max
    k_indices <- which(x$K_range %in% k)
    x$K_range <- x$K_range[k_indices]
    # find the entries in $aligned_modes that = paste0("K",K_range[])
    matched_indices <- c()
    for (i in k) {
      this_k <- paste0("K", i)
      this_k_indices <- grep(this_k, names(x$aligned_modes))
      matched_indices <- c(matched_indices, this_k_indices)
    }
    # remove matched_indices from x$aligned_modes
    x$aligned_modes <- x$aligned_modes[matched_indices]
    x$mode_replicates <- x$mode_replicates[matched_indices]

    # find the entries in $cost_acrossK that need removing
    matched_indices <- c()
    for (i in k_rm) {
      this_k <- paste0("K", i)
      rm_indices <- grep(this_k, names(x$cost_acrossK))
      matched_indices <- c(matched_indices, rm_indices)
    }
    x$cost_acrossK <- x$cost_acrossK[-matched_indices]
  }

  return(x)
}

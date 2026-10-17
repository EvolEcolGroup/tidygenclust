#' run clumppling
#'
#' This function runs the clumppling algorithm by calling the python module. If
#' you use this function, make sure that you cite the relevant paper by Liu,
#' Kopelman, and Rosenberg (2024).
#'
#' If you would like to generate an annotated autoplot from your gt_clumppling
#' object, ensure that all individuals from the same population are adjacent to
#' one another in the Q-matrix or gt_admix object supplied to gt_clumppling.
#' Autoplot 'group' argument requires that all individuals from the same group
#' are adjacent.
#'
#' Q matrices from a `gt_admix` object (or a list of Q matrices) are passed
#' directly to python. Q files in a directory (or a zip archive) are parsed by
#' clumppling itself (using a temporary directory), which allows the use of the
#' different formats supported by clumppling. Note that, for 'admixture' files
#' with extension '.indivq', and for 'structure' files, clumppling uses the
#' population labels in the files to reorder individuals so that individuals
#' from the same population are adjacent.
#'
#' @references X. Liu, N. M. Kopelman, N. A. Rosenberg (2024) Clumppling:
#'   cluster matching and permutation program with integer linear programming.
#'   Bioinformatics 40(1): btad751. doi:
#'   https://doi.org/10.1093/bioinformatics/btad751
#'
#' @param input_path either a `gt_admix` object, a list of Q matrices, or the
#'   path where the Q files are stored, either a directory or a zip archive
#' @param input_format a string defining the format of the input files, one of
#'   'admixture' (default),'structure','fastStructure' or 'generalQ'. Only used
#'   if `input_path` is a path.
#' @param cd_method the community detection method to use, one of 'louvain'
#'  (default), 'leiden', 'infomap', 'markov_clustering', 'label_propagation',
#'   'walktrap', 'custom'
#' @param use_rep boolean, whether to use representative modes (alternative:
#'  average), defaults to TRUE
#' @param merge boolean, whether to merge two clusters when aligning K+1 to K,
#'  defaults to TRUE
#' @param use_best_pair boolean, whether to use best pair as anchor for across-K
#'  alignment (alternative: major), defaults to TRUE
#' @param extension (optional) if loading from files rather than a `gt_admix`
#'   object, specify the extension e.g. ".Q" or ".indivq"
#' @param seed (optional) an integer used to set the random seed in python, so
#'   that community detection (e.g. 'louvain') is reproducible. Defaults to
#'   NULL (no seed is set).
#' @param output_path deprecated and ignored. Clumppling functions are now
#'   called directly from R, and their outputs are returned as R objects rather
#'   than written to files.
#' @returns a list of class `gt_clumppling` containing:
#' - N: number of individuals
#' - K_range: vector of K values analyzed
#' - mode_replicates: a list of replicate indices for each mode
#' - cost_acrossK: a named list of costs for each pairwise K alignment
#' - aligned_modes: a list of data.frames, each data.frame is a Q-matrix
#' - mode_stats: a data.frame with the size, cost and performance (average
#'   similarity among replicates) of each mode
#' @export

gt_clumppling <- function(
  input_path,
  input_format = "admixture",
  use_rep = TRUE,
  merge = TRUE,
  cd_method = "louvain",
  use_best_pair = TRUE,
  extension = ".Q",
  seed = NULL,
  output_path = NULL
) {
  .ensure_clumppling_python()

  if (!is.null(output_path)) {
    warning(
      "'output_path' is deprecated and ignored: clumppling results are ",
      "returned as R objects and are no longer written to files"
    )
  }

  clump_args <- list(
    cd_method = cd_method,
    cd_res = 1.0, # default
    test_comm = FALSE,
    comm_min = 1e-6, # default
    comm_max = 1e-2, # default
    merge = merge,
    use_rep = use_rep,
    use_best_pair = use_best_pair,
    random_seed = if (is.null(seed)) NULL else as.integer(seed)
  )

  if (inherits(input_path, "gt_admix")) {
    q_mat_list <- input_path$Q
  } else if (is.list(input_path) && !is.data.frame(input_path)) {
    q_mat_list <- input_path
  } else if (is.character(input_path) && length(input_path) == 1) {
    q_mat_list <- NULL
  } else {
    stop(
      "input_path must be a gt_admix object, a list of Q matrices, or the ",
      "path to a directory or zip archive of Q files"
    )
  }

  if (!is.null(q_mat_list)) {
    if (length(q_mat_list) < 2) {
      stop("Input list must contain more than one Q-matrix")
    }
    # strip any class (e.g. q_matrix) so that reticulate sends plain numpy
    # arrays to python
    clump_args$Q_list <- lapply(q_mat_list, function(q) {
      q <- as.matrix(q)
      matrix(as.numeric(q), nrow = nrow(q), ncol = ncol(q))
    })
  } else {
    # Check if input_path is a zip file and unzip it
    if (tools::file_ext(input_path) == "zip") {
      temp_zip_path <- tempfile()
      utils::unzip(input_path, exdir = temp_zip_path)
      on.exit(unlink(temp_zip_path, recursive = TRUE), add = TRUE)
      input_path <- temp_zip_path
    }
    # clumppling parses the input files and writes the processed Q matrices
    # to a directory, which we then read back in python
    processed_dir <- tempfile("clumppling_input")
    on.exit(unlink(processed_dir, recursive = TRUE), add = TRUE)
    clump_args <- c(
      clump_args,
      list(
        input_dir = normalizePath(input_path, mustWork = TRUE),
        processed_dir = processed_dir,
        fmt = input_format,
        extension = extension
      )
    )
  }

  clump_py <- do.call(.py_rclumppling$clumppling_run, clump_args)

  # now create the output object
  clump_res <- list()
  clump_res$aligned_modes <- lapply(clump_py$aligned_modes, as.data.frame)
  names(clump_res$aligned_modes) <- clump_py$mode_names
  clump_res$cost_acrossK <- clump_py$cost_acrossK

  # for each mode, get the replicates (the replicate id is the number after R
  # in the name of the replicate, e.g. 3_K2R1)
  mode_alignment <- clump_py$mode_alignment
  replicate_id <- as.integer(sub(".*R(\\d+)$", "\\1", mode_alignment$Replicate))
  clump_res$mode_replicates <- split(
    replicate_id,
    factor(mode_alignment$Mode, levels = clump_py$mode_names)
  )

  clump_res$K_range <- as.integer(clump_py$K_range)
  clump_res$N <- as.integer(clump_py$N)
  clump_res$mode_stats <- clump_py$mode_stats
  # drop the reference to the python index, which would not survive saving
  attr(clump_res$mode_stats, "pandas.index") <- NULL

  class(clump_res) <- "gt_clumppling"
  return(clump_res)
}

.ensure_clumppling_python <- function() {
  reticulate::use_condaenv("ctidygenclust", required = TRUE)
  # get info on the current python environment
  cfg <- reticulate::py_config()

  if (!reticulate::py_module_available("clumppling")) {
    stop(
      "Python environment selected but 'clumppling' is not importable.\n",
      "Current python: ", cfg$python
    )
  }
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
    if (!is.null(x$mode_stats)) {
      x$mode_stats <- x$mode_stats[
        x$mode_stats$Mode %in% names(x$aligned_modes), ,
        drop = FALSE
      ]
    }

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

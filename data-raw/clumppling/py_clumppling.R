# Run the python implementation of clumppling (2.0.10) from R, returning the
# same gt_clumppling objects as the python-based gt_clumppling() of previous
# versions of tidygenclust (this is that function, which called
# data-raw/clumppling/py_rclumppling.py through reticulate).
#
# Sourced by regenerate_python_fixtures.R and validate_vs_python.R, from the
# root of the package, in the conda environment created by
# setup_clumppling_env.sh (set CLUMPPLING_REFERENCE_ENV to use a differently
# named environment).

reticulate::use_condaenv(
  Sys.getenv("CLUMPPLING_REFERENCE_ENV", "clumppling_reference"),
  required = TRUE
)
py_rclumppling <- reticulate::import_from_path("py_rclumppling",
  path = file.path("data-raw", "clumppling")
)
# clumppling processes Q files in the arbitrary order returned by
# os.listdir(); we sort them, as the example gt_admix objects of tidygenclust
# do (see data-raw/clumppling/create_example_gt_admix.R)
reticulate::py_run_string("
import os, types
import clumppling.parseInput.utils as _pu
_pu.os = types.SimpleNamespace(path=os.path, makedirs=os.makedirs,
                               listdir=lambda d: sorted(os.listdir(d)))
")

# x is a gt_admix object, a list of Q matrices, or the path of a directory or
# zip archive of Q files (read by clumppling with `input_format` and
# `extension`); seed is the seed of numpy's random number generator
py_clumppling <- function(x, input_format = "admixture", use_rep = TRUE,
                          merge = TRUE, cd_method = "louvain",
                          use_best_pair = TRUE, extension = ".Q",
                          seed = NULL, test_comm = FALSE) {
  clump_args <- list(
    cd_method = cd_method,
    cd_res = 1.0,
    test_comm = test_comm,
    comm_min = 1e-6,
    comm_max = 1e-2,
    merge = merge,
    use_rep = use_rep,
    use_best_pair = use_best_pair,
    random_seed = if (is.null(seed)) NULL else as.integer(seed)
  )
  if (inherits(x, "gt_admix")) {
    x <- x$Q
  }
  if (is.list(x)) {
    # strip any class (e.g. q_matrix) so that reticulate sends plain numpy
    # arrays to python
    clump_args$Q_list <- lapply(x, function(q) {
      q <- as.matrix(q)
      matrix(as.numeric(q), nrow = nrow(q), ncol = ncol(q))
    })
  } else {
    if (tools::file_ext(x) == "zip") {
      temp_zip_path <- tempfile()
      utils::unzip(x, exdir = temp_zip_path)
      on.exit(unlink(temp_zip_path, recursive = TRUE), add = TRUE)
      x <- temp_zip_path
    }
    # clumppling parses the input files and writes the processed Q matrices
    # to a directory, which are then read back in python
    processed_dir <- tempfile("clumppling_input")
    on.exit(unlink(processed_dir, recursive = TRUE), add = TRUE)
    clump_args <- c(clump_args, list(
      input_dir = normalizePath(x, mustWork = TRUE),
      processed_dir = processed_dir,
      fmt = input_format,
      extension = extension
    ))
  }

  # python prints a message for each matrix with K = 1
  reticulate::py_capture_output(
    clump_py <- do.call(py_rclumppling$clumppling_run, clump_args)
  )

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
  clump_res
}

# the versions of python and of the main python packages
py_clumppling_versions <- function() {
  metadata <- reticulate::import("importlib.metadata")
  pkgs <- c(
    "clumppling", "cdlib", "networkx", "numpy", "python-louvain",
    "python-igraph", "cvxpy", "scipy"
  )
  c(
    python = as.character(reticulate::py_config()$version),
    vapply(pkgs, function(pkg) metadata$version(pkg), character(1))
  )
}

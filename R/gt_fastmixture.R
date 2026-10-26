#' fastmixture algorithm for population genetics clustering
#'
#' This function implements the fastmixture algorithm for population genetics
#' clustering by calling the python module. If you use this function, make sure
#' that you cite the relevant paper by Santander, Refoyo-Martínez, and Meisner
#' (2024).
#'
#' This function returns a `gt_admix` object that can be plotted with
#' `autoplot`, and tidied with `tidy` methods from the `tidypopgen` package.
#' Cross-validation is not performed by default (`cv = NULL`); if you want to
#' include cross-validation you can set `cv` to the number of folds, a value
#' greater than 1 (ADMIXTURE performs 5-fold cv as default).
#'
#' `fastmixture` is an optional tool, which is not needed by the rest of
#' `tidygenclust` (see [gt_openadmixture()] for a native implementation of
#' the ADMIXTURE model). Install it with
#' `tgc_tools_install(tools = "fastmixture")`, or add it to an existing
#' `ctidygenclust` conda environment with [tgc_fastmixture_install()].
#'
#' @references C. G. Santander, A. Refoyo Martinez, J. Meisner (2024) Faster
#'   model-based estimation of ancestry proportions. bioRxiv 2024.07.08.602454;
#'   doi: https://doi.org/10.1101/2024.07.08.602454
#'
#' @param x either a [`tidypopgen::gen_tibble`], or the name of the binary plink
#'   file (without the .bed extension)
#' @param k the number of ancestral components (clusters), either a single value
#'   or a vector
#' @param n_runs the number of repeats for each k value
#' @param threads the number of threads to use (1)
#' @param seed the random seed (defaults to 42): a vector of length
#'   `n_runs`, with one seed per run (the same seeds are used for each value
#'   of `k`), or a single seed, used for the first run and incremented for
#'   the following ones
#' @param iter the maximum number of iterations (1000)
#' @param tole the tolerance in log-likelihood units between iterations (1e-9)
#' @param batches the number of maximum mini-batches (32)
#' @param supervised the name of the file with the supervised labels (NULL)
#' @param check the number of iterations to check for convergence (5)
#' @param power number of power iterations in randomised SVD (11)
#' @param chunk the number of SNPs in chunk operations (8192)
#' @param subsample Fraction of SNPs to subsample in SVD/ALS (0.7)
#' @param min_subsample Minimum number of SNPs to subsample in SVD/ALS (50000)
#' @param max_subsample Maximum number of SNPs to subsample in SVD/ALS (500000)
#' @param als_iter the maximum number of iterations in the ALS algorithm (1000)
#' @param als_tole the tolerance for the RMSE of P between iterations (1e-4)
#' @param no_freqs do not save P-matrix (TRUE)
#' @param random_init random initialisation of parameters (TRUE)
#' @param safety deprecated and ignored: fastmixture turns on its safety
#'   updates by itself when an update is unstable
#' @param cv the number of cross-validation folds, greater than 1 (NULL, no
#'   cross-validation)
#' @param cv_tole the tolerance for the cross-validation error
#' in scaled log-likelihood units (1e-7)
#' @return an object of class `gt_admix`. See [tidypopgen::gt_admixture()] for
#'   details.
#' @export

gt_fastmixture <- function(
  x,
  k,
  n_runs = 1,
  threads = 1,
  seed = 42,
  iter = 1000,
  tole = 1e-9,
  batches = 32,
  supervised = NULL,
  check = 5,
  power = 11,
  chunk = 8192,
  subsample = 0.7,
  min_subsample = 50000,
  max_subsample = 500000,
  als_iter = 1000,
  als_tole = 1e-4,
  no_freqs = TRUE,
  random_init = TRUE,
  safety = NULL,
  cv = NULL,
  cv_tole = 1e-7
) {
  if (!is.null(safety)) {
    warning(
      "'safety' is deprecated and ignored: fastmixture turns on its safety ",
      "updates by itself when an update is unstable"
    )
  }
  # check the arguments before starting python (fastmixture needs K > 1)
  check_count(k, "k", min = 2, scalar = FALSE,
    message = "'k' should be a vector of integers greater than 1"
  )
  check_count(n_runs, "n_runs")
  check_count(threads, "threads")
  # one seed per run: the same seeds are reused for every value of k. A
  # single seed (e.g. the default) is used for the first run, and
  # incremented for the following ones
  if (!is.numeric(seed) || length(seed) == 0 || any(!is.finite(seed))) {
    stop("'seed' should only contain numbers")
  }
  if (length(seed) == 1 && n_runs > 1) {
    seed <- seed + seq_len(n_runs) - 1
  }
  if (length(seed) != n_runs) {
    stop("'seed' should be a single number or a vector of length 'n_runs'")
  }

  # the python wrapper reads either a set of PLINK files (from their prefix)
  # or the backing file of a gen_tibble. The latter is the raw .bk file of
  # the FBM, which has no header, so its dimensions are passed as well. The
  # FBM can hold more individuals and loci than the gen_tibble, in a
  # different order: filtering, reordering and select_loci() only change the
  # indices of the gen_tibble into the FBM (its rows for the individuals, in
  # the genotypes column, and its columns for the loci, in `big_index`), so
  # these indices are passed too (0-based, for python)
  fbm_dim <- NULL
  rows <- NULL
  cols <- NULL
  if (inherits(x, "character")) {
    bfile <- x
    plink <- TRUE
  } else if (inherits(x, "gen_tbl")) {
    bfile <- tidypopgen::gt_get_file_names(x)[2]
    fbm <- attr(x$genotypes, "fbm")
    fbm_dim <- as.integer(c(fbm$nrow, fbm$ncol))
    rows <- as.integer(vctrs::vec_data(x$genotypes)) - 1L
    cols <- as.integer(tidypopgen::show_loci(x)$big_index) - 1L
    plink <- FALSE
  } else {
    stop(
      paste0(
        "data must be a gen_tibble, or a character string ",
        "with the prefix of the plink files"
      )
    )
  }

  if (!is.null(cv)) {
    check_count(cv, "cv", min = 2,
      message = "cv must be a single integer greater than 1"
    )
  }

  if (!is.null(supervised) && !is.null(cv)) {
    stop("Cross-validation only works with unsupervised mode")
  }

  # activate the conda environment and import the python wrapper (this
  # fails with an informative error if fastmixture is not installed)
  py_rfastmixture <- .ensure_fastmixture_python()

  # create a namespace object with all the inputs: the python wrapper is
  # derived from the main function of fastmixture, which takes the parsed
  # command line arguments (an argparse Namespace)
  argparse <- reticulate::import("argparse")

  # initialise list to store results
  adm_list <- list(
    k = NULL,
    Q = list()
  )

  if (!is.null(cv)) {
    adm_list$cv <- numeric()
  }

  if (!no_freqs) {
    adm_list$P <- list()
  }

  # Populate adm_list object: one run for each combination of k and seed,
  # with all the runs of a given k stored consecutively (index is the
  # position of the run in the lists of Q and P matrices)
  index <- 1
  for (this_k in as.integer(k)) {
    for (this_rep in seq_len(n_runs)) {
      # the two namespaces only differ in the cross-validation arguments;
      # integer arguments are converted explicitly, as python's argparse
      # would do (R numbers are passed to python as floats)
      if (is.null(cv)) {
        rfastmixture_args <- argparse$Namespace(
          bfile = bfile,
          K = this_k,
          threads = as.integer(threads),
          seed = as.integer(seed[this_rep]),
          iter = as.integer(iter),
          tole = tole,
          batches = as.integer(batches),
          supervised = supervised,
          check = as.integer(check),
          power = as.integer(power),
          chunk = as.integer(chunk),
          subsample = subsample,
          min_subsample = as.integer(min_subsample),
          max_subsample = as.integer(max_subsample),
          als_iter = as.integer(als_iter),
          als_tole = als_tole,
          no_freqs = no_freqs,
          random_init = random_init,
          plink = plink,
          fbm_dim = fbm_dim,
          rows = rows,
          cols = cols,
          projection = NULL,
          cv = NULL,
          cv_tole = NULL
        )
        fastmixture_res <- py_rfastmixture$fastmixture_run(
          args = rfastmixture_args
        )
      } else if (cv > 1) {
        rfastmixture_args <- argparse$Namespace(
          bfile = bfile,
          K = this_k,
          threads = as.integer(threads),
          seed = as.integer(seed[this_rep]),
          iter = as.integer(iter),
          tole = tole,
          batches = as.integer(batches),
          supervised = supervised,
          check = as.integer(check),
          power = as.integer(power),
          chunk = as.integer(chunk),
          subsample = subsample,
          min_subsample = as.integer(min_subsample),
          max_subsample = as.integer(max_subsample),
          als_iter = as.integer(als_iter),
          als_tole = als_tole,
          no_freqs = no_freqs,
          random_init = random_init,
          plink = plink,
          fbm_dim = fbm_dim,
          rows = rows,
          cols = cols,
          projection = NULL,
          cv = as.integer(cv),
          cv_tole = cv_tole
        )
        fastmixture_res <- py_rfastmixture$fastmixture_run(
          args = rfastmixture_args
        )
      }

      # the python function returns Q alone, or a tuple (a list in R) of Q,
      # P (unless no_freqs) and the cross-validation results (if cv, a dict
      # with the average and standard deviation of the error over folds)
      if (no_freqs) {
        if (is.null(cv)) {
          q_matrix <- tidypopgen::q_matrix(fastmixture_res)
          adm_list$Q[[index]] <- q_matrix
          adm_list$k <- sapply(adm_list$Q, ncol)
        } else {
          q_matrix <- tidypopgen::q_matrix(fastmixture_res[[1]])
          adm_list$Q[[index]] <- q_matrix
          adm_list$k <- sapply(adm_list$Q, ncol)
          adm_list$cv <- c(
            adm_list$cv,
            fastmixture_res[[2]]$avg
          )
        }
      } else {
        if (is.null(cv)) {
          names(fastmixture_res) <- c("Q", "P")
          q_matrix <- tidypopgen::q_matrix(fastmixture_res$Q)
          p_matrix <- as.matrix(fastmixture_res$P)
          adm_list$Q[[index]] <- q_matrix
          adm_list$P[[index]] <- p_matrix
          adm_list$k <- sapply(adm_list$Q, ncol)
        } else {
          names(fastmixture_res) <- c("Q", "P", "cv")
          q_matrix <- tidypopgen::q_matrix(fastmixture_res$Q)
          p_matrix <- as.matrix(fastmixture_res$P)
          adm_list$Q[[index]] <- q_matrix
          adm_list$P[[index]] <- p_matrix
          adm_list$k <- sapply(adm_list$Q, ncol)
          adm_list$cv <- c(
            adm_list$cv,
            fastmixture_res$cv$avg
          )
        }
      }
      index <- index + 1
    }
  }

  # return list, as a gt_admix object (like tidypopgen::gt_admixture())
  class(adm_list) <- c("gt_admix", class(adm_list))

  # add metadata if x is a gen_tibble
  if (inherits(x, "gen_tbl")) {
    adm_list$id <- x$id
    # if it is grouped, add the group
    if (inherits(x, "grouped_gen_tbl")) {
      adm_list$group <- x[[dplyr::group_vars(x)]]
    }
  }

  # add info on algorithm
  adm_list$algorithm <- "fastmixture"

  return(adm_list)
}

# Check that fastmixture is installed, activate its conda environment, and
# return the python wrapper module (inst/python/py_rfastmixture.py)
.ensure_fastmixture_python <- function() {
  envname <- tgc_conda_env()
  if (!tgc_tool_installed("fastmixture")) {
    stop(
      "gt_fastmixture() needs the python package 'fastmixture', an optional ",
      "tool that is not installed. Install it with ",
      "tgc_tools_install(tools = \"fastmixture\") or, if the conda ",
      "environment '", envname, "' already exists, with ",
      "tgc_fastmixture_install()"
    )
  }
  # use the python of the conda environment (an error if python has already
  # been started with a different environment in this session)
  reticulate::use_condaenv(envname, required = TRUE)
  # get info on the current python environment
  cfg <- reticulate::py_config()

  if (!reticulate::py_module_available("fastmixture")) {
    stop(
      "Python environment selected but 'fastmixture' is not importable.\n",
      "Current python: ", cfg$python
    )
  }
  # import the python module only now, as fastmixture is an optional tool
  # (importing it when the package is loaded would fail if python had
  # already been started without fastmixture)
  reticulate::import_from_path(
    module = "py_rfastmixture",
    path = system.file("python", package = "tidygenclust")
  )
}

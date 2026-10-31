#' ADAMIXTURE algorithm for population genetics clustering
#'
#' This function runs ADAMIXTURE, a python package for biobank-scale genetic
#' clustering with the model of ADMIXTURE. It combines a spectral
#' initialisation (SVD and ALS) with accelerated block relaxation (BR-QN) to
#' reach high likelihoods with fast convergence. If you use this function, make
#' sure that you cite the paper by Saurina-i-Ricos, Mas Montserrat and
#' Ioannidis (2026).
#'
#' This function returns a `gt_admix` object that can be plotted with
#' `autoplot`, and tidied with `tidy` methods from the `tidypopgen` package.
#' Cross-validation is not performed by default (`cv = NULL`); if you want to
#' include cross-validation you can set `cv` to the number of folds, a value
#' greater than 1 (e.g. 5 for 5-fold cross-validation).
#'
#' `ADAMIXTURE` is an optional tool, which is not needed by the rest of
#' `tidygenclust` (see [gt_radmix()] for a native implementation of the
#' ADMIXTURE model). It is installed by [tgc_tools_install()], or it can be
#' added to an existing `ctidygenclust` conda environment with
#' [tgc_adamixture_install()]. Its CPU version, used here, does not need
#' PyTorch.
#'
#' @references Saurina-i-Ricos, J., Mas Montserrat, D., & Ioannidis, A. G.
#'   (2026). ADAMIXTURE: adaptive first-order optimization for biobank-scale
#'   genetic clustering. Bioinformatics, 42(Supplement_1). doi:
#'   https://doi.org/10.1093/bioinformatics/btag236
#'
#' @param x either a [`tidypopgen::gen_tibble`], or the name of the binary plink
#'   file (without the .bed extension)
#' @param k the number of ancestral components (clusters), either a single value
#'   or a vector
#' @param n_runs the number of repeats for each k value (1)
#' @param threads the number of threads to use (1)
#' @param seed the random seed (defaults to 42): a vector of length
#'   `n_runs`, with one seed per run (the same seeds are used for each value
#'   of `k`), or a single seed, used for the first run and incremented for
#'   the following ones
#' @param tol the tolerance: the fit stops when the log-likelihood changes by
#'   less than `tol` between iterations (0.1)
#' @param max_iter the maximum number of iterations (10000)
#' @param no_freqs do not save P-matrix (TRUE)
#' @param cv the number of cross-validation folds, greater than 1 (NULL, no
#'   cross-validation)
#' @param verbose boolean; should the progress of ADAMIXTURE be printed
#'   (TRUE)?
#' @return an object of class `gt_admix`. See [tidypopgen::gt_admixture()] for
#'   details. As for [gt_radmix()], the P matrices give the frequency of the
#'   alternate allele of each locus (for PLINK files, of the first allele in
#'   the bim file).
#' @export

gt_adamixture <- function(
  x,
  k,
  n_runs = 1,
  threads = 1,
  seed = 42,
  tol = 0.1,
  max_iter = 10000,
  no_freqs = TRUE,
  cv = NULL,
  verbose = TRUE
) {
  # check the arguments before starting python
  check_count(k, "k", scalar = FALSE,
    message = "'k' should be a vector of positive integers"
  )
  check_count(n_runs, "n_runs")
  check_count(threads, "threads")
  check_count(max_iter, "max_iter")
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol <= 0) {
    stop("'tol' should be a single positive number")
  }
  if (!is.null(cv)) {
    check_count(cv, "cv", min = 2,
      message = "cv must be a single integer greater than 1"
    )
  }
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
  group <- NULL
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
    # the group of each individual, stored in the result: it is taken now,
    # so that a gen_tibble grouped by several variables fails before the
    # fits rather than after them
    group <- gt_single_group(x)
  } else {
    stop(
      paste0(
        "data must be a gen_tibble, or a character string ",
        "with the prefix of the plink files"
      )
    )
  }

  # activate the conda environment and import the python wrapper (this
  # fails with an informative error if ADAMIXTURE is not installed)
  py_radamixture <- .ensure_adamixture_python()

  adm_list <- list(
    k = NULL,
    Q = list()
  )
  if (!no_freqs) {
    adm_list$P <- list()
  }
  if (!is.null(cv)) {
    adm_list$cv <- numeric()
  }

  # one run for each combination of k and seed, with all the runs of a given
  # k stored consecutively (index is the position of the run in the lists
  # of Q and P matrices); integers are converted explicitly, as R numbers
  # are passed to python as floats
  index <- 1
  for (this_k in as.integer(k)) {
    for (this_rep in seq_len(n_runs)) {
      res <- py_radamixture$adamixture_fit(
        bfile = bfile,
        plink = plink,
        fbm_dim = fbm_dim,
        rows = rows,
        cols = cols,
        K = this_k,
        seed = as.integer(seed[this_rep]),
        threads = as.integer(threads),
        tol = tol,
        max_iter = as.integer(max_iter),
        cv = if (is.null(cv)) NULL else as.integer(cv),
        verbose = verbose
      )
      adm_list$Q[[index]] <- tidypopgen::q_matrix(res$Q)
      if (!no_freqs) {
        adm_list$P[[index]] <- as.matrix(res$P)
      }
      if (!is.null(cv)) {
        adm_list$cv <- c(adm_list$cv, res$cv)
      }
      index <- index + 1
    }
  }
  adm_list$k <- sapply(adm_list$Q, ncol)

  # return list, as a gt_admix object (like tidypopgen::gt_admixture())
  class(adm_list) <- c("gt_admix", class(adm_list))

  # add metadata if x is a gen_tibble (the group was taken before the fits)
  if (inherits(x, "gen_tbl")) {
    adm_list$id <- x$id
    if (!is.null(group)) {
      adm_list$group <- group
    }
  }

  # add info on algorithm
  adm_list$algorithm <- "adamixture"

  return(adm_list)
}

# Check that ADAMIXTURE is installed, activate its conda environment, and
# return the python wrapper module (inst/python/py_radamixture.py)
.ensure_adamixture_python <- function() {
  envname <- tgc_conda_env()
  if (!tgc_tool_installed("adamixture")) {
    stop(
      "gt_adamixture() needs the python package 'adamixture', an optional ",
      "tool that is not installed. Install it with ",
      "tgc_tools_install(tools = \"adamixture\") or, if the conda ",
      "environment '", envname, "' already exists, with ",
      "tgc_adamixture_install()"
    )
  }
  # use the python of the conda environment (an error if python has already
  # been started with a different environment in this session)
  reticulate::use_condaenv(envname, required = TRUE)
  cfg <- reticulate::py_config()
  if (!reticulate::py_module_available("adamixture")) {
    stop(
      "Python environment selected but 'adamixture' is not importable.\n",
      "Current python: ", cfg$python
    )
  }
  # import the python module only now, as ADAMIXTURE is an optional tool
  reticulate::import_from_path(
    module = "py_radamixture",
    path = system.file("python", package = "tidygenclust")
  )
}

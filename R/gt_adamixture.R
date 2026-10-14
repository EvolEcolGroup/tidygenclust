#' ADAMIXTURE algorithm for population genetics clustering
#'
#' This function implements the ADAMIXTURE algorithm for biobank-scale
#' population genetics clustering by calling the python module.
#'
#' This function returns a q_matrix that can be plotted with `autoplot`, and
#' tidied with `tidy` methods from the `tidypopgen` package. Cross-validation
#' is set to NULL as default, if you want to include cross-validation you can
#' set `cv` to a value greater than 1 (e.g. 5 for 5-fold cross-validation).
#'
#' `ADAMIXTURE` is not installed by default by [tgc_tools_install()]; use
#' [tgc_adamixture_install()] to add it to the `ctidygenclust` environment.
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
#' @param seed the random seed (defaults to 42); can be a single integer or a
#'   vector of length `n_runs`
#' @param tole the tolerance in log-likelihood units between iterations (0.1)
#' @param no_freqs do not save P-matrix (TRUE)
#' @param cv the number of cross-validation folds (NULL disables CV)
#' @return an object of class `gt_admix`. See [tidypopgen::gt_admixture()] for
#'   details.
#' @export

gt_adamixture <- function(
  x,
  k,
  n_runs = 1,
  threads = 1,
  seed = 42,
  tole = 0.1,
  no_freqs = TRUE,
  cv = NULL
) {
  .ensure_adamixture_python()

  if (length(seed) == 1L && n_runs > 1L) {
    seed <- seed + seq_len(n_runs) - 1L
  } else if (length(seed) != n_runs) {
    stop("'seed' should be a vector of length 'n_runs'")
  }

  if (inherits(x, "character")) {
    bfile <- x
    n_indiv <- NULL
    n_loci <- NULL
    plink <- TRUE
  } else if (inherits(x, "gen_tbl")) {
    if (!requireNamespace("tidypopgen", quietly = TRUE)) {
      stop("Package 'tidypopgen' is required when passing a gen_tbl object.")
    }
    bfile <- tidypopgen::gt_get_file_names(x)[2]
    n_indiv <- nrow(x)
    n_loci <- nrow(tidypopgen::show_loci(x))
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
    if (length(cv) != 1L || !is.numeric(cv) || cv <= 1) {
      stop("cv must be a single integer greater than 1")
    }
  }

  argparse <- reticulate::import("argparse")

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

  index <- 1
  for (this_k in as.integer(k)) {
    for (this_rep in seq_len(n_runs)) {
      radamixture_args <- argparse$Namespace(
        bfile = bfile,
        K = this_k,
        threads = as.integer(threads),
        seed = as.integer(seed[this_rep]),
        tole = tole,
        no_freqs = no_freqs,
        plink = plink,
        n_indiv = n_indiv,
        n_loci = n_loci,
        cv = if (is.null(cv)) NULL else as.integer(cv)
      )

      adamixture_res <- .py_radamixture$adamixture_run(
        args = radamixture_args
      )

      if (no_freqs) {
        if (is.null(cv)) {
          q_matrix <- .as_q_matrix(adamixture_res)
          adm_list$Q[[index]] <- q_matrix
        } else {
          q_matrix <- .as_q_matrix(adamixture_res[[1]])
          adm_list$Q[[index]] <- q_matrix
          adm_list$cv <- c(
            adm_list$cv,
            adamixture_res[[2]]$avg
          )
        }
      } else {
        if (is.null(cv)) {
          names(adamixture_res) <- c("Q", "P")
          q_matrix <- .as_q_matrix(adamixture_res$Q)
          p_matrix <- as.matrix(adamixture_res$P)
          adm_list$Q[[index]] <- q_matrix
          adm_list$P[[index]] <- p_matrix
        } else {
          names(adamixture_res) <- c("Q", "P", "cv")
          q_matrix <- .as_q_matrix(adamixture_res$Q)
          p_matrix <- as.matrix(adamixture_res$P)
          adm_list$Q[[index]] <- q_matrix
          adm_list$P[[index]] <- p_matrix
          adm_list$cv <- c(
            adm_list$cv,
            adamixture_res$cv$avg
          )
        }
      }
      adm_list$k <- sapply(adm_list$Q, ncol)
      index <- index + 1
    }
  }

  # return list
  class(adm_list) <- c("gt_admix", class(adm_list))

  # add metadata if x is a gen_tibble
  if (inherits(x, "gen_tbl")) {
    adm_list$id <- x$id
    if (inherits(x, "grouped_gen_tbl")) {
      adm_list$group <- x[[dplyr::group_vars(x)]]
    }
  }

  # add info on algorithm
  adm_list$algorithm <- "adamixture"

  return(adm_list)
}

.as_q_matrix <- function(m) {
  if (requireNamespace("tidypopgen", quietly = TRUE)) {
    tidypopgen::q_matrix(m)
  } else {
    class(m) <- c("q_matrix", "matrix", "array")
    m
  }
}

.ensure_adamixture_python <- function() {
  reticulate::use_condaenv("ctidygenclust", required = TRUE)
  if (!reticulate::py_module_available("adamixture")) {
    cfg <- reticulate::py_config()
    stop(
      "Python environment selected but 'adamixture' is not importable.\n",
      "Current python: ", cfg$python, "\n",
      "Install it with 'tgc_adamixture_install()'"
    )
  }
}

#' Run clumppling
#'
#' This function aligns multiple runs of a clustering algorithm (e.g.
#' ADMIXTURE or STRUCTURE) with the method of clumppling (Liu, Kopelman and
#' Rosenberg 2024). Its R and C++ code is ported, with modifications, from the
#' python package clumppling (version 2.0.10, by Xiran Liu;
#' <https://github.com/PopGenClustering/Clumppling>), and gives the same
#' results without the need for Python. The main modifications are that the
#' alignment of each pair of runs, an integer linear program in clumppling,
#' is solved exactly as an equivalent assignment problem, and that community
#' detection uses ports of the implementations that clumppling uses through
#' cdlib (python-louvain, in C++, and markov_clustering, in R) or the
#' `igraph` package; see `data-raw/clumppling/README.md` in the source of the
#' package for all the differences. If you use this function, make sure that
#' you cite the paper on clumppling by Liu, Kopelman, and Rosenberg (2024).
#'
#' The example data (`capeverde_admix.rds` and `chicken_admix.rds` in
#' `extdata`) are the example datasets of clumppling: runs of ADMIXTURE on the
#' Cape Verde data of Verdu et al. (2017), and runs of STRUCTURE on the
#' chicken data of Rosenberg et al. (2001).
#'
#' The algorithm has three steps:
#' 1. all pairs of replicates with the same K are aligned, finding the
#'    matching of clusters that minimises the squared distance between their
#'    memberships;
#' 2. replicates with the same K are grouped into modes by community
#'    detection on the network of alignment costs; each mode is summarised by
#'    a representative replicate (the one with the smallest cost to the
#'    others) or by the average of its replicates;
#' 3. modes are aligned across K (merging two clusters of the larger K when
#'    two values of K differ by one), and their clusters are reordered so that
#'    they are consistent across K.
#'
#' The replicates are the Q matrices of a `gt_admix` object, as returned by
#' the clustering functions of 'tidypopgen' (e.g.
#' `tidypopgen::gt_admixture()`) and 'tidygenclust' (e.g.
#' [gt_radmix()]); to align replicates
#' stored as Q files, read them into a `gt_admix` object with
#' `tidypopgen::read_q_files()`. Replicates with the same K are numbered in the
#' order in which they appear in the `gt_admix` object. The rows of each Q
#' matrix must sum to 1 (to within 1e-6), so memberships that were rounded
#' (e.g. by STRUCTURE) should be rescaled first.
#'
#' If you would like to generate an annotated autoplot from your
#' `gt_clumppling` object, ensure that all individuals from the same
#' population are adjacent to one another in the `gt_admix` object supplied to
#' `gt_clumppling` (see `tidypopgen::gt_admix_reorder_q()`). Autoplot 'group'
#' argument requires that all individuals from the same group are adjacent.
#'
#' The only random step is community detection. For the default 'louvain'
#' method, the C++ implementation follows python-louvain (used by
#' clumppling), but random numbers are drawn from R's generator: the modes
#' are the same as in python unless the community structure is ambiguous, but
#' the same seed does not give the same random choices in R and python. If
#' `seed` is given, the state of R's random number generator is restored on
#' exit. When a mode includes replicates with identical Q matrices, the
#' representative of the mode can also differ from the python implementation
#' (it is then an identical replicate), as their costs only differ by
#' rounding errors.
#'
#' @references X. Liu, N. M. Kopelman, N. A. Rosenberg (2024) Clumppling:
#'   cluster matching and permutation program with integer linear programming.
#'   Bioinformatics 40(1): btad751. doi:
#'   https://doi.org/10.1093/bioinformatics/btad751
#'
#'   N. A. Rosenberg, T. Burke, K. Elo, M. W. Feldman, P. J. Freidlin,
#'   M. A. M. Groenen, J. Hillel, A. Mäki-Tanila, M. Tixier-Boichard,
#'   A. Vignal, K. Wimmers, S. Weigend (2001) Empirical evaluation of genetic
#'   clustering methods using multilocus genotypes from 20 chicken breeds.
#'   Genetics 159(2): 699-713.
#'
#'   P. Verdu, E. M. Jewett, T. J. Pemberton, N. A. Rosenberg, M. Baptista
#'   (2017) Parallel trajectories of genetic and linguistic admixture in a
#'   genetically admixed creole population. Current Biology 27(16):
#'   2529-2535.
#'
#' @param x a `gt_admix` object with at least two Q matrices
#' @param use_rep boolean, whether to use representative modes (alternative:
#'  average), defaults to TRUE
#' @param merge boolean, whether to merge two clusters when aligning K+1 to K,
#'  defaults to TRUE
#' @param cd_method the community detection method to use, one of 'louvain'
#'   (default), 'leiden', 'infomap', 'markov_clustering', 'label_propagation',
#'   'walktrap', 'custom'. 'leiden', 'infomap', 'label_propagation' and
#'   'walktrap' require the 'igraph' package; as in clumppling (which calls
#'   them through 'cdlib'), only 'louvain', 'infomap' and 'markov_clustering'
#'   use the similarity between replicates as edge weights, while 'leiden',
#'   'label_propagation' and 'walktrap' only use the presence of edges (all
#'   pairs of replicates are connected, except for the most dissimilar
#'   pair). 'custom' reproduces the example
#'   custom method of clumppling (Markov clustering with a pruning threshold
#'   of 0.2); alternatively, `cd_method` can be a function that takes the
#'   adjacency matrix between replicates (similarities rescaled to \[0, 1\],
#'   with ones on the diagonal) and returns a vector of community labels.
#' @param use_best_pair boolean, whether to use best pair as anchor for across-K
#'  alignment (alternative: major), defaults to TRUE
#' @param seed (optional) an integer used to set the random seed, so that
#'   community detection (e.g. 'louvain') is reproducible. Defaults to NULL
#'   (no seed is set).
#' @param output_path deprecated and ignored. Results are returned as R objects
#'   rather than written to files.
#' @param cd_res the resolution of the 'louvain' community detection method,
#'   defaults to 1
#' @param test_comm boolean, whether to test for community structure (with a
#'   Tracy-Widom test) before running community detection; if there is no
#'   structure, all replicates form a single mode. Defaults to FALSE.
#' @param comm_min,comm_max cost thresholds for community detection: if all
#'   pairwise costs between replicates are below `comm_min`, they form a
#'   single mode; if they are all above `comm_max`, each replicate is a mode.
#'   Defaults to 1e-6 and 1e-2.
#' @returns a list of class `gt_clumppling` containing:
#' - aligned_modes: a list of data.frames, each data.frame is a Q-matrix
#' - cost_acrossK: a named list of costs for each pairwise K alignment
#' - mode_replicates: a list of replicate indices for each mode
#' - K_range: vector of K values analyzed
#' - N: number of individuals
#' - mode_stats: a data.frame with the size, cost and performance (average
#'   similarity among replicates) of each mode
#' @export
#' @examples
#' capeverde <- readRDS(
#'   system.file("extdata", "capeverde_admix.rds", package = "tidygenclust")
#' )
#' clump_res <- gt_clumppling(capeverde, seed = 1)
#' clump_res$mode_stats
gt_clumppling <- function(
  x,
  use_rep = TRUE,
  merge = TRUE,
  cd_method = "louvain",
  use_best_pair = TRUE,
  seed = NULL,
  output_path = NULL,
  cd_res = 1,
  test_comm = FALSE,
  comm_min = 1e-6,
  comm_max = 1e-2
) {
  if (!is.null(output_path)) {
    warning(
      "'output_path' is deprecated and ignored: clumppling results are ",
      "returned as R objects and are no longer written to files"
    )
  }
  # cd_method is either the name of a supported method or a user function
  # (see clumppling_cd_default())
  is_method_name <- is.character(cd_method) && length(cd_method) == 1 &&
    cd_method %in% clumppling_cd_methods
  valid_method <- is.function(cd_method) || is_method_name
  if (!valid_method) {
    stop(
      "'cd_method' must be one of ",
      paste(clumppling_cd_methods, collapse = ", "), ", or a function"
    )
  }
  # the other options are single values: NA or non-finite values would make
  # community detection fail in obscure ways (a non-finite resolution made
  # louvain loop forever)
  for (arg in c("use_rep", "merge", "use_best_pair", "test_comm")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1 || is.na(value)) {
      stop("'", arg, "' must be TRUE or FALSE")
    }
  }
  is_number <- function(value) {
    is.numeric(value) && length(value) == 1 && is.finite(value)
  }
  if (!is_number(cd_res) || cd_res <= 0) {
    stop("'cd_res' must be a single positive number")
  }
  thresholds_ok <- is_number(comm_min) && is_number(comm_max) &&
    comm_min >= 0 && comm_max >= 0
  if (!thresholds_ok) {
    stop("'comm_min' and 'comm_max' must be single non-negative numbers")
  }

  if (!inherits(x, "gt_admix")) {
    stop(
      "'x' must be a gt_admix object (Q files can be read with ",
      "tidypopgen::read_q_files())"
    )
  }
  if (length(x$Q) < 2) {
    stop("'x' must contain more than one Q matrix")
  }
  # validate and rescale the Q matrices, drop those with K = 1, and name them
  # as clumppling does (e.g. "3_K2R1"), grouped by K
  q_data <- clumppling_prepare_q_list(x$Q)

  # the seed is set only now (and the RNG state restored when this function
  # exits), so that it only affects community detection
  if (!is.null(seed)) {
    clumppling_local_seed(seed)
  }
  clump_run <- clumppling_run(q_data,
    cd_method = cd_method, cd_res = cd_res, test_comm = test_comm,
    comm_min = comm_min, comm_max = comm_max, merge = merge,
    use_rep = use_rep, use_best_pair = use_best_pair
  )

  # now create the output object, with the same structure as the one
  # returned by previous versions of tidygenclust (which called python)
  clump_res <- list()
  clump_res$aligned_modes <- lapply(clump_run$aligned_modes, as.data.frame)
  names(clump_res$aligned_modes) <- clump_run$mode_names
  clump_res$cost_acrossK <- clump_run$cost_acrossK

  # for each mode, get the replicates (the replicate id is the number after R
  # in the name of the replicate, e.g. 3_K2R1)
  mode_alignment <- clump_run$mode_alignment
  replicate_id <- as.integer(sub(".*R(\\d+)$", "\\1", mode_alignment$Replicate))
  clump_res$mode_replicates <- split(
    replicate_id,
    factor(mode_alignment$Mode, levels = clump_run$mode_names)
  )

  clump_res$K_range <- clump_run$K_range
  clump_res$N <- as.integer(clump_run$N)
  clump_res$mode_stats <- clump_run$mode_stats

  class(clump_res) <- "gt_clumppling"
  return(clump_res)
}

#' Set the random seed until the calling function exits
#'
#' Sets the seed of R's random number generator, and restores its previous
#' state (or absence of state) when the calling function exits.
#' @param seed an integer
#' @param envir the frame of the calling function
#' @returns NULL (invisibly)
#' @keywords internal
#' @noRd
clumppling_local_seed <- function(seed, envir = parent.frame()) {
  global_env <- globalenv()
  old_seed <- global_env[[".Random.seed"]]
  restore <- function() {
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = global_env, inherits = FALSE)) {
        rm(".Random.seed", envir = global_env)
      }
    } else {
      global_env[[".Random.seed"]] <- old_seed
    }
  }
  # on.exit() registers its expression in the frame it is evaluated in, so
  # it is called (with do.call()) in the calling frame, as withr does. The
  # call embeds the function itself, as `restore` is not visible from the
  # calling frame
  do.call(on.exit, list(as.call(list(restore)), add = TRUE), envir = envir)
  set.seed(seed)
  invisible(NULL)
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
#' @param k a vector of values of K to subset to: consecutive values of
#'   `K_range` (e.g. `c(17, 19)` if `K_range` is 17, 19 and 21), in increasing
#'   order
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
    # (indexing beyond the individuals would silently give rows of NA)
    if (min(indivs) < 1 || max(indivs) > x$N) {
      stop("indivs values must be between 1 and N (", x$N, ")")
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
    # k must be values of K_range that follow each other in it, in increasing
    # order (with gaps in K_range, e.g. 17, 19, 21, k = c(17, 19) is
    # consecutive), so that the modes that are kept are aligned to each other
    k_pos <- match(k, x$K_range)
    if (length(k) == 0 || anyNA(k_pos)) {
      stop(
        "k values must be values of K in K_range (",
        paste(x$K_range, collapse = ", "), ")"
      )
    }
    if (any(diff(k_pos) != 1)) {
      stop(
        "k values must be consecutive from min to max, among the values of ",
        "K in K_range"
      )
    }
    # the K of each mode, from its name (e.g. 12 for "K12M3"); selecting modes
    # by matching the text "K1" would also select K10, K11, ...
    mode_k <- function(mode) as.integer(sub("^K(\\d+)M\\d+$", "\\1", mode))
    keep_modes <- mode_k(names(x$aligned_modes)) %in% k
    x$K_range <- x$K_range[x$K_range %in% k]
    x$aligned_modes <- x$aligned_modes[keep_modes]
    x$mode_replicates <- x$mode_replicates[keep_modes]
    if (!is.null(x$mode_stats)) {
      x$mode_stats <- x$mode_stats[mode_k(x$mode_stats$Mode) %in% k, ,
        drop = FALSE
      ]
    }
    # costs are named after the two modes they compare (e.g. "K4M1-K5M2"):
    # keep those between modes that are kept
    cost_modes <- strsplit(names(x$cost_acrossK), "-", fixed = TRUE)
    keep_costs <- vapply(cost_modes, function(m) all(mode_k(m) %in% k),
      logical(1)
    )
    x$cost_acrossK <- x$cost_acrossK[keep_costs]
  }

  return(x)
}

#' Estimate ancestry proportions with the ADMIXTURE method
#'
#' This function estimates individual ancestry proportions (Q) and the allele
#' frequencies of the ancestral populations (P) with the method of ADMIXTURE
#' (Alexander et al. 2009), i.e. by maximum likelihood under the ADMIXTURE
#' model. It uses a modified version of the algorithm of OpenADMIXTURE (Ko et
#' al. 2023): EM priming, block-wise Newton updates solved as constrained
#' quadratic programs, and quasi-Newton acceleration. Much of its C++ code is
#' ported from the Julia implementation of OpenADMIXTURE,
#' [OpenADMIXTURE.jl](https://github.com/OpenMendel/OpenADMIXTURE.jl), by Ko
#' et al. The function has a different name because it differs from
#' OpenADMIXTURE: it does not include the selection of ancestry informative
#' markers (sparse K-means) that is a central part of OpenADMIXTURE, it
#' supports pseudohaploid data and cross-validation, and it changes some
#' details of the computations (see "Differences from OpenADMIXTURE.jl"
#' below). If you use this function, please cite Ko et al. (2023) for the
#' algorithm, as well as Alexander et al. (2009) for the method.
#'
#' Genotypes are read directly from the backing file of the `gen_tibble`
#' (respecting any filtering of individuals and selection of loci), so no
#' intermediate PLINK file is written. Computation is multithreaded with
#' OpenMP.
#'
#' @details Missing genotypes are ignored in the likelihood. Pseudohaploid
#'   individuals (see [tidypopgen::gt_pseudohaploid()]) contribute a single
#'   allele per locus. Only diploid and pseudohaploid `gen_tibble`s are
#'   supported.
#'
#'   `P` gives the frequency of the **alternate** allele (`allele_alt` in
#'   [tidypopgen::show_loci()]) in each ancestral population.
#'
#'   Supplying `init_p` together with `fix_p = TRUE` estimates the ancestry of
#'   new individuals given the allele frequencies of a reference panel (a
#'   "projection" analysis); the loci of `x` must then match the rows of
#'   `init_p`.
#'
#'   ## Algorithm
#'
#'   As in OpenADMIXTURE.jl, each run consists of:
#'   1. a random initialisation followed by a few EM (FRAPPE) iterations
#'      (`n_em_iter`);
#'   2. block-wise Newton updates of Q (one individual at a time, on the
#'      simplex) and P (one locus at a time, in the unit interval), each
#'      solved exactly as a small quadratic program with K. Lange's
#'      sweep-operator active-set method;
#'   3. quasi-Newton acceleration (Zhou, Alexander and Lange 2011) using the
#'      last `n_qn` secant pairs, accepted only when it improves the
#'      log-likelihood.
#'
#'   Iterations stop when the relative change in log-likelihood falls below
#'   `rtol`. All estimates are kept within \[1e-5, 1 - 1e-5\], as in ADMIXTURE
#'   and OpenADMIXTURE.jl.
#'
#'   ## Cross-validation
#'
#'   With `crossval = TRUE`, the cross-validation error of ADMIXTURE
#'   (Alexander and Lange 2011) is computed for each run, and can be used to
#'   choose `k` (lower is better; see `autoplot(type = "cv")`). The observed
#'   genotypes are split into `n_folds` folds of roughly equal size; for each
#'   fold, the model is refitted with the genotypes of that fold treated as
#'   missing, and each of them is predicted by its expected value
#'   \eqn{\mu = n \sum_k q_{ik} p_{kj}}{mu = n sum_k q_ik p_kj}, where
#'   \eqn{n} is the number of allele copies of the individual. The error is
#'   the binomial deviance
#'   \deqn{d(g, \mu) = g \log(g / \mu) + (n - g) \log((n - g) / (n - \mu))}{
#'   d(g, mu) = g log(g / mu) + (n - g) log((n - g) / (n - mu))}
#'   averaged over all the masked genotypes of all folds. Where ADMIXTURE does
#'   not document the details, we made the following choices:
#'   - each fold is fitted starting from the estimates for the full data, with
#'     the same convergence criteria (`max_iter` and `rtol`);
#'   - genotypes are assigned to folds by a hash of their position and of a
#'     seed drawn from the random number generator (reproducible with
#'     `seed`), so the folds differ from those of ADMIXTURE and the errors
#'     agree with ADMIXTURE's closely but not exactly;
#'   - pseudohaploid individuals use \eqn{n = 1}.
#'
#'   Cross-validation does not change the estimates of Q and P, but requires
#'   `n_folds` additional fits for each run (each usually converging in a few
#'   iterations).
#'
#'   ## Comparison with ADMIXTURE
#'
#'   The model and the log-likelihood are those of ADMIXTURE (Alexander et al.
#'   2009), and the optimisation is closely related (block updates with
#'   constrained Newton steps and quasi-Newton acceleration), so both programs
#'   reach the same maximum likelihood estimates, up to local optima and the
#'   convergence tolerance. The test suite checks this on the `EUR_subset`
#'   data used in the OpenADMIXTURE.jl tests. Note that:
#'   - ADMIXTURE reports the frequency of the second allele in the PLINK bim
#'     file, which for a `gen_tibble` is the reference allele, so its P
#'     matrices are 1 - P of those returned here;
#'   - cross-validation follows ADMIXTURE (see above), and the test suite
#'     checks that the errors agree with those of ADMIXTURE on `EUR_subset`;
#'   - bootstrap standard errors and supervised analyses with individuals of
#'     known ancestry, available in ADMIXTURE, are not implemented (a fixed
#'     reference panel can be used with `init_p` and `fix_p`).
#'
#'   ## Differences from OpenADMIXTURE.jl
#'
#'   `gt_radmix()` uses a modified version of the algorithm of
#'   OpenADMIXTURE.jl. Starting from the same values, it reproduces the
#'   computations of OpenADMIXTURE.jl to numerical precision (this is checked
#'   by the test suite), except for the following:
#'   - there is no GPU (CUDA) support;
#'   - cross-validation, which OpenADMIXTURE.jl does not offer, is available
#'     (see above);
#'   - sparse K-means selection of ancestry informative markers (SKFR) is not
#'     implemented; loci can be selected beforehand with
#'     [tidypopgen::select_loci()];
#'   - pseudohaploid individuals are supported (they contribute one allele
#'     per locus, whereas OpenADMIXTURE.jl assumes every individual is
#'     diploid);
#'   - with missing data, the EM update of Q is normalised by the number of
#'     loci observed in each individual (OpenADMIXTURE.jl divides by the
#'     total number of loci, so that Q is pushed off the simplex);
#'   - the log-likelihood is only computed when the algorithm needs it: not
#'     after each EM step (unless `verbose = TRUE`), and, in the main
#'     iterations, the log-likelihood after two plain updates only when the
#'     quasi-Newton step is rejected (OpenADMIXTURE.jl computes both every
#'     time). This saves about one pass over the genotypes in seven, without
#'     changing the estimates; the log then reports "basic = not computed";
#'   - a few numerical safeguards that the Julia code lacks: the conditioning
#'     penalty of the quadratic program for Q is never negative, the
#'     active-set solver never pivots on a vanishing element, Newton steps
#'     that are not finite are discarded (the individual or locus keeps its
#'     current value for that iteration), and the projections onto the
#'     feasible set tolerate non-finite values. None of these changes the
#'     results when the computations are well behaved; each is documented in
#'     the C++ source.
#'
#' @references Ko S, Chu BB, Peterson D, Okenwa C, Papp JC, Alexander DH,
#'   Sobel EM, Zhou H, Lange KL (2023) Unsupervised discovery of ancestry
#'   informative markers and genetic admixture proportions in biobank-scale
#'   datasets. The American Journal of Human Genetics, 110(2), 314-325.
#'
#'   Alexander DH, Novembre J, Lange K (2009) Fast model-based estimation of
#'   ancestry in unrelated individuals. Genome Research, 19(9), 1655-1664.
#'
#'   Alexander DH, Lange K (2011) Enhancements to the ADMIXTURE algorithm for
#'   individual ancestry estimation. BMC Bioinformatics, 12, 246.
#'
#'   Zhou H, Alexander D, Lange K (2011) A quasi-Newton acceleration for
#'   high-dimensional optimization algorithms. Statistics and Computing,
#'   21(2), 261-273.
#'
#' @param x a [`tidypopgen::gen_tibble`] (possibly grouped)
#' @param k an integer vector giving the number(s) of ancestral populations
#' @param n_runs the number of runs for each value of `k` (defaults to 1)
#' @param seed the seed for the random number generator (defaults to `NULL`).
#'   Either a vector of length `n_runs` (the same seeds are used for each `k`)
#'   or of length `n_runs * length(k)`. The global random number generator
#'   state is restored on exit.
#' @param max_iter the maximum number of iterations (defaults to 1000)
#' @param rtol the convergence threshold on the relative change in
#'   log-likelihood between iterations (defaults to 1e-7)
#' @param n_em_iter the number of EM iterations used to prime the estimates
#'   (defaults to 5)
#' @param n_qn the number of secant pairs used by the quasi-Newton
#'   acceleration (defaults to 3)
#' @param init_q optional matrix (individuals x k) of starting ancestry
#'   proportions; only allowed with a single value of `k`
#' @param init_p optional matrix (loci x k) of starting allele frequencies
#'   (of the alternate allele, as in the P matrices returned by this
#'   function; those of [tidypopgen::gt_admixture()] give the frequency of
#'   the reference allele, so use `1 - P`); only allowed with a single value
#'   of `k`. If it has row names (as the P matrices returned by this
#'   function), they are used to match its rows to the loci of `x` (it can
#'   then include other loci, which are ignored); otherwise its rows must be
#'   in the order of the loci of `x`
#' @param fix_q,fix_p booleans; if `TRUE`, Q (respectively P) is kept at the
#'   value given in `init_q` (respectively `init_p`)
#' @param use_imputed boolean; if `TRUE`, imputed genotypes are used as if
#'   they were observed; otherwise they are treated as missing (the default)
#' @param crossval boolean; should the cross-validation error be computed
#'   (defaults to `FALSE`)? See the details.
#' @param n_folds the number of folds used for cross-validation (defaults to
#'   5, as in ADMIXTURE)
#' @param n_cores the number of cores to use (defaults to 1)
#' @param verbose boolean; should progress be printed (defaults to `FALSE`)?
#'   With `verbose = TRUE`, the log-likelihood is also computed (and logged)
#'   every 5 EM steps during the initial EM iterations; otherwise only their
#'   number is logged, as these values are not needed by the algorithm.
#' @returns an object of class `gt_admix`, a list with elements:
#' - `k` the number of clusters of each run
#' - `Q` a list of [tidypopgen::q_matrix()] objects with the ancestry
#'   proportions
#' - `P` a list of matrices (loci x k) with the alternate allele frequencies
#' - `log` a list with the log of each run
#' - `loglik` the log-likelihood of each run
#' - `converged` whether each run converged
#' - `cv` the cross-validation error of each run (if `crossval` is `TRUE`)
#' - `id` the id column of `x`
#' - `group` the group column of `x` (if `x` is grouped)
#' - `algorithm` "radmix"
#' @export
#' @examples
#' bed_file <-
#'   system.file("extdata", "pop_a.bed", package = "tidypopgen")
#' pop_a <- tidypopgen::gen_tibble(bed_file,
#'   backingfile = tempfile("pop_a"),
#'   quiet = TRUE
#' )
#' adm <- gt_radmix(pop_a, k = 2, seed = 1)
#' adm$Q[[1]]
#'
#' # choose k by cross-validation
#' adm_cv <- gt_radmix(pop_a, k = 1:3, seed = 1, crossval = TRUE)
#' adm_cv$cv
gt_radmix <- function(x, k, n_runs = 1, seed = NULL, max_iter = 1000,
                      rtol = 1e-7, n_em_iter = 5, n_qn = 3,
                      init_q = NULL, init_p = NULL, fix_q = FALSE,
                      fix_p = FALSE, use_imputed = FALSE,
                      crossval = FALSE, n_folds = 5, n_cores = 1,
                      verbose = FALSE) {
  # This function is a thin wrapper: it checks the arguments, prepares the
  # data once, and then calls radmix_fit() (R/radmix_fit.R, which
  # runs the fit in C++, src/radmix_fit.cpp) for each run,
  # collecting the results into a gt_admix object with the same structure as
  # the one returned by tidypopgen::gt_admixture().

  # Check the arguments: counts must be whole numbers (rather than being
  # truncated, e.g. k = 2.5), flags TRUE or FALSE (NA would be taken as TRUE
  # by the C++ code), and rtol a positive number (NA or a negative value
  # would silently stop the fit from converging)
  check_count(k, "k", scalar = FALSE,
    message = "'k' should be a vector of positive integers"
  )
  check_count(n_runs, "n_runs",
    message = "'n_runs' should be a positive integer"
  )
  check_count(max_iter, "max_iter", min = 0)
  check_count(n_em_iter, "n_em_iter", min = 0)
  check_count(n_qn, "n_qn")
  check_count(n_cores, "n_cores")
  check_count(n_folds, "n_folds", min = 2,
    message = "'n_folds' should be an integer of at least 2"
  )
  for (flag in c("fix_q", "fix_p", "use_imputed", "crossval", "verbose")) {
    check_flag(get(flag), flag)
  }
  rtol_ok <- is.numeric(rtol) && length(rtol) == 1 && is.finite(rtol) &&
    rtol > 0
  if (!rtol_ok) {
    stop("'rtol' should be a positive number", call. = FALSE)
  }
  k <- as.integer(k)
  n_runs <- as.integer(n_runs)
  n_folds <- as.integer(n_folds)
  n_seeds_ok <- length(seed) %in% c(n_runs, n_runs * length(k))
  if (!is.null(seed) && !n_seeds_ok) {
    stop(
      "'seed' should be a vector of length 'n_runs' OR 'n_runs' * ",
      "'length(k)'",
      call. = FALSE
    )
  }
  # seeds are given to set.seed() for each run: check them now, rather than
  # after the first fits
  if (!is.null(seed) && (!is.numeric(seed) || any(!is.finite(seed)))) {
    stop("'seed' should only contain numbers", call. = FALSE)
  }
  # as in gt_admixture(), seeds can be given per run (and reused for each k)
  # or per combination of k and run
  if (!is.null(seed) && length(seed) == n_runs) {
    seed <- rep(seed, length(k))
  }
  if ((!is.null(init_q) || !is.null(init_p)) && length(k) > 1) {
    stop("'init_q' and 'init_p' can only be used with a single 'k'",
      call. = FALSE
    )
  }
  # fixing both blocks would leave nothing to estimate (and make the
  # quasi-Newton system singular); fixing a block without starting values
  # would fix it at a random value
  if (fix_q && fix_p) {
    stop("'fix_q' and 'fix_p' cannot both be TRUE", call. = FALSE)
  }
  if (fix_q && is.null(init_q)) {
    stop("'fix_q = TRUE' requires 'init_q'", call. = FALSE)
  }
  if (fix_p && is.null(init_p)) {
    stop("'fix_p = TRUE' requires 'init_p'", call. = FALSE)
  }
  # cross-validation refits both Q and P for each fold, so it is not
  # available when one of them is fixed
  if (crossval && (fix_q || fix_p)) {
    stop("'crossval' cannot be used with 'fix_q' or 'fix_p'", call. = FALSE)
  }

  # Prepare the data: indices into the backing file, decoding table and
  # ploidy of each individual (R/radmix_data.R). This is done once and
  # shared by all runs.
  data <- radmix_prepare(x, use_imputed = use_imputed)
  # with as many clusters as individuals, each individual could form its own
  # ancestral population, so the model would not be identifiable
  # (with P fixed, the ancestry of each individual is estimated on its own,
  # so a few individuals, or a single one, can be projected onto a panel)
  if (!fix_p && any(k >= data$n_indiv)) {
    stop("'k' should be smaller than the number of individuals",
      call. = FALSE
    )
  }
  # The starting values can only be checked against the data now that the
  # numbers of individuals and loci are known. They are only allowed with a
  # single k (checked above), hence k[1]. The validated matrices are plain
  # double matrices, moved strictly inside the feasible set, and are shared
  # by all runs.
  init_q <- radmix_check_init(init_q, data$n_indiv, k[1], "init_q", TRUE)
  # a P with row names (e.g. estimated on another dataset, whose loci may be
  # in a different order, or include other loci) is matched to the loci by
  # name: matching by position would silently give wrong results
  if (!is.null(init_p) && !is.null(rownames(init_p))) {
    p_rows <- match(tidypopgen::show_loci(x)$name, rownames(init_p))
    if (anyNA(p_rows)) {
      stop(
        "the row names of 'init_p' do not include all the loci of 'x' ",
        "(remove them to match the rows of 'init_p' to the loci by position)",
        call. = FALSE
      )
    }
    init_p <- init_p[p_rows, , drop = FALSE]
  }
  init_p <- radmix_check_init(init_p, data$n_loci, k[1], "init_p", FALSE)

  # use our own seeds without disturbing the global RNG state
  # (the random starting values are drawn with R's generator, so set.seed()
  # makes each run reproducible; on exit we restore the user's RNG state, or
  # remove .Random.seed if there was none before). The handler is registered
  # before the first set.seed(), so the state is also restored if a fit
  # fails. Without seeds, the fits draw from (and advance) the global
  # generator, like any other R function.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv())) # nolint
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()))
    }
  }

  # The group of each individual, stored in the result. It is taken now, so
  # that a gen_tibble grouped by several variables fails before the fits
  # rather than after them
  group <- gt_single_group(x)

  # Run the fits, one per combination of k and run, filling in the elements
  # of the gt_admix object in the same order as gt_admixture() (by k, then by
  # run)
  # show_loci() lists the loci in the same order as the columns prepared
  # above (`data$cols`), so these names label the rows of P
  loci_names <- tidypopgen::show_loci(x)$name
  adm_list <- list(
    k = integer(0), Q = list(), P = list(), log = list(),
    loglik = numeric(0), converged = logical(0)
  )
  if (crossval) adm_list$cv <- numeric(0)
  # position of the current combination of k and run, both in the output
  # and in `seed` (which was expanded above to one seed per combination)
  index <- 1
  for (this_k in k) {
    for (this_rep in seq_len(n_runs)) {
      if (!is.null(seed)) set.seed(seed[index])
      if (verbose) {
        message("Running k = ", this_k, ", run ", this_rep)
      }
      # Each run draws its own random start inside the C++ code. When both
      # init_q and init_p are given, that start is overwritten and EM
      # priming is skipped, so all the runs give the same result.
      fit <- radmix_fit(
        data,
        k = this_k, init_q = init_q, init_p = init_p,
        fix_q = fix_q, fix_p = fix_p, max_iter = max_iter, rtol = rtol,
        n_em_iter = n_em_iter, n_qn = n_qn, n_cores = n_cores,
        verbose = verbose
      )
      # P has one row per locus, in the order of show_loci(x), and gives the
      # frequency of the alternate allele in each ancestral population
      p_mat <- fit$p
      dimnames(p_mat) <- list(loci_names, paste0(".P", seq_len(this_k)))
      adm_list$k[index] <- this_k
      # q_matrix() gives Q the class and column names (.Q1, .Q2, ...) that
      # tidypopgen's methods (autoplot(), tidy(), augment()) expect
      adm_list$Q[[index]] <- tidypopgen::q_matrix(fit$q)
      adm_list$P[[index]] <- p_mat
      adm_list$log[[index]] <- fit$log
      adm_list$loglik[index] <- fit$loglik
      adm_list$converged[index] <- fit$converged
      if (crossval) {
        # The folds of each run get their own seed, drawn from R's random
        # number generator after the main fit (so that the main fit, and
        # hence Q and P, are the same with or without cross-validation), and
        # therefore reproducible with `seed`. The seed is an integer in
        # [0, 2^31), stored exactly as a double, as radmix_set_cv()
        # expects. Cross-validation also draws random numbers in the fits of
        # the folds, so the state of the generator is restored afterwards:
        # otherwise, without `seed`, it would shift the random starts of the
        # following runs, and change their estimates.
        rng_state <- get(".Random.seed", envir = globalenv())
        cv_seed <- floor(stats::runif(1) * 2^31)
        cv_res <- radmix_cv(data,
          fit = fit, k = this_k, n_folds = n_folds, cv_seed = cv_seed,
          max_iter = max_iter, rtol = rtol, n_qn = n_qn, n_cores = n_cores
        )
        assign(".Random.seed", rng_state, envir = globalenv()) # nolint
        adm_list$cv[index] <- cv_res$cv
        # the log of the fit is followed by one line per fold and by the
        # "CV error (K=...)" line, in the format of ADMIXTURE's output
        adm_list$log[[index]] <- c(fit$log, cv_res$log)
        if (verbose) message(utils::tail(cv_res$log, 1))
      }
      index <- index + 1
    }
  }

  # metadata used by tidypopgen's methods to label and group individuals in
  # plots and tables (the group was taken from x before the fits); both
  # follow the rows of x, as do the rows of Q
  adm_list$id <- x$id
  if (!is.null(group)) {
    adm_list$group <- group
  }
  adm_list$algorithm <- "radmix"
  class(adm_list) <- c("gt_admix", "list")
  adm_list
}

#' Validate a matrix of starting values
#' @param init the matrix (or `NULL`)
#' @param n_rows the expected number of rows
#' @param k the expected number of columns
#' @param arg the name of the argument, for error messages
#' @param simplex boolean; should rows sum to one
#' @returns a numeric matrix, or `NULL`
#' @keywords internal
#' @noRd
radmix_check_init <- function(init, n_rows, k, arg, simplex) {
  if (is.null(init)) {
    return(NULL)
  }
  # as.matrix() also accepts data frames, q_matrix objects (e.g. a Q matrix
  # from a previous run) and, for k = 1, plain vectors; non-numeric values
  # become NA here and are rejected below
  init <- as.matrix(init)
  storage.mode(init) <- "double"
  if (nrow(init) != n_rows || ncol(init) != k) {
    stop("'", arg, "' should be a matrix with ", n_rows, " rows and ", k,
      " columns",
      call. = FALSE
    )
  }
  if (any(!is.finite(init)) || any(init < 0) || any(init > 1)) {
    stop("'", arg, "' should only contain values between 0 and 1",
      call. = FALSE
    )
  }
  if (simplex && any(abs(rowSums(init) - 1) > 1e-6)) {
    stop("the rows of '", arg, "' should sum to 1", call. = FALSE)
  }
  # keep starting values strictly inside the feasible set
  # (the same bounds used for the estimates: entries of Q of at least 1e-5,
  # with rows summing to one, and P in [1e-5, 1 - 1e-5]); this avoids
  # log(0) in the log-likelihood and keeps the Newton steps well defined
  eps <- 1e-5
  if (simplex) {
    # the projection works on the columns of a K x I matrix (the layout of
    # the C++ code), hence the transposes; rows that are already feasible
    # are left unchanged
    init <- t(radmix_project_q_cpp(t(init)))
  } else {
    init <- pmin(pmax(init, eps), 1 - eps)
  }
  # drop row and column names (e.g. the loci names of a previous P), which
  # the C++ code does not use
  unname(init)
}

# Check that an argument is a whole number of at least `min` (or, with
# scalar = FALSE, a vector of them), with an error naming the argument
check_count <- function(x, name, min = 1, scalar = TRUE, message = NULL) {
  ok <- is.numeric(x) && length(x) >= 1 && (!scalar || length(x) == 1) &&
    all(is.finite(x)) && all(x == round(x)) && all(x >= min)
  if (!ok) {
    if (is.null(message)) {
      message <- paste0("'", name, "' should be an integer of at least ", min)
    }
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

# Check that an argument is a single TRUE or FALSE
check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("'", name, "' should be TRUE or FALSE", call. = FALSE)
  }
  invisible(TRUE)
}

# The group of each individual of a grouped gen_tibble (NULL if it is not
# grouped), for the `group` element of a gt_admix object. As in
# gt_admixture(), a single grouping variable is supported
gt_single_group <- function(x) {
  if (!inherits(x, "grouped_gen_tbl")) {
    return(NULL)
  }
  group_vars <- dplyr::group_vars(x)
  if (length(group_vars) != 1) {
    stop(
      "'x' should be grouped by a single variable (it is grouped by ",
      paste(group_vars, collapse = ", "), ")",
      call. = FALSE
    )
  }
  x[[group_vars]]
}

#' Estimate ancestry proportions with the ADMIXTURE method
#'
#' This function estimates individual ancestry proportions (Q) and the allele
#' frequencies of the ancestral populations (P) with the method of ADMIXTURE
#' (Alexander et al. 2009), i.e. by maximum likelihood under the ADMIXTURE
#' model. It uses the optimisation strategy of ADMIXTURE: a few EM
#' iterations to prime a random start, then block-wise Newton updates of Q
#' and P, each solved as a constrained quadratic program, accelerated by the
#' quasi-Newton scheme of Zhou, Alexander and Lange (2011). It also offers
#' the cross-validation of ADMIXTURE to choose K (Alexander and Lange 2011),
#' a start from a truncated singular value decomposition (SVD) of the
#' genotypes, and support for pseudohaploid data. If you use this function,
#' please cite Alexander et al. (2009) for the method, and, if you use
#' cross-validation, Alexander and Lange (2011).
#'
#' Genotypes are read directly from the backing file of the `gen_tibble`
#' (respecting any filtering of individuals and selection of loci), so no
#' intermediate PLINK file is written. Computation is in C++, multithreaded
#' with OpenMP. The C++ code started as a port of the Julia package
#' [OpenADMIXTURE.jl](https://github.com/OpenMendel/OpenADMIXTURE.jl), and
#' has since been extensively rewritten (see the file `COPYRIGHTS` of the
#' package, `system.file("COPYRIGHTS", package = "tidygenclust")`).
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
#'   As in ADMIXTURE (Alexander et al. 2009; Alexander and Lange 2011), each
#'   run consists of:
#'   1. a random initialisation followed by a few EM iterations
#'      (`n_em_iter`; the EM algorithm of FRAPPE, Tang et al. 2005), or,
#'      with `init = "svd"`, a start from the SVD of the genotypes (see
#'      below);
#'   2. block-wise Newton updates (block relaxation): all of Q given P (one
#'      individual at a time, on the simplex), then all of P given the new
#'      Q (one locus at a time, in the unit interval). Each Newton step is
#'      solved exactly as a small quadratic program, by an active-set method
#'      based on the sweep operator (Lange 2010);
#'   3. quasi-Newton acceleration (Zhou, Alexander and Lange 2011): each
#'      iteration applies two rounds of block updates, and the changes that
#'      they made in the last `n_qn` iterations (secant pairs) give a
#'      quasi-Newton step, which is accepted only when it improves the
#'      log-likelihood of the previous iteration (otherwise the iteration
#'      ends with the plain block updates).
#'
#'   Iterations stop when the relative change in log-likelihood falls below
#'   `rtol`. All estimates are kept within \[1e-5, 1 - 1e-5\], as in
#'   ADMIXTURE. The log of each run (element `log` of the result) reports the
#'   log-likelihood at each iteration, together with that of the
#'   quasi-Newton point ("QN") and that of the plain block updates
#'   ("basic"). The latter is only computed when the quasi-Newton point is
#'   rejected (otherwise it is reported as "not computed"), which saves a
#'   pass over the genotypes in most iterations.
#'
#'   A few numerical safeguards handle degenerate cases, without changing
#'   the results when the computations are well behaved: for example, a
#'   Newton step that is not finite (because a genotype is impossible under
#'   the current estimates) is discarded, and the individual or locus keeps
#'   its current value for that iteration.
#'
#'   ## Starting values from the SVD
#'
#'   With `init = "svd"`, each run starts from the principal components of
#'   the genotypes rather than from random values, as in ADAMIXTURE
#'   (Saurina-i-Ricos et al. 2026) and fastmixture (Santander et al. 2024):
#'   1. a partial SVD of the genotypes, centred and scaled as by
#'      [tidypopgen::gt_pca_partialSVD()], with `max(k) + 5` components, is
#'      computed once and shared by all the runs (monomorphic loci, which
#'      the default scaling of tidypopgen rejects, are left unscaled);
#'   2. for each run, alternating least squares (ALS) find Q and P whose
#'      product approximates the allele frequencies predicted by the first
#'      `k + 5` components, starting from random allele frequencies. As in
#'      ADAMIXTURE, the least squares are bounded, the loci are weighted by
#'      the inverse of their binomial variance when solving for Q, and Q is
#'      extrapolated (momentum) when solving for P. The ALS only uses the
#'      SVD, not the genotypes, so it is fast;
#'   3. the best of 3 such solutions, by log-likelihood, starts the fit,
#'      without EM priming (`n_em_iter` is not used).
#'
#'   The SVD cannot use missing genotypes, so they must be imputed first
#'   (e.g. with [tidypopgen::gt_impute_simple()]); the imputed genotypes are
#'   only used for the SVD, and the fit ignores them (unless
#'   `use_imputed = TRUE`). The SVD also needs more individuals and loci than
#'   `max(k) + 5`.
#'
#'   Fits from the SVD need fewer iterations (often less than half as many as
#'   from random starts in our tests), and usually reach a good fit (though,
#'   as from any start, not always the best one). They are therefore useful
#'   for a single good fit of a large dataset, where many random runs would
#'   be costly, or as a first exploration before a more extensive search
#'   with random starts. However, the starts of different runs vary little,
#'   and favour some solutions over others: to explore the modes of the
#'   likelihood with several runs (e.g. with [gt_clumppling()]), use random
#'   starts (the default), and do not mix runs started in the two ways, as
#'   the runs started from the SVD would over-represent the solutions that
#'   they favour.
#'
#'   ## Starting values from a subset of the loci
#'
#'   With `init = "subset"`, each run starts with a fit to a random subset
#'   of the loci (a fraction `init_subset` of them, 5% by default), from a
#'   random start primed by `n_em_iter` EM iterations, to a relative
#'   tolerance of 1e-6. The ancestry proportions of that fit, and the
#'   allele frequencies that best fit each locus given them (one EM and two
#'   Newton steps over all the loci), then start the fit of all the loci,
#'   without EM priming. Passes over the subset cost a fraction of a pass
#'   over all the loci, and the fit of all the loci starts close to a mode,
#'   so it needs fewer iterations: in our tests on the HGDP (1,043
#'   individuals, K = 5), 9-11 iterations instead of 21-27, and 0.4-0.55
#'   times the passes over the genotypes of a random start, on 130,943 or
#'   642,794 loci. But the fit of the subset chooses the mode: on
#'   EUR_subset (379 individuals, 54,051 loci, K = 3, 10% of the loci), 4
#'   of 5 runs ended in a lower mode, against 1 of 5 from random starts,
#'   and on the HGDP at K = 10 (a single run, 5% of the loci) the run
#'   ended 67,600 log-likelihood units below the random start of the same
#'   seed (in 64 passes, against 206).
#'   This start is meant for a fast fit of a large dataset; to explore the
#'   modes of the likelihood with several runs (e.g. with [gt_clumppling()]),
#'   use random starts, and do not mix runs started in different ways.
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
#'     the same convergence criteria (`max_iter` and `rtol`) by default
#'     (see below for `cv_max_iter` and `cv_rtol`);
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
#'   The refits of the folds can be shortened with `cv_max_iter` (at most
#'   that many iterations per fold; ADAMIXTURE refits each fold with 5) or
#'   `cv_rtol` (a looser convergence threshold). This saves time, but biases
#'   the errors downwards: each refit starts from the estimates for the full
#'   data, which have seen the masked genotypes, and a refit stopped early
#'   still predicts them too well, the more so at larger `k`. In our tests
#'   (one seed per value of `k`), on the HGDP (1,043 individuals, 130,943
#'   loci, `k` = 2, 5 and 8) the refits converged in 6-7 iterations, and 5
#'   iterations per fold changed the errors by less than 1e-5. On EUR_subset
#'   (379 individuals, 54,051 loci, `k` = 1-5) they took 7-34 iterations: 5
#'   iterations saved 60% of the passes over the genotypes of
#'   cross-validation but lowered the error at `k` = 5 by 4e-4, and
#'   `cv_rtol = 1e-6` saved 50% and lowered it by 3e-4, while the errors of
#'   successive values of `k` differed by 1e-3 to 3e-3. Neither changed the
#'   value of `k` with the lowest error.
#'
#'   ## Reproducibility
#'
#'   The random starting values (and the folds of cross-validation) are
#'   drawn from R's random number generator, so runs are reproducible with
#'   `seed` (or [set.seed()]). The work is divided between threads in a way
#'   that depends only on the dimensions of the data, so the results are
#'   identical for any `n_cores`.
#'
#'   On x86-64 Linux, when the package is compiled with GCC 12 or later, the
#'   main computations are also compiled for processors with the AVX2 and
#'   FMA instructions (most processors made since 2015), and the faster
#'   version is chosen when the package is loaded. Results are reproducible
#'   on a given machine, but can differ in the last digits between machines
#'   with and without these instructions. Compiling the package with
#'   `-DRADMIX_NO_MULTIVERSION` (e.g. added to `CXXFLAGS` in
#'   `~/.R/Makevars`) builds only the portable version, whose results are the
#'   same on every x86-64 machine.
#'
#'   ## Comparison with ADMIXTURE
#'
#'   The model, the log-likelihood and the optimisation strategy are those of
#'   ADMIXTURE (Alexander et al. 2009; Alexander and Lange 2011), so both
#'   programs reach the same maximum likelihood estimates, up to local optima
#'   and the convergence criterion (ADMIXTURE stops when the log-likelihood
#'   increases by less than 1e-4 between iterations). The test suite checks
#'   this on the `EUR_subset` data of SnpArrays.jl. Note that:
#'   - ADMIXTURE reports the frequency of the second allele in the PLINK bim
#'     file, which for a `gen_tibble` is the reference allele, so its P
#'     matrices are 1 - P of those returned here;
#'   - cross-validation follows ADMIXTURE (see above), and the test suite
#'     checks that the errors agree with those of ADMIXTURE on `EUR_subset`;
#'   - ADMIXTURE treats the autosomal genotypes of every individual as
#'     diploid, whereas `gt_radmix()` counts a single allele per locus for
#'     pseudohaploid individuals (see above);
#'   - bootstrap standard errors and supervised analyses with individuals of
#'     known ancestry, available in ADMIXTURE, are not implemented (a fixed
#'     reference panel can be used with `init_p` and `fix_p`).
#'
#' @references Alexander DH, Novembre J, Lange K (2009) Fast model-based
#'   estimation of ancestry in unrelated individuals. Genome Research, 19(9),
#'   1655-1664.
#'
#'   Alexander DH, Lange K (2011) Enhancements to the ADMIXTURE algorithm for
#'   individual ancestry estimation. BMC Bioinformatics, 12, 246.
#'
#'   Zhou H, Alexander D, Lange K (2011) A quasi-Newton acceleration for
#'   high-dimensional optimization algorithms. Statistics and Computing,
#'   21(2), 261-273.
#'
#'   Tang H, Peng J, Wang P, Risch NJ (2005) Estimation of individual
#'   admixture: analytical and study design considerations. Genetic
#'   Epidemiology, 28(4), 289-301.
#'
#'   Lange K (2010) Numerical Analysis for Statisticians, 2nd edition.
#'   Springer, New York.
#'
#'   Saurina-i-Ricos J, Mas Montserrat D, Ioannidis AG (2026) ADAMIXTURE:
#'   adaptive first-order optimization for biobank-scale genetic clustering.
#'   Bioinformatics, 42(Supplement_1). doi:
#'   https://doi.org/10.1093/bioinformatics/btag236
#'
#'   Santander CG, Refoyo Martinez A, Meisner J (2024) Faster model-based
#'   estimation of ancestry proportions. Peer Community Journal, 4, e115.
#'   doi: https://doi.org/10.24072/pcjournal.503
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
#'   log-likelihood between iterations (defaults to 1e-8). A threshold of
#'   1e-7 takes 10-20% less time, but can stop a few log-likelihood units
#'   short of convergence, with ancestry proportions up to about 0.05 away
#'   from the converged ones
#' @param n_em_iter the number of EM iterations used to prime the estimates
#'   (defaults to 5, as in ADMIXTURE); not used with `init = "svd"`
#' @param n_qn the number of secant pairs used by the quasi-Newton
#'   acceleration (defaults to 3, as in ADMIXTURE)
#' @param init how the starting values of each run are chosen: `"random"`
#'   (the default), random values followed by `n_em_iter` EM iterations, as
#'   in ADMIXTURE; or `"svd"`, from a truncated SVD of the
#'   genotypes, which needs fewer iterations but is not suited to exploring
#'   the modes of the likelihood with several runs (see "Starting values
#'   from the SVD" in the details); or `"subset"`, a random start fitted to
#'   a random subset of the loci and extended to all of them, which needs
#'   fewer passes over the genotypes on large data (see "Starting values
#'   from a subset of the loci" in the details). `"svd"` and `"subset"`
#'   cannot be combined with `init_q` or `init_p`.
#' @param init_subset the fraction of the loci used by `init = "subset"`
#'   (defaults to 0.05)
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
#' @param cv_max_iter,cv_rtol the maximum number of iterations and the
#'   convergence threshold of the refit of each fold in cross-validation
#'   (default to `max_iter` and `rtol`, as for the full data). Smaller
#'   `cv_max_iter` or larger `cv_rtol` make cross-validation faster but bias
#'   its errors downwards (see "Cross-validation" in the details).
#' @param n_cores the number of cores to use (defaults to 1); the results
#'   do not depend on it
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
                      rtol = 1e-8, n_em_iter = 5, n_qn = 3,
                      init = c("random", "svd", "subset"),
                      init_subset = 0.05, init_q = NULL, init_p = NULL,
                      fix_q = FALSE, fix_p = FALSE, use_imputed = FALSE,
                      crossval = FALSE, n_folds = 5,
                      cv_max_iter = max_iter, cv_rtol = rtol, n_cores = 1,
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
  init <- match.arg(init)
  check_count(cv_max_iter, "cv_max_iter", min = 0)
  # rtol and cv_rtol must be positive numbers
  for (tol in c("rtol", "cv_rtol")) {
    value <- get(tol)
    tol_ok <- is.numeric(value) && length(value) == 1 && is.finite(value) &&
      value > 0
    if (!tol_ok) {
      stop("'", tol, "' should be a positive number", call. = FALSE)
    }
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
  if (init != "random" && (!is.null(init_q) || !is.null(init_p))) {
    stop("'init = \"", init, "\"' cannot be combined with 'init_q' or ",
      "'init_p'",
      call. = FALSE
    )
  }
  subset_ok <- is.numeric(init_subset) && length(init_subset) == 1 &&
    is.finite(init_subset) && init_subset > 0 && init_subset <= 1
  if (!subset_ok) {
    stop("'init_subset' should be a number in (0, 1]", call. = FALSE)
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
  # the SVD start uses K + 5 components, which needs more individuals and
  # loci than that (it is meant for large datasets anyway)
  svd_rank <- max(k) + 5
  if (init == "svd" && min(data$n_indiv, data$n_loci) <= svd_rank) {
    stop(
      "init = \"svd\" needs more than max(k) + 5 = ", svd_rank,
      " individuals and loci",
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

  # The SVD of the genotypes for the SVD starts (R/radmix_svd.R), with the
  # components needed by the largest k, computed once and shared by all the
  # runs (each k uses its first k + 5 components). It is computed after the
  # RNG state is saved above, in case it draws random numbers.
  if (init == "svd") {
    if (verbose) {
      message(
        "Computing the SVD of the genotypes (", svd_rank, " components)"
      )
    }
    svd_factors <- radmix_svd_factors(x, rank = svd_rank, n_cores = n_cores)
  }

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
      # priming is skipped, so all the runs give the same result. An SVD
      # start gives both, from ALS solutions that each run draws from its
      # own random allele frequencies; a subset start gives both, from the
      # fit of a random start to its own random subset of the loci
      # (R/radmix_subset.R).
      run_q <- init_q
      run_p <- init_p
      if (init == "subset") {
        start <- radmix_subset_start(data,
          k = this_k, frac = init_subset, max_iter = max_iter,
          n_em_iter = n_em_iter, n_qn = n_qn, n_cores = n_cores,
          verbose = verbose
        )
        run_q <- start$q
        run_p <- start$p
      }
      if (init == "svd") {
        start <- radmix_svd_start(data,
          factors = svd_factors, k = this_k, n_cores = n_cores,
          verbose = verbose
        )
        run_q <- start$q
        run_p <- start$p
      }
      fit <- radmix_fit(
        data,
        k = this_k, init_q = run_q, init_p = run_p,
        fix_q = fix_q, fix_p = fix_p, max_iter = max_iter, rtol = rtol,
        n_em_iter = n_em_iter, n_qn = n_qn, n_cores = n_cores,
        verbose = verbose
      )
      if (init != "random") fit$log <- c(start$log, fit$log)
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
          max_iter = cv_max_iter, rtol = cv_rtol, n_qn = n_qn,
          n_cores = n_cores
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
    init <- t(radmix_project_cpp(t(init), TRUE))
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

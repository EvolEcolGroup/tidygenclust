# Exploration of a start for gt_radmix() from a randomised SVD followed by
# alternating least squares (ALS) on its factors, as the default start of
# fastmixture (Santander et al. 2024; fastmixture/utils.py, randomSVD() and
# factorALS()), and ADAMIXTURE in a refined form, against the default random
# start of gt_radmix() (random values and 5 EM steps).
#
# The SVD is that of tidypopgen (gt_pca_randomSVD(), with the default binomial
# scaling of bigsnpr, i.e. the Patterson weights), on genotypes imputed with
# their (rounded) mean; the ALS follows factorALS() of fastmixture. Both
# starts are then fitted by gt_radmix() with its default settings (with both
# init_q and init_p given, gt_radmix() skips its EM steps).
#
# Run from the root of the package, with tidygenclust installed with
# R CMD INSTALL --preclean . (not loaded with load_all(), which compiles
# without optimisation), on an otherwise idle machine:
#
#   Rscript data-raw/radmix/explore_svd_start.R <output csv>

suppressPackageStartupMessages({
  library(tidypopgen)
  library(tidygenclust)
})

# A start for gt_radmix() at K = k, from a gen_tibble x without missing
# genotypes, or with imputed ones (imputed genotypes are only used by the
# SVD: gt_radmix() treats them as missing). Returns the starting Q
# (individuals x k) and P (loci x k, frequencies of the alternate allele, as
# counted by the FBM), the number of ALS iterations, and the time taken by
# the SVD and by the ALS.
svd_als_start <- function(x, k, seed = 1, max_iter = 1000, tol = 1e-4) {
  # the randomised SVD of the genotypes, centred by 2f and scaled by
  # sqrt(2f(1 - f)), with K - 1 components, as in fastmixture (no extra pass
  # for the total variance)
  t_svd <- system.time(
    pca <- gt_pca_randomSVD(x, k = k - 1, total_var = FALSE)
  )[["elapsed"]]
  t_als <- system.time({
    # the expected frequencies of the alternate allele, Pi (loci x
    # individuals) = 0.5 (G - 2f)' + f 1', with the centred G ~ U D V' S
    # (S the diagonal of the scales), i.e. Pi ~ 0.5 Z V' + f 1', with
    # Z = S V D (loci x (k - 1)) and V = U (individuals x (k - 1)), in the
    # notation of fastmixture; ALS looks for P Q' ~ Pi
    f <- pca$center / 2
    z <- sweep(pca$v, 2, pca$d, "*") * pca$scale
    v <- pca$u
    n <- nrow(v)
    m <- nrow(z)
    clip <- function(a) pmin(pmax(a, 1e-5), 1 - 1e-5)
    # Q given P, and P given Q (least squares, then mapped back to their
    # domains: P clipped, the rows of Q clipped and normalised)
    q_step <- function(p) {
      h <- p %*% MASS::ginv(crossprod(p))
      q <- 0.5 * v %*% crossprod(z, h) +
        matrix(colSums(h * f), n, k, byrow = TRUE)
      q <- clip(q)
      q / rowSums(q)
    }
    p_step <- function(q) {
      h <- q %*% MASS::ginv(crossprod(q))
      clip(0.5 * z %*% crossprod(v, h) + outer(f, colSums(h)))
    }
    set.seed(seed)
    p <- clip(matrix(stats::runif(m * k), m, k))
    q <- q_step(p)
    for (als_iter in seq_len(max_iter)) {
      q_old <- q
      p <- p_step(q)
      q <- q_step(p)
      if (sqrt(mean((q - q_old)^2)) < tol) {
        break
      }
    }
  })[["elapsed"]]
  stopifnot(
    all(abs(rowSums(q) - 1) < 1e-10),
    all(p >= 1e-5 & p <= 1 - 1e-5)
  )
  list(
    q = unname(q), p = unname(p), als_iter = als_iter,
    svd_seconds = t_svd, als_seconds = t_als
  )
}

# the number of iterations of a fit, from its log (NA if it did not
# converge)
n_iter <- function(fit) {
  line <- grep("Converged in", fit$log[[1]], value = TRUE)
  if (length(line) == 0) {
    return(NA_integer_)
  }
  as.integer(sub(".*Converged in ([0-9]+) iterations.*", "\\1", line[1]))
}

# the largest absolute difference between two Q matrices, after matching
# their clusters (trying all the permutations)
permutations <- function(n) {
  if (n == 1) {
    return(matrix(1L))
  }
  smaller <- permutations(n - 1)
  do.call(rbind, lapply(seq_len(n), function(first) {
    cbind(first, matrix(setdiff(seq_len(n), first)[smaller], ncol = n - 1))
  }))
}
max_q_diff <- function(q, q_ref) {
  perms <- permutations(ncol(q))
  min(apply(perms, 1, function(perm) max(abs(q[, perm] - q_ref))))
}

# one fit of gt_radmix() from a start (seed for a random start, or the
# matrices of an SVD start), with the log-likelihood of the start (after the
# EM steps, for a random start) and of the fit
run_fit <- function(x, k, seed = NULL, start = NULL) {
  args <- list(x, k = k)
  if (is.null(start)) {
    args$seed <- seed
  } else {
    args$init_q <- start$q
    args$init_p <- start$p
  }
  loglik_start <- do.call(gt_radmix, c(args, max_iter = 0))$loglik
  seconds <- system.time(fit <- do.call(gt_radmix, args))[["elapsed"]]
  list(
    loglik_start = loglik_start, loglik = fit$loglik, iterations = n_iter(fit),
    seconds = seconds, q = unclass(fit$Q[[1]])
  )
}

# the datasets, and the values of K
vcf_path <- system.file(
  "/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
  package = "tidypopgen"
)
datasets <- list(
  hgdp = list(
    file = "vignettes/articles/data/hgdp/hgdp_bench_quick.bed", k = c(3, 5, 7)
  ),
  eur = list(
    file = file.path(
      tools::R_user_dir("tidygenclust", which = "cache"), "EUR_subset.bed"
    ),
    k = 2:4
  ),
  anolis = list(file = vcf_path, k = 2:4)
)
random_seeds <- 1:5
svd_seeds <- 1:3

out_file <- commandArgs(trailingOnly = TRUE)[1]
results <- list()
for (dataset in names(datasets)) {
  file <- datasets[[dataset]]$file
  x <- if (grepl("\\.bed$", file)) {
    gen_tibble(file, backingfile = tempfile(dataset), quiet = TRUE)
  } else {
    gen_tibble(file,
      backingfile = tempfile(dataset), quiet = TRUE, parser = "cpp"
    )
  }
  # missing genotypes are imputed for the SVD (once per backing file)
  impute_seconds <- 0
  if (any(loci_missingness(x) > 0)) {
    impute_seconds <- system.time(
      x <- gt_impute_simple(x, method = "mean0")
    )[["elapsed"]]
  }
  message(
    dataset, ": ", nrow(x), " individuals, ", count_loci(x), " loci, ",
    "imputation ", impute_seconds, " s"
  )
  for (k in datasets[[dataset]]$k) {
    fits <- list()
    for (seed in random_seeds) {
      fit <- run_fit(x, k, seed = seed)
      fits[[length(fits) + 1]] <- c(
        list(start = "random", seed = seed, als_iter = NA, start_seconds = 0),
        fit
      )
      Sys.sleep(2)
    }
    for (seed in svd_seeds) {
      start <- svd_als_start(x, k, seed = seed)
      fit <- run_fit(x, k, start = start)
      fits[[length(fits) + 1]] <- c(
        list(
          start = "svd", seed = seed, als_iter = start$als_iter,
          start_seconds = start$svd_seconds + start$als_seconds
        ),
        fit
      )
      Sys.sleep(2)
    }
    # the differences from the best run of this K
    best <- which.max(vapply(fits, function(f) f$loglik, 0))
    for (fit in fits) {
      results[[length(results) + 1]] <- data.frame(
        dataset = dataset, k = k, start = fit$start, seed = fit$seed,
        als_iter = fit$als_iter, start_seconds = fit$start_seconds,
        impute_seconds = impute_seconds, fit_seconds = fit$seconds,
        iterations = fit$iterations, loglik_start = fit$loglik_start,
        loglik = fit$loglik, gap = fits[[best]]$loglik - fit$loglik,
        max_q_diff = max_q_diff(fit$q, fits[[best]]$q)
      )
    }
    message("  K = ", k, " done")
  }
}
results <- do.call(rbind, results)
utils::write.csv(results, out_file, row.names = FALSE)

# summary: per dataset, K and start, the median iterations and times, and
# how many runs reached the best log-likelihood (to within 0.1)
summary_tbl <- do.call(rbind, lapply(
  split(results, list(results$dataset, results$k, results$start), drop = TRUE),
  function(r) {
    data.frame(
      dataset = r$dataset[1], k = r$k[1], start = r$start[1], runs = nrow(r),
      iterations = stats::median(r$iterations),
      start_s = round(stats::median(r$start_seconds), 2),
      fit_s = round(stats::median(r$fit_seconds), 2),
      at_best = sum(r$gap < 0.1),
      max_gap = round(max(r$gap), 1),
      max_q_diff = signif(max(r$max_q_diff), 2)
    )
  }
))
summary_tbl <- summary_tbl[
  order(summary_tbl$dataset, summary_tbl$k, summary_tbl$start),
]
rownames(summary_tbl) <- NULL
print(summary_tbl)

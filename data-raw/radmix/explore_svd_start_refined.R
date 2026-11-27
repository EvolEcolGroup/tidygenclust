# Exploration of refined SVD+ALS starts for gt_radmix() (following
# explore_svd_start.R, in the same directory): a cheaper SVD, a better start,
# and enough variability between runs to explore several modes.
#
# - The SVD is computed once per dataset, with as many components as the
#   largest K (as ADAMIXTURE does), with tidypopgen::gt_pca_partialSVD()
#   (Patterson scaling), which gives the same subspace as
#   gt_pca_randomSVD() in a fraction of the time; missing genotypes are
#   imputed with their (rounded) mean for the SVD only.
# - Starts from it, by ALS on its factors (no pass over the genotypes):
#   "simple", as fastmixture (factorALS(): K - 1 components, unconstrained
#   least squares, then clipping); "refined", as ADAMIXTURE (model/als.py:
#   K components, least squares bounded in [0, 1], Fisher-weighted Q step,
#   momentum 0.8 on Q).
# - Variability: a refined start perturbed towards a random start (Q mixed
#   with Dirichlet(1) rows, P with uniform values), by a fraction `perturb`;
#   or a refined start whose ALS uses a random subset of the loci (a
#   fraction `subsample`), P then being estimated for all the loci from its
#   Q (one P step).
# - Selection: the best of `best_of` refined starts (from different ALS
#   seeds), by their starting log-likelihood (one pass over the genotypes
#   each); run r uses ALS seeds (r - 1) best_of + 1 to r best_of, so that
#   runs still start from different ALS solutions.
# Each start is then fitted by gt_radmix() with its default settings, and
# compared with random starts (random values and 5 EM steps), for how many
# iterations and how much time they take, how often they reach the best
# log-likelihood, and how many distinct modes repeated runs find.
#
# Run from the root of the package, with tidygenclust installed with
# R CMD INSTALL --preclean . (not loaded with load_all()):
#
#   Rscript data-raw/radmix/explore_svd_start_refined.R <output csv> [seeds] \
#     [strategies] [datasets]

suppressPackageStartupMessages({
  library(tidypopgen)
  library(tidygenclust)
})

# the SVD of a gen_tibble (imputed for missing genotypes) with `rank`
# components, as the factors of the expected frequencies of the alternate
# allele: Pi (loci x individuals) = 0.5 (G - 2f)' + f 1', with the centred
# G ~ U D V' S (S the diagonal of the Patterson scales), i.e.
# Pi ~ 0.5 Z U' + f 1', Z = S V D (loci x rank), U (individuals x rank)
svd_factors <- function(x, rank) {
  seconds <- system.time(
    pca <- gt_pca_partialSVD(x, k = rank, total_var = FALSE)
  )[["elapsed"]]
  list(
    z = sweep(pca$v, 2, pca$d, "*") * pca$scale, u = pca$u,
    f = pca$center / 2, seconds = seconds
  )
}

clip <- function(a) pmin(pmax(a, 1e-5), 1 - 1e-5)
map_q <- function(q) {
  q <- clip(q)
  q / rowSums(q)
}

# the rows of x that minimise 0.5 x' A x - b' x within [0, 1], by cyclic
# coordinate descent over all the rows at once, starting from x
box_ls <- function(a, b, x, sweeps = 20, tol = 1e-7) {
  for (s in seq_len(sweeps)) {
    change <- 0
    for (k in seq_len(ncol(b))) {
      new <- (b[, k] - x %*% a[, k] + x[, k] * a[k, k]) / a[k, k]
      new <- pmin(pmax(new, 0), 1)
      change <- max(change, abs(new - x[, k]))
      x[, k] <- new
    }
    if (change < tol) break
  }
  x
}

# a start at K = k from the SVD factors; method "simple" (fastmixture) or
# "refined" (ADAMIXTURE); perturb: the weight of a random start mixed in
# extra: the number of components beyond K used by the refined ALS (a
# number, or "k" for 2K components in all)
als_start <- function(fac, k, seed, method = c("simple", "refined"),
                      perturb = 0, subsample = 1, extra = 0, max_iter = 1000,
                      tol = 1e-4) {
  method <- match.arg(method)
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  if (identical(extra, "k")) {
    extra <- k
  }
  rank <- if (method == "simple") k - 1 else k + extra
  # Pi = za %*% t(ua), loci x individuals
  za_all <- cbind(0.5 * fac$z[, seq_len(rank), drop = FALSE], fac$f)
  f_all <- fac$f
  # the loci used by the ALS (all, or a random subset)
  loci <- if (subsample < 1) {
    sort(sample(nrow(za_all), round(subsample * nrow(za_all))))
  } else {
    seq_len(nrow(za_all))
  }
  za <- za_all[loci, , drop = FALSE]
  f <- f_all[loci]
  ua <- cbind(fac$u[, seq_len(rank), drop = FALSE], 1)
  m <- nrow(za)
  n <- nrow(ua)
  p <- clip(matrix(stats::runif(m * k), m, k))
  if (method == "simple") {
    q_step <- function(p, q) {
      h <- p %*% MASS::ginv(crossprod(p))
      map_q(ua %*% crossprod(za, h))
    }
    p_step <- function(q, p) {
      h <- q %*% MASS::ginv(crossprod(q))
      clip(za %*% crossprod(ua, h))
    }
  } else {
    # Fisher weights of the loci, for the Q step
    w <- 1 / pmax(f * (1 - f), 1e-3)
    w <- w / mean(w)
    q_step <- function(p, q) {
      wp <- p * w
      map_q(box_ls(crossprod(p, wp), ua %*% crossprod(za, wp), q))
    }
    p_step <- function(q, p) {
      clip(box_ls(crossprod(q), za %*% crossprod(ua, q), p))
    }
  }
  q <- q_step(p, matrix(1 / k, n, k))
  q_prev <- q
  for (iter in seq_len(max_iter)) {
    q_use <- if (method == "refined" && iter > 1) {
      map_q(q + 0.8 * (q - q_prev))
    } else {
      q
    }
    p <- p_step(q_use, p)
    q_new <- q_step(p, q)
    done <- sqrt(mean((q_new - q)^2)) < tol
    q_prev <- q
    q <- q_new
    if (done) break
  }
  # with a subset of the loci, P of all the loci, from Q (one P step, from
  # the frequencies of the alleles)
  if (subsample < 1) {
    za <- za_all
    m <- nrow(za)
    p <- clip(box_ls(crossprod(q), za %*% crossprod(ua, q),
      matrix(f_all, m, k)
    ))
  }
  if (perturb > 0) {
    d <- matrix(stats::rexp(n * k), n, k)
    q <- map_q((1 - perturb) * q + perturb * d / rowSums(d))
    p <- clip((1 - perturb) * p + perturb * stats::runif(m * k))
  }
  list(
    q = unname(q), p = unname(p), als_iter = iter,
    seconds = proc.time()[["elapsed"]] - t0
  )
}

# the number of iterations of a fit, from its log
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

# the number of distinct modes among Q matrices (greedily: a Q matrix within
# 0.05 of a mode found earlier belongs to it)
n_modes <- function(qs) {
  modes <- list()
  for (q in qs) {
    if (!any(vapply(modes, function(m) max_q_diff(q, m) < 0.05, NA))) {
      modes[[length(modes) + 1]] <- q
    }
  }
  length(modes)
}

# one fit from a start (seed for a random start, or the matrices of an SVD
# start)
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

vcf_path <- system.file(
  "/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
  package = "tidypopgen"
)
datasets <- list(
  hgdp = list(
    file = "vignettes/articles/data/hgdp/hgdp_bench_quick.bed", k = c(5, 7)
  ),
  eur = list(
    file = file.path(
      tools::R_user_dir("tidygenclust", which = "cache"), "EUR_subset.bed"
    ),
    k = 3:4
  ),
  anolis = list(file = vcf_path, k = 3:4)
)
all_strategies <- list(
  random = list(),
  simple = list(method = "simple"),
  refined = list(method = "refined"),
  refined_p30 = list(method = "refined", perturb = 0.3),
  refined_p60 = list(method = "refined", perturb = 0.6),
  refined_sub50 = list(method = "refined", subsample = 0.5),
  refined_sub25 = list(method = "refined", subsample = 0.25),
  refined_x2 = list(method = "refined", extra = 2),
  refined_x5 = list(method = "refined", extra = 5),
  refined_2k = list(method = "refined", extra = "k"),
  refined_x5_best3 = list(method = "refined", extra = 5, best_of = 3),
  refined_x5_best5 = list(method = "refined", extra = 5, best_of = 5)
)

# the output file, the number of seeds, and the strategies (all if not
# given), e.g. random,refined,refined_sub50
args <- commandArgs(trailingOnly = TRUE)
out_file <- args[1]
seeds <- if (length(args) >= 2) seq_len(as.integer(args[2])) else 1:5
strategies <- if (length(args) >= 3) {
  all_strategies[strsplit(args[3], ",")[[1]]]
} else {
  all_strategies
}
# and the datasets (all if not given), e.g. hgdp,eur
if (length(args) >= 4) {
  datasets <- datasets[strsplit(args[4], ",")[[1]]]
}
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
  impute_seconds <- 0
  if (any(loci_missingness(x) > 0)) {
    impute_seconds <- system.time(
      x <- gt_impute_simple(x, method = "mean0")
    )[["elapsed"]]
  }
  # one SVD for all the values of K and all the runs (with enough
  # components for the largest K, and the extra ones: up to K + 5 or 2K)
  k_max <- max(datasets[[dataset]]$k)
  fac <- svd_factors(x, max(2 * k_max, k_max + 5))
  message(
    dataset, ": imputation ", impute_seconds, " s, SVD ", fac$seconds, " s"
  )
  for (k in datasets[[dataset]]$k) {
    fits <- list()
    for (strategy in names(strategies)) {
      st <- strategies[[strategy]]
      for (seed in seeds) {
        best_of <- if (is.null(st$best_of)) 1 else st$best_of
        start <- NULL
        if (strategy != "random") {
          # the best of best_of starts, by their starting log-likelihood
          t0 <- proc.time()[["elapsed"]]
          for (i in seq_len(best_of)) {
            candidate <- als_start(fac, k, (seed - 1) * best_of + i,
              method = st$method,
              perturb = if (is.null(st$perturb)) 0 else st$perturb,
              subsample = if (is.null(st$subsample)) 1 else st$subsample,
              extra = if (is.null(st$extra)) 0 else st$extra
            )
            if (best_of > 1) {
              candidate$loglik <- gt_radmix(x,
                k = k, init_q = candidate$q, init_p = candidate$p,
                max_iter = 0
              )$loglik
            }
            if (is.null(start) || candidate$loglik > start$loglik) {
              start <- candidate
            }
          }
          start$seconds <- proc.time()[["elapsed"]] - t0
        }
        fit <- run_fit(x, k, seed = seed, start = start)
        fits[[length(fits) + 1]] <- c(
          list(
            strategy = strategy, seed = seed,
            als_iter = if (is.null(start)) NA else start$als_iter,
            als_seconds = if (is.null(start)) 0 else start$seconds
          ),
          fit
        )
        Sys.sleep(1)
      }
    }
    best <- max(vapply(fits, function(f) f$loglik, 0))
    best_q <- fits[[which.max(vapply(fits, function(f) f$loglik, 0))]]$q
    for (strategy in names(strategies)) {
      sf <- Filter(function(f) f$strategy == strategy, fits)
      modes <- n_modes(lapply(sf, function(f) f$q))
      for (f in sf) {
        results[[length(results) + 1]] <- data.frame(
          dataset = dataset, k = k, strategy = strategy, seed = f$seed,
          svd_seconds = fac$seconds, als_iter = f$als_iter,
          als_seconds = f$als_seconds, fit_seconds = f$seconds,
          iterations = f$iterations, loglik_start = f$loglik_start,
          loglik = f$loglik, gap = best - f$loglik,
          q_diff_best = max_q_diff(f$q, best_q), modes = modes
        )
      }
    }
    message("  K = ", k, " done")
    # saved after each K, so that an interruption keeps what was done
    utils::write.csv(do.call(rbind, results), out_file, row.names = FALSE)
  }
}
results <- do.call(rbind, results)
utils::write.csv(results, out_file, row.names = FALSE)

# summary per dataset, K and strategy
summary_tbl <- do.call(rbind, lapply(
  split(results, list(results$dataset, results$k, results$strategy),
    drop = TRUE
  ),
  function(r) {
    data.frame(
      dataset = r$dataset[1], k = r$k[1], strategy = r$strategy[1],
      runs = nrow(r), iterations = stats::median(r$iterations),
      als_s = round(stats::median(r$als_seconds), 2),
      fit_s = round(stats::median(r$fit_seconds), 1),
      at_best = sum(r$gap < 0.1), near_best = sum(r$gap < 5),
      modes = r$modes[1],
      max_gap = round(max(r$gap), 1)
    )
  }
))
summary_tbl <- summary_tbl[
  order(summary_tbl$dataset, summary_tbl$k,
    match(summary_tbl$strategy, names(strategies))),
]
rownames(summary_tbl) <- NULL
print(summary_tbl)

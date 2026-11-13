# Part B: starts and warm-ups of the main fit, compared in passes over the
# genotypes (full-pass equivalents), final log-likelihood and iterations.
#
# Methods (`methods`, comma-separated):
#   random        the default: random start, 5 EM steps, fit (radmix_fit())
#   svd           init = "svd": SVD factors, best of 3 ALS starts, fit
#                 (the SVD is not counted in passes: see its time)
#   sub<pct>      warm-up on <pct>% of the loci (random start, 5 EM, fit to
#                 rtol 1e-6 on the subset), P extended to all the loci given
#                 Q (1 EM + 2 Newton P passes), then the full fit
#   mb            the random start and 5 EM steps of `random`, admixer's
#                 mini-batch warm-up (B = min(32, M / 500), halved when an
#                 epoch does not raise the log-likelihood), then the full fit
#   mbtol<x>      as mb, also halving B when an epoch gains less than 10^-x
#                 (relative)
# The final fit is always radmix_fit() at rtol 1e-8 from the warm-up's
# estimates. The same seed gives the same random start to `random` and `mb`.
#
# Usage:
#   Rscript warmup_compare.R <lib> <dataset> <k> <seeds> <methods> <n_cores> <out_dir>
# e.g. Rscript warmup_compare.R $LIB hgdp 5 1:5 random,sub10,sub25,mb 1 $OUT/warmup
# One .rds per (dataset, K, seed, method); existing files are skipped.
# Summarise with warmup_summary.R.
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
dataset <- args[2]
k <- as.integer(args[3])
seeds <- eval(parse(text = args[4]))
methods <- strsplit(args[5], ",")[[1]]
n_cores <- as.integer(args[6])
out_dir <- args[7]
library(tidygenclust, lib.loc = lib)
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(),
  value = TRUE
)))
source(file.path(script_dir, "algo_lib.R"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
x <- load_gt(dataset)
data <- tg$radmix_prepare(x)
svd_factors <- NULL

for (seed in seeds) {
  for (method in methods) {
    out_file <- file.path(out_dir, sprintf(
      "%s_k%02d_s%d_%s.rds", dataset, k, seed, method
    ))
    if (file.exists(out_file)) next
    set.seed(seed)
    res <- list(dataset = dataset, k = k, seed = seed, method = method)
    t0 <- proc.time()[["elapsed"]]
    if (method == "random") {
      fit <- tg$radmix_fit(data, k = k, n_cores = n_cores)
      res$warm_passes <- 0
      res$fit_passes <- fit_passes(fit)
    } else if (method == "svd") {
      if (is.null(svd_factors)) {
        # (the SVD needs the missing genotypes imputed: a separate imputed
        # copy of the data, with the same individuals and loci; the fit
        # itself uses the plain data, ignoring the imputed genotypes, as
        # gt_radmix() does)
        x_svd <- x
        if (any(tidypopgen::loci_missingness(x, n_cores = n_cores) > 0)) {
          x_svd <- load_gt(paste0(dataset, "_imputed"))
        }
        t_svd <- system.time(svd_factors <- tg$radmix_svd_factors(x_svd,
          rank = k + 5, n_cores = n_cores
        ))[["elapsed"]]
        res$time_svd <- t_svd
      }
      # (factors of rank k + 5 serve this K only)
      start <- tg$radmix_svd_start(data,
        factors = svd_factors, k = k, n_cores = n_cores
      )
      res$warm_passes <- 3
      t1 <- proc.time()[["elapsed"]]
      fit <- tg$radmix_fit(data, k = k, init_q = start$q, init_p = start$p,
        n_cores = n_cores)
      res$time_fit <- proc.time()[["elapsed"]] - t1
      res$fit_passes <- fit_passes(fit, em_ran = FALSE)
    } else if (grepl("^sub[0-9]+$", method)) {
      frac <- as.numeric(sub("^sub", "", method)) / 100
      w <- warmup_subset(data, k, frac, n_cores = n_cores)
      res$warm_passes <- w$passes
      res$sub_iter <- w$sub_iter
      t1 <- proc.time()[["elapsed"]]
      res$time_warm <- t1 - t0
      fit <- tg$radmix_fit(data, k = k, init_q = w$q, init_p = w$p,
        n_cores = n_cores)
      res$fit_passes <- fit_passes(fit, em_ran = FALSE)
      res$warm_loglik <- fit$loglik_history[1]
    } else if (grepl("^mb", method)) {
      rtol_halve <- if (grepl("^mbtol", method)) {
        10^-as.numeric(sub("^mbtol", "", method))
      } else {
        0
      }
      st <- random_start(k, data$n_indiv, data$n_loci)
      for (i in 1:5) st <- tg$radmix_em_cpp(data, st$q, st$p, n_cores)
      w <- warmup_minibatch(data, st$q, st$p, rtol_halve = rtol_halve,
        n_cores = n_cores)
      res$warm_passes <- 5 + w$passes
      res$epochs <- w$epochs
      res$warm_history <- w$history
      t1 <- proc.time()[["elapsed"]]
      res$time_warm <- t1 - t0
      fit <- tg$radmix_fit(data, k = k, init_q = w$q, init_p = w$p,
        n_cores = n_cores)
      res$fit_passes <- fit_passes(fit, em_ran = FALSE)
      res$warm_loglik <- fit$loglik_history[1]
    } else {
      stop("unknown method ", method)
    }
    res$time <- proc.time()[["elapsed"]] - t0
    res$loglik <- fit$loglik
    res$n_iter <- fit$n_iter
    res$converged <- fit$converged
    res$passes <- res$warm_passes + res$fit_passes
    res$q <- fit$q
    res$load <- readLines("/proc/loadavg")
    saveRDS(res, out_file)
    cat(sprintf(
      "%s K=%d seed=%d %-8s passes %6.1f (warm-up %5.1f, fit %3d it) loglik %.3f time %.1f s | load %s\n",
      dataset, k, seed, method, res$passes, res$warm_passes, res$n_iter,
      res$loglik, res$time, res$load
    ))
  }
}

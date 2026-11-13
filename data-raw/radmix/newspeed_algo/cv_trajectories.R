# Part A: the trajectory of each cross-validation refit, scored on the
# held-out genotypes after every iteration, so that every refit policy
# (converged, a fixed number of iterations, a looser tolerance) can be read
# from one run per fold (the policies all stop somewhere along the same
# deterministic trajectory: a refit with max_iter = m, or with a looser
# rtol, gives exactly the iterate at which this trajectory would stop).
#
# For each K and seed: the main fit as gt_radmix(seed = seed) runs it (random
# start, 5 EM steps, rtol 1e-8), the fold seed drawn as gt_radmix() draws it
# (so the converged policy reproduces gt_radmix(crossval = TRUE) exactly),
# then each of the 5 folds refitted from the full-data estimates to rtol
# 1e-8 by the R replica of the loop (algo_lib.R, bit-identical to
# radmix_fit(), see check_replica.R), with the held-out deviance of the
# iterate after every iteration.
#
# Usage:
#   Rscript cv_trajectories.R <lib> <dataset> <ks> <seeds> <n_cores> <out_dir>
# e.g. Rscript cv_trajectories.R $LIB hgdp 2:8 1 1 $OUT/cv
# One .rds per (dataset, K, seed) in out_dir; existing files are skipped,
# so an interrupted job can be restarted. Summarise with cv_policies.R.
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
dataset <- args[2]
ks <- eval(parse(text = args[3]))
seeds <- eval(parse(text = args[4]))
n_cores <- as.integer(args[5])
out_dir <- args[6]
n_folds <- 5
library(tidygenclust, lib.loc = lib)
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(),
  value = TRUE
)))
source(file.path(script_dir, "algo_lib.R"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
x <- load_gt(dataset)
data <- tg$radmix_prepare(x)

for (seed in seeds) {
  for (k in ks) {
    out_file <- file.path(out_dir, sprintf("%s_k%02d_s%d.rds", dataset, k, seed))
    if (file.exists(out_file)) next
    # the main fit, exactly as gt_radmix(seed = seed) runs it
    set.seed(seed)
    t_fit <- system.time(
      fit <- tg$radmix_fit(data, k = k, n_cores = n_cores)
    )[["elapsed"]]
    # the fold seed, drawn after the main fit as in gt_radmix()
    cv_seed <- floor(stats::runif(1) * 2^31)
    traj <- list()
    t_cv <- system.time(for (fold in seq_len(n_folds)) {
      cv <- list(fold = fold, n_folds = n_folds, seed = cv_seed)
      train <- tg$radmix_set_cv(data, c(list(mode = "train"), cv))
      held <- tg$radmix_set_cv(data, c(list(mode = "held_out"), cv))
      rows <- list()
      score <- function(iter, q, p, ll, reldiff, passes) {
        dev <- tg$radmix_deviance_cpp(held, q, p, n_cores)
        rows[[length(rows) + 1]] <<- data.frame(
          fold = fold, iter = iter, loglik = ll, reldiff = reldiff,
          dev = dev[1], n = dev[2], passes = passes
        )
      }
      r <- fit_replica(train, t(fit$q), t(fit$p),
        rtol = 1e-8, max_iter = 1000, n_cores = n_cores, on_iter = score
      )
      traj[[fold]] <- do.call(rbind, rows)
      cat(dataset, "K", k, "seed", seed, "fold", fold, ":", r$n_iter,
        "iterations\n")
    })[["elapsed"]]
    res <- list(
      dataset = dataset, k = k, seed = seed, cv_seed = cv_seed,
      fit = list(
        loglik = fit$loglik, n_iter = fit$n_iter, converged = fit$converged,
        passes = fit_passes(fit), time = t_fit
      ),
      traj = do.call(rbind, traj), time_cv = t_cv,
      load = readLines("/proc/loadavg")
    )
    saveRDS(res, out_file)
    # the estimates, for later reuse (e.g. as reference modes)
    saveRDS(list(q = fit$q, p = fit$p, loglik = fit$loglik),
      sub("\\.rds$", "_fit.rds", out_file)
    )
    cat(dataset, "K", k, "seed", seed, "main fit", fit$n_iter, "iterations",
      round(t_fit, 1), "s, CV", round(t_cv, 1), "s, load", res$load, "\n")
  }
}

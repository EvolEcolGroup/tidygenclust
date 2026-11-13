# Wall-clock times of the cross-validation refit policies, run through the
# package's own radmix_cv() (no replica): for each K, the main fit as in
# gt_radmix(seed = seed), then radmix_cv() with each policy, timed, with the
# CV error and the iterations of each fold. For the desktop (n_cores = 8,
# threads bound to the P-cores; see README.md).
#
# Policies (`policies`, comma-separated): conv (max_iter 1000, rtol 1e-8,
# the current one), it<m> (max_iter = m), rtol<r> (rtol = r).
#
# Usage:
#   Rscript cv_policy_timing.R <lib> <dataset> <ks> <seed> <policies> <n_cores> <out_csv>
# e.g. Rscript cv_policy_timing.R $LIB hgdp 2:8 1 conv,it5,it10,rtol1e-6,rtol1e-5 8 cv_timing_hgdp.csv
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
dataset <- args[2]
ks <- eval(parse(text = args[3]))
seed <- as.integer(args[4])
policies <- strsplit(args[5], ",")[[1]]
n_cores <- as.integer(args[6])
out_csv <- args[7]
library(tidygenclust, lib.loc = lib)
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(),
  value = TRUE
)))
source(file.path(script_dir, "algo_lib.R"))
x <- load_gt(dataset)
data <- tg$radmix_prepare(x)
rows <- list()
for (k in ks) {
  set.seed(seed)
  t_fit <- system.time(fit <- tg$radmix_fit(data, k = k, n_cores = n_cores))
  cv_seed <- floor(stats::runif(1) * 2^31)
  for (policy in policies) {
    max_iter <- 1000
    rtol <- 1e-8
    if (grepl("^it", policy)) max_iter <- as.integer(sub("^it", "", policy))
    if (grepl("^rtol", policy)) rtol <- as.numeric(sub("^rtol", "", policy))
    t_cv <- system.time(cv <- tg$radmix_cv(data,
      fit = fit, k = k, n_folds = 5, cv_seed = cv_seed,
      max_iter = max_iter, rtol = rtol, n_cores = n_cores
    ))[["elapsed"]]
    fold_iter <- as.integer(sub(".* ([0-9]+) iterations.*", "\\1",
      grep("^CV fold", cv$log, value = TRUE)))
    rows[[length(rows) + 1]] <- data.frame(
      dataset = dataset, k = k, seed = seed, policy = policy,
      cv = cv$cv, fold_iter = paste(fold_iter, collapse = "/"),
      time_cv = t_cv, time_fit = t_fit[["elapsed"]], fit_iter = fit$n_iter,
      n_cores = n_cores, load1 = as.numeric(strsplit(readLines(
        "/proc/loadavg"), " ")[[1]][1])
    )
    print(utils::tail(rows, 1)[[1]], row.names = FALSE)
    utils::write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
  }
}

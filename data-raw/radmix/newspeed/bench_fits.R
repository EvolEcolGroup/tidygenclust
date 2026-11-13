# Fits of gt_radmix() with one build of tidygenclust: random start (with the
# default 5 EM steps), one seed per fit, as gt_radmix(x, k, seed = seed)
# does (the fit is run through the internal radmix_fit(), so that the full
# history of log-likelihoods is kept; the time excludes radmix_prepare(),
# which costs milliseconds).
#
#   Rscript data-raw/radmix/newspeed/bench_fits.R <lib> <label> <out dir> \
#     <cases> [seeds] [n_cores] [max_iter] [rtol]
#
# cases: dataset:K pairs, e.g. "hgdp:5,eur:4" (datasets: see common.R);
# seeds: e.g. "1,2,3" (default 1); max_iter: default 1000 (a small value,
# e.g. 8, times a fixed number of iterations, the same work for every build);
# rtol: default 1e-8 (the default of gt_radmix()).
# Appends one row per fit to <out dir>/fits.csv (seconds, iterations,
# passes, log-likelihood, 1-minute load average before and after the fit),
# and saves the log-likelihood history and Q of each fit in
# <out dir>/fits/<label>_<data>_k<K>_s<seed>_c<n_cores>_i<max_iter>.rds,
# which summarise.R compares between builds.
args <- commandArgs(TRUE)
if (length(args) < 4) stop("usage: bench_fits.R <lib> <label> <out> <cases> ...")
lib <- args[1]
label <- args[2]
out <- args[3]
cases <- args[4]
seeds <- if (length(args) >= 5) args[5] else "1"
n_cores <- if (length(args) >= 6) as.integer(args[6]) else 1L
max_iter <- if (length(args) >= 7) as.integer(args[7]) else 1000L
rtol <- if (length(args) >= 8) as.numeric(args[8]) else 1e-8
pause <- as.numeric(Sys.getenv("PAUSE", "2"))
source("data-raw/radmix/newspeed/common.R")
ns <- tgc_load(lib)
dir.create(file.path(out, "fits"), showWarnings = FALSE, recursive = TRUE)

for (cs in tgc_cases(cases)) {
  x <- tgc_read(cs$data)
  dat <- ns$radmix_prepare(x)
  for (seed in tgc_ints(seeds)) {
    Sys.sleep(pause)
    load_before <- tgc_load_avg()
    t <- system.time({
      set.seed(seed)
      fit <- ns$radmix_fit(dat,
        k = cs$k, max_iter = max_iter, rtol = rtol, n_cores = n_cores
      )
    })[["elapsed"]]
    counts <- tgc_log_counts(fit$log)
    # written after each fit, so that an interrupted run keeps its results
    tgc_append(data.frame(
      label = label, data = cs$data, k = cs$k, seed = seed,
      n_cores = n_cores, max_iter = max_iter, rtol = rtol, seconds = t,
      iterations = fit$n_iter, converged = fit$converged,
      rejected = counts[["rejected"]],
      passes_before = counts[["passes_before"]],
      passes_fused = counts[["passes_fused"]],
      seconds_per_iter = t / fit$n_iter,
      loglik = sprintf("%.10f", fit$loglik),
      load_before = load_before, load_after = tgc_load_avg(),
      time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ), file.path(out, "fits.csv"))
    saveRDS(
      list(
        loglik_history = fit$loglik_history, q = fit$q, log = fit$log
      ),
      file.path(out, "fits", sprintf(
        "%s_%s_k%d_s%d_c%d_i%d.rds", label, cs$data, cs$k, seed, n_cores,
        max_iter
      ))
    )
  }
}

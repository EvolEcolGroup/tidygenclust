# gt_radmix() fits of a build of tidygenclust: seed 1, random start.
#
#   Rscript bench_fits.R <lib> <label> <out dir> [n_cores] [cv|nocv] [rtol] [cases]
#
# cases: dataset:K pairs; default "eur:3,eur:4,hgdp_quick:5,hgdp_quick:7"
# without cross-validation, "hgdp_quick:5" with it (those of the laptop).
# Appends to <out>/fits.csv (seconds, iterations, log-likelihood, CV error)
# and saves Q, P and the log-likelihood of each fit in
# <out>/fits_<label>_<n_cores>[_cv]_<rtol>.rds (summarise.R compares them
# between builds).
args <- commandArgs(TRUE)
lib <- args[1]
label <- args[2]
out <- args[3]
n_cores <- if (length(args) >= 4) as.integer(args[4]) else 1L
crossval <- length(args) >= 5 && args[5] == "cv"
rtol <- if (length(args) >= 6) as.numeric(args[6]) else 1e-7
cases <- if (length(args) >= 7 && nzchar(args[7])) {
  args[7]
} else if (crossval) {
  "hgdp_quick:5"
} else {
  "eur:3,eur:4,hgdp_quick:5,hgdp_quick:7"
}
pause <- as.numeric(Sys.getenv("PAUSE", "3"))
source("data-raw/radmix/perf_verification/scripts/common.R")
ns <- tgc_load(lib)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

cases <- tgc_cases(cases)
fits <- list()
rds <- file.path(out, paste0(
  "fits_", label, "_", n_cores, if (crossval) "_cv", "_", rtol, ".rds"
))
for (d in unique(vapply(cases, `[[`, "", "data"))) {
  x <- tgc_read(d)
  for (cs in Filter(function(z) z$data == d, cases)) {
    Sys.sleep(pause)
    t <- system.time(a <- gt_radmix(x,
      k = cs$k, seed = 1, n_cores = n_cores,
      crossval = crossval, rtol = rtol
    ))[["elapsed"]]
    fits[[paste(d, cs$k)]] <- list(
      q = unclass(a$Q[[1]]), p = a$P[[1]],
      loglik = a$loglik, cv = a$cv
    )
    # written after each fit, so that an interrupted run keeps its results
    tgc_append(data.frame(
      label = label, data = d, k = cs$k, n_cores = n_cores,
      crossval = crossval, rtol = rtol, seconds = t,
      iterations = tgc_iterations(a), loglik = a$loglik,
      cv = if (crossval) a$cv else NA
    ), file.path(out, "fits.csv"))
    saveRDS(fits, rds)
  }
}

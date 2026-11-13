# Times single passes over the genotypes (the kernels of gt_radmix()) with
# one build of tidygenclust, from the same (q, p) for every build: a random
# start followed by 5 EM steps (seed 1). The passes: loglik, update_q,
# update_p, and update_q_loglik (the fused pass, in builds that have it); on
# the plain view of the genotypes ("none") and, optionally, on the training
# view of cross-validation (fold 1 of 5 hidden, "train").
#
#   Rscript data-raw/radmix/newspeed/bench_kernels.R <lib> <label> <out dir> \
#     <cases> [n_cores] [reps] [views]
#
# cases: dataset:K pairs, e.g. "hgdp:5"; reps: timings per pass (median
# reported, default 5); views: "none" (default) or "none,train".
# Appends to <out dir>/kernels.csv, with the 1-minute load average.
args <- commandArgs(TRUE)
lib <- args[1]
label <- args[2]
out <- args[3]
cases <- args[4]
n_cores <- if (length(args) >= 5) as.integer(args[5]) else 1L
reps <- if (length(args) >= 6) as.integer(args[6]) else 5L
views <- if (length(args) >= 7) strsplit(args[7], ",")[[1]] else "none"
source("data-raw/radmix/newspeed/common.R")
ns <- tgc_load(lib)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
has_fused <- exists("radmix_update_q_loglik_cpp", envir = ns)

for (cs in tgc_cases(cases)) {
  x <- tgc_read(cs$data)
  dat <- ns$radmix_prepare(x)
  k <- cs$k
  set.seed(1)
  q <- matrix(stats::runif(k * dat$n_indiv), k)
  q <- sweep(q, 2, colSums(q), "/")
  p <- matrix(stats::runif(k * dat$n_loci), k)
  for (i in 1:5) {
    em <- ns$radmix_em_cpp(dat, q, p, n_cores)
    q <- em$q
    p <- em$p
  }
  q <- ns$radmix_project_cpp(q, TRUE)
  p <- ns$radmix_project_cpp(p, FALSE)
  for (view in views) {
    d <- switch(view,
      none = dat,
      train = ns$radmix_set_cv(dat, list(
        mode = "train", fold = 1, n_folds = 5, seed = 1
      ))
    )
    kernels <- list(
      loglik = function() ns$radmix_loglik_cpp(d, q, p, n_cores),
      update_q = function() ns$radmix_update_q_cpp(d, q, p, n_cores),
      update_p = function() ns$radmix_update_p_cpp(d, q, p, n_cores)
    )
    if (has_fused) {
      kernels$update_q_loglik <- function() {
        ns$radmix_update_q_loglik_cpp(d, q, p, n_cores)
      }
    }
    for (kernel in names(kernels)) {
      invisible(kernels[[kernel]]())
      times <- replicate(reps, system.time(kernels[[kernel]]())[["elapsed"]])
      tgc_append(data.frame(
        label = label, data = cs$data, k = k, n_cores = n_cores,
        view = view, kernel = kernel, reps = reps,
        seconds = stats::median(times), min_seconds = min(times),
        load = tgc_load_avg(),
        time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      ), file.path(out, "kernels.csv"))
    }
  }
}

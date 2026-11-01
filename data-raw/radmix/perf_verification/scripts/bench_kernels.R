# Per-kernel times of a build of tidygenclust: median of `reps` calls of
# each kernel, with fixed Q and P (random, seed 1).
#
#   Rscript bench_kernels.R <lib> <label> <out dir> [n_cores] [cases] [reps]
#
# cases: dataset:K pairs, default "eur:4,hgdp_quick:5,hgdp_quick:10" (those
# of the laptop). Appends to <out>/kernels.csv and saves the outputs of the
# kernels in <out>/kernels_<label>.rds (with n_cores = 1; otherwise
# kernels_<label>_<n_cores>.rds), to compare the builds (summarise.R).
args <- commandArgs(TRUE)
lib <- args[1]
label <- args[2]
out <- args[3]
n_cores <- if (length(args) >= 4) as.integer(args[4]) else 1L
cases <- if (length(args) >= 5 && nzchar(args[5])) {
  args[5]
} else {
  "eur:4,hgdp_quick:5,hgdp_quick:10"
}
reps <- if (length(args) >= 6) as.integer(args[6]) else 5L
source("data-raw/radmix/perf_verification/scripts/common.R")
ns <- tgc_load(lib)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

cases <- tgc_cases(cases)
rows <- list()
outputs <- list()
for (d in unique(vapply(cases, `[[`, "", "data"))) {
  x <- tgc_read(d)
  data <- ns$radmix_prepare(x)
  data_cv <- ns$radmix_set_cv(
    data,
    list(mode = "train", fold = 1, n_folds = 5, seed = 1)
  )
  for (cs in Filter(function(z) z$data == d, cases)) {
    k <- cs$k
    n_i <- data$n_indiv
    n_j <- data$n_loci
    set.seed(1)
    q <- matrix(runif(k * n_i), k, n_i)
    q <- sweep(q, 2, colSums(q), "/")
    p <- matrix(runif(k * n_j, 0.05, 0.95), k, n_j)
    kernels <- list(
      loglik = function() ns$radmix_loglik_cpp(data, q, p, n_cores),
      em = function() ns$radmix_em_cpp(data, q, p, n_cores),
      update_q = function() ns$radmix_update_q_cpp(data, q, p, n_cores),
      update_p = function() ns$radmix_update_p_cpp(data, q, p, n_cores),
      update_q_cv = function() ns$radmix_update_q_cpp(data_cv, q, p, n_cores),
      update_p_cv = function() ns$radmix_update_p_cpp(data_cv, q, p, n_cores)
    )
    for (kn in names(kernels)) {
      outputs[[paste(d, k, kn)]] <- kernels[[kn]]()
      rows[[length(rows) + 1]] <- data.frame(
        label = label, data = d, k = k, kernel = kn, n_cores = n_cores,
        seconds = tgc_median_time(kernels[[kn]], reps)
      )
    }
  }
}
tgc_append(do.call(rbind, rows), file.path(out, "kernels.csv"))
saveRDS(outputs, file.path(out, paste0(
  "kernels_", label, if (n_cores > 1) paste0("_", n_cores), ".rds"
)))

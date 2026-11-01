# Multi-core scaling of the kernels of a build of tidygenclust: the time of
# one pass of each kernel (median of `reps`), with fixed Q and P (random,
# seed 1), for each number of cores; and a check that the outputs are
# identical to those on 1 core (the kernels are meant to be
# thread-invariant).
#
#   Rscript scaling_kernels.R <lib> <label> <out dir> <dataset> <ks> <cores> [reps]
#
# e.g. Rscript scaling_kernels.R $TGC_WORK/libs/comb comb out hgdp 5,10 \
#        1,2,4,8,16,32 3
# Appends to <out>/scaling_kernels.csv: seconds, speed-up and efficiency
# against 1 core, the largest difference from the 1-core output (0 expected),
# and the OpenMP binding settings (OMP_PROC_BIND, OMP_PLACES) of the run.
args <- commandArgs(TRUE)
lib <- args[1]
label <- args[2]
out <- args[3]
dataset <- args[4]
source("data-raw/radmix/perf_verification/scripts/common.R")
ks <- tgc_ints(args[5])
cores <- tgc_ints(args[6])
reps <- if (length(args) >= 7) as.integer(args[7]) else 3L
ns <- tgc_load(lib)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
binding <- paste0(
  "bind=", Sys.getenv("OMP_PROC_BIND", "unset"),
  ";places=", Sys.getenv("OMP_PLACES", "unset")
)

x <- tgc_read(dataset)
data <- ns$radmix_prepare(x)
data_cv <- ns$radmix_set_cv(
  data,
  list(mode = "train", fold = 1, n_folds = 5, seed = 1)
)
cores <- sort(unique(c(1L, cores)))
for (k in ks) {
  set.seed(1)
  q <- matrix(runif(k * data$n_indiv), k, data$n_indiv)
  q <- sweep(q, 2, colSums(q), "/")
  p <- matrix(runif(k * data$n_loci, 0.05, 0.95), k, data$n_loci)
  kernels <- list(
    loglik = function(n) ns$radmix_loglik_cpp(data, q, p, n),
    em = function(n) ns$radmix_em_cpp(data, q, p, n),
    update_q = function(n) ns$radmix_update_q_cpp(data, q, p, n),
    update_p = function(n) ns$radmix_update_p_cpp(data, q, p, n),
    update_q_cv = function(n) ns$radmix_update_q_cpp(data_cv, q, p, n),
    update_p_cv = function(n) ns$radmix_update_p_cpp(data_cv, q, p, n)
  )
  for (kn in names(kernels)) {
    ref <- kernels[[kn]](1L)
    t1 <- NA
    for (n in cores) {
      diff <- max(abs(unlist(kernels[[kn]](n)) - unlist(ref)))
      t <- tgc_median_time(function() kernels[[kn]](n), reps)
      if (n == 1) t1 <- t
      tgc_append(data.frame(
        label = label, data = dataset, n_indiv = data$n_indiv,
        n_loci = data$n_loci, k = k, kernel = kn, n_cores = n, seconds = t,
        speedup = t1 / t, efficiency = t1 / t / n, diff_vs_1core = diff,
        binding = binding
      ), file.path(out, "scaling_kernels.csv"))
      cat(sprintf(
        "%s %s K=%d %-12s %3d cores: %.4f s (x%.2f), diff %.3g\n",
        label, dataset, k, kn, n, t, t1 / t, diff
      ))
    }
  }
}

# Diagnostics of the scaling (phase E): times `passes` calls of one kernel
# of a build on one dataset (fixed Q and P, as in scaling_kernels.R), after
# loading the data, so that a profiler wrapped around this script (perf
# stat) sees mostly the kernel. Prints the median time of a pass.
#
#   Rscript kernel_loop.R <lib> <dataset> <k> <kernel> <n_cores> <passes>
args <- commandArgs(TRUE)
source("data-raw/radmix/perf_verification/scripts/common.R")
ns <- tgc_load(args[1])
x <- tgc_read(args[2])
data <- ns$radmix_prepare(x)
k <- as.integer(args[3])
kernel <- args[4]
n <- as.integer(args[5])
passes <- as.integer(args[6])
set.seed(1)
q <- matrix(runif(k * data$n_indiv), k, data$n_indiv)
q <- sweep(q, 2, colSums(q), "/")
p <- matrix(runif(k * data$n_loci, 0.05, 0.95), k, data$n_loci)
f <- switch(kernel,
  loglik = function() ns$radmix_loglik_cpp(data, q, p, n),
  em = function() ns$radmix_em_cpp(data, q, p, n),
  update_q = function() ns$radmix_update_q_cpp(data, q, p, n),
  update_p = function() ns$radmix_update_p_cpp(data, q, p, n)
)
t <- replicate(passes, system.time(f())[["elapsed"]])
cat(sprintf("%s %s K=%d %s n=%d: median %.4f s (min %.4f, max %.4f)\n",
  basename(args[1]), args[2], k, kernel, n, median(t), min(t), max(t)))

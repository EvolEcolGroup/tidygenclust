# Instructions and cycles of the passes over the data of gt_radmix(), with
# perf, to compare builds when the machine is not quiet enough for timings
#
# On the data of benchmark_gl.R (200 individuals and 50,000 loci from K = 3
# ancestral populations, with genotype likelihoods from reads at a mean
# depth of 4), runs one pass (PASS: loglik, em, update_q or update_p) `reps`
# times, or a full fit (PASS = fit, K = 3, seed 1), with one thread, on the
# genotypes or on the likelihoods (DATA: genotypes or likelihoods), and has
# perf count only that part (not the start of R, nor the simulation),
# through perf's control FIFO. Instruction counts do not depend on the load
# of the machine, and cycles are only counted while the process runs (they
# still depend somewhat on the caches, which other processes share).
#
# The likelihoods are given by the example source of tidygenclust
# (src/radmix_gl_source.cpp, which keeps them in memory as doubles), if the
# build of tidygenclust has it; `dev` reads genotypes only.
#
# Needs kernel.perf_event_paranoid <= 2 (sudo sysctl
# kernel.perf_event_paranoid=2), and an optimised build of tidygenclust in a
# library of its own (TIDYGENCLUST_LIB). From the root of the package:
#   mkfifo ctl ack
#   PASS=loglik DATA=genotypes PERF_CTL=ctl PERF_ACK=ack \
#     TIDYGENCLUST_LIB=<lib> perf stat -e instructions:u,cycles:u -x, \
#     -o <out> --delay=-1 --control fifo:ctl,ack \
#     Rscript data-raw/radmix/gl/perf_gl.R
# perf writes the counts to <out>.

pass <- Sys.getenv("PASS", "loglik")
data_kind <- Sys.getenv("DATA", "genotypes")
reps <- if (pass == "fit") 1 else 20
n_indiv <- 200
n_loci <- 50000
k <- 3

lib <- Sys.getenv("TIDYGENCLUST_LIB")
if (nzchar(lib)) {
  library(tidygenclust, lib.loc = lib)
} else {
  library(tidygenclust)
}
ns <- asNamespace("tidygenclust")

# the data of benchmark_gl.R (the same seeds)
source(file.path("tests", "testthat", "helper-radmix-data.R"))
sim <- radmix_simulate_admixed(
  n = n_indiv, m = n_loci, k = k, missing = 0.01, seed = 1
)
if (data_kind == "genotypes") {
  dat <- ns$radmix_prepare(radmix_test_gt(sim$genotypes))
} else {
  # the example source of tidygenclust, through the helper of the tests
  # (which calls an internal function of the package)
  radmix_gl_dense_source_cpp <- ns$radmix_gl_dense_source_cpp
  dat <- radmix_test_gl(radmix_simulate_gl(sim$genotypes, depth = 4))
}
rm(sim)
invisible(gc())

# the starting values of benchmark_gl.R
set.seed(42)
q <- matrix(stats::runif(k * n_indiv), k, n_indiv)
q <- sweep(q, 2, colSums(q), "/")
p <- matrix(stats::runif(k * n_loci, 0.05, 0.95), k, n_loci)
f <- switch(pass,
  loglik = function() ns$radmix_loglik_cpp(dat, q, p, 1),
  em = function() ns$radmix_em_cpp(dat, q, p, 1),
  update_q = function() ns$radmix_update_q_cpp(dat, q, p, 1),
  update_p = function() ns$radmix_update_p_cpp(dat, q, p, 1),
  fit = function() {
    set.seed(1)
    ns$radmix_fit(dat, k = k, n_cores = 1)
  }
)
# a warm-up (the pages of the backing file, the caches)
if (pass != "fit") invisible(f())

# the measured part, between perf's enable and disable
ctl <- Sys.getenv("PERF_CTL")
ack <- Sys.getenv("PERF_ACK")
perf_command <- function(cmd) {
  if (!nzchar(ctl)) return(invisible())
  con <- file(ctl, "w", raw = TRUE)
  writeLines(cmd, con)
  close(con)
  con <- file(ack, "r", raw = TRUE)
  readLines(con, n = 1)
  close(con)
}
perf_command("enable")
for (r in seq_len(reps)) res <- f()
perf_command("disable")
if (pass == "fit") message("iterations: ", res$n_iter)

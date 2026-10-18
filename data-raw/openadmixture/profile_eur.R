# Time the individual C++ kernels on EUR_subset, K = 4.
# Requires the EUR_subset data: run data-raw/openadmixture/download_eur_subset.R first.
# Run from the root of the package, with tidygenclust installed (timings
# from devtools::load_all() are misleading, as it compiles without
# optimisation).
library(tidygenclust)
ns <- asNamespace("tidygenclust")
x <- tidypopgen::gt_load("data-raw/openadmixture/data/eur.gt")
dat <- ns$openadmixture_attach(ns$openadmixture_prepare(x))
set.seed(1)
k <- 4
q <- matrix(runif(k * dat$n_indiv), k)
q <- sweep(q, 2, colSums(q), "/")
p <- matrix(runif(k * dat$n_loci), k)
for (i in 1:5) {
  em <- ns$openadmixture_em_cpp(dat$geno_ptr, q, p, 1)
  q <- em$q
  p <- em$p
}
for (nc in c(1, 4)) {
  cat("n_cores", nc, "\n")
  print(system.time(ns$openadmixture_loglik_cpp(dat$geno_ptr, q, p, nc)))
  print(system.time(ns$openadmixture_grad_q_cpp(dat$geno_ptr, q, p, nc)))
  print(system.time(q1 <- ns$openadmixture_update_q_cpp(dat$geno_ptr, q, p, nc)))
  print(system.time(ns$openadmixture_grad_p_cpp(dat$geno_ptr, q1, p, nc)))
  print(system.time(ns$openadmixture_update_p_cpp(dat$geno_ptr, q1, p, nc)))
}

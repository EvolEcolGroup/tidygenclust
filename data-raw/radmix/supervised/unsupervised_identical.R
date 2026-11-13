# Check that the unsupervised path of gt_radmix() is bit-identical to a
# baseline build (the supervised mode must not change it).
#
# Usage (each library in its own R process):
#   Rscript unsupervised_identical.R <library> <output.rds> <scratch dir>
# then compare the two outputs with identical() (see the end of the file).
# The fits cover random starts with EM, SVD starts, init_q/init_p with
# fix_p, cross-validation, missing data and pseudohaploids, on the example
# data of tidypopgen and on EUR_subset (K = 3), with 1 and 2 threads.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
out_file <- args[2]
scratch <- args[3]
library(tidygenclust, lib.loc = lib)

eur_bed <- "/home/andrea/.cache/R/tidygenclust/EUR_subset.bed"
eur <- tidypopgen::gen_tibble(eur_bed,
  backingfile = tempfile("eur_", tmpdir = scratch), quiet = TRUE
)
pop_a <- tidypopgen::gen_tibble(
  system.file("extdata", "pop_a.bed", package = "tidypopgen"),
  backingfile = tempfile("pop_a_", tmpdir = scratch), quiet = TRUE
)
# a version of pop_a with some individuals pseudohaploid
pop_ph <- tidypopgen::gen_tibble(
  system.file("extdata", "pop_a.bed", package = "tidypopgen"),
  backingfile = tempfile("pop_ph_", tmpdir = scratch), quiet = TRUE
)
pop_ph <- tidypopgen::gt_pseudohaploid(pop_ph, test_n_loci = 10)
eur_imp <- tidypopgen::gt_impute_simple(
  tidypopgen::gen_tibble(eur_bed,
    backingfile = tempfile("eur_imp_", tmpdir = scratch), quiet = TRUE
  ),
  method = "mode"
)

res <- list()
res$eur_k3 <- gt_radmix(eur, k = 3, seed = 1, n_cores = 1)
res$eur_k3_2cores <- gt_radmix(eur, k = 3, seed = 1, n_cores = 2)
res$eur_k2_cv <- gt_radmix(eur, k = 2, seed = 2, crossval = TRUE,
  n_cores = 2)
res$pop_a <- gt_radmix(pop_a, k = 2:3, n_runs = 2, seed = 1:4)
res$pop_a_cv <- gt_radmix(pop_a, k = 2, seed = 3, crossval = TRUE)
res$pop_ph <- gt_radmix(pop_ph, k = 2, seed = 4)
res$eur_imp_svd <- gt_radmix(eur_imp, k = 2, seed = 5, init = "svd")
p_ref <- res$pop_a$P[[1]]
res$pop_a_fixp <- gt_radmix(pop_a, k = 2, init_p = p_ref, fix_p = TRUE,
  seed = 6)
res$pop_a_initq <- gt_radmix(pop_a, k = 2, init_q = res$pop_a$Q[[1]],
  seed = 7)
saveRDS(res, out_file)

# To compare:
# a <- readRDS("base.rds"); b <- readRDS("new.rds")
# for (n in names(a)) print(c(n, identical(a[[n]], b[[n]])))

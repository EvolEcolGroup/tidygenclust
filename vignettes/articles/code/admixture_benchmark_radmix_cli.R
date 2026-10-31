# gt_radmix() run from the command line by admixture_benchmark_cli.sh (in
# the same directory as this script), for the benchmark shown in the article
# admixture_benchmark.Rmd, as:
#
#   Rscript code/admixture_benchmark_radmix_cli.R <bed file> <K> <seed> \
#     <cores> <cross-validation (TRUE or FALSE)> <folds> <run directory>
#
# Like ADMIXTURE and fastmixture, each run is a new process: it starts R,
# loads tidygenclust, reads the bed file into a gen_tibble (with its backing
# file in the run directory), runs gt_radmix(), and saves the gt_admix object
# to an RDS file. It also writes the Q matrix and the log of the run, and
# prints the log-likelihood and the cross-validation error.

## ---- r-cli ----
args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]
k <- as.integer(args[2])
seed <- as.integer(args[3])
n_cores <- as.integer(args[4])
crossval <- as.logical(args[5])
n_folds <- as.integer(args[6])
run_dir <- args[7]
library(tidygenclust)
bigparallelr::set_blas_ncores(1)
hgdp <- tidypopgen::gen_tibble(
  bed_file,
  backingfile = file.path(run_dir, "hgdp_bench"),
  quiet = TRUE
)
fit <- gt_radmix(
  hgdp,
  k = k,
  seed = seed,
  crossval = crossval,
  n_folds = n_folds,
  n_cores = n_cores
)
saveRDS(fit, file.path(run_dir, "radmix.rds"))
# the Q matrix (one row per individual, with 6 decimals, as ADMIXTURE writes
# it) and the log
q <- fit$Q[[1]]
write.table(
  matrix(sprintf("%.6f", q), nrow = nrow(q)),
  file.path(run_dir, "radmix.Q"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
writeLines(fit$log[[1]], file.path(run_dir, "radmix.log"))
cat(sprintf("Loglikelihood: %.6f\n", fit$loglik))
if (crossval) {
  cat(sprintf("CV error (K=%d): %.6f\n", k, fit$cv))
}

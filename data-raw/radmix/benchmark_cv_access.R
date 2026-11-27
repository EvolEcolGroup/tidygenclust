# Time the C++ kernels of gt_radmix() with and without a cross-validation
# view (training and held-out), on one core, at K = 5, on the HGDP data of
# the benchmark article: the quick subset (209 x 25,000; also a full fit with
# 5-fold cross-validation) or the full data (1043 x 130,943; kernels only).
# The estimates are saved, so that two versions of the code can be checked
# to give identical results.
#
# Run from the root of the package, with tidygenclust installed with
# R CMD INSTALL --preclean . (not loaded with load_all()), on an otherwise
# idle machine:
#
#   Rscript data-raw/radmix/benchmark_cv_access.R <label> <output directory> \
#     [quick|full]
#
# It writes <output directory>/<label>_<data>_times.csv and _fits.rds.

suppressPackageStartupMessages({
  library(tidypopgen)
  library(tidygenclust)
})
ns <- asNamespace("tidygenclust")
args <- commandArgs(trailingOnly = TRUE)
label <- args[1]
out_dir <- args[2]
data_set <- if (length(args) >= 3) args[3] else "quick"
quick <- data_set == "quick"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

bed <- if (quick) "hgdp_bench_quick.bed" else "hgdp_bench.bed"
x <- gen_tibble(file.path("vignettes/articles/data/hgdp", bed),
  backingfile = tempfile("hgdp"), quiet = TRUE
)
dat <- ns$radmix_prepare(x)
k <- 5
set.seed(1)
q <- matrix(stats::runif(k * dat$n_indiv), k)
q <- sweep(q, 2, colSums(q), "/")
p <- matrix(stats::runif(k * dat$n_loci), k)
em <- ns$radmix_em_cpp(dat, q, p, 1)
q <- em$q
p <- em$p

# the views: all the genotypes, fold 1 of 5 hidden (training), and only
# fold 1 of 5 (held out)
views <- list(
  none = dat,
  train = ns$radmix_set_cv(dat, list(mode = "train", fold = 1, n_folds = 5,
    seed = 1)),
  held_out = ns$radmix_set_cv(dat, list(mode = "held_out", fold = 1,
    n_folds = 5, seed = 1))
)
# the median of 5 (3 on the full data) timings of a kernel
time_of <- function(f) {
  stats::median(replicate(if (quick) 5 else 3, system.time(f())[["elapsed"]]))
}
times <- list()
results <- list()
for (view in names(views)) {
  d <- views[[view]]
  kernels <- list(
    loglik = function() ns$radmix_loglik_cpp(d, q, p, 1),
    update_q = function() ns$radmix_update_q_cpp(d, q, p, 1),
    update_p = function() ns$radmix_update_p_cpp(d, q, p, 1)
  )
  for (kernel in names(kernels)) {
    results[[paste(view, kernel)]] <- kernels[[kernel]]()
    times[[length(times) + 1]] <- data.frame(
      label = label, data = data_set, view = view, kernel = kernel,
      seconds = time_of(kernels[[kernel]])
    )
    Sys.sleep(1)
  }
}
# a full fit with 5-fold cross-validation (on the quick subset only)
if (quick) {
  fit_seconds <- system.time(
    fit <- gt_radmix(x, k = k, seed = 1, crossval = TRUE)
  )[["elapsed"]]
  times[[length(times) + 1]] <- data.frame(
    label = label, data = data_set, view = "fit", kernel = "gt_radmix_cv",
    seconds = fit_seconds
  )
  results$fit <- fit[c("Q", "P", "loglik", "cv")]
}
times <- do.call(rbind, times)
prefix <- file.path(out_dir, paste0(label, "_", data_set))
utils::write.csv(times, paste0(prefix, "_times.csv"), row.names = FALSE)
saveRDS(results, paste0(prefix, "_fits.rds"))
print(times)

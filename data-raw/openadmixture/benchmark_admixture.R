# Speed of gt_openadmixture() against ADMIXTURE 1.3.0 on EUR_subset
#
# K = 2 to 4 on the EUR_subset data (379 individuals, 54,051 SNPs), seed 1:
# - ADMIXTURE, single-threaded (on the laptop used for development, its
#   multithreading made it slower: e.g. 113 s with 1 thread against 177 s
#   with 4 threads for K = 3), about 5 minutes in total;
# - gt_openadmixture() with 1 thread, with its default convergence threshold
#   (relative change in log-likelihood below 1e-7) and with a threshold
#   matched to ADMIXTURE's (absolute change below 1e-4), and with 4 threads
#   at the default threshold.
# Each run is preceded by a pause, to limit the effect of thermal throttling
# on laptops. The results are written to data-raw/openadmixture/benchmark_admixture.csv.
#
# Requirements: the EUR_subset data (Rscript data-raw/openadmixture/download_eur_subset.R),
# ADMIXTURE (on the PATH, or set ADMIXTURE_BIN), and tidygenclust installed
# (timings from devtools::load_all() are misleading, as it compiles without
# optimisation).
#
# Run from the root of the package:
#   Rscript data-raw/openadmixture/benchmark_admixture.R

library(tidygenclust)

admixture_bin <- Sys.getenv("ADMIXTURE_BIN", Sys.which("admixture"))
if (!nzchar(admixture_bin)) stop("ADMIXTURE not found; set ADMIXTURE_BIN")
bed_file <- normalizePath(file.path("data-raw", "openadmixture", "data", "EUR_subset.bed"))
eur <- tidypopgen::gt_load(file.path("data-raw", "openadmixture", "data", "eur.gt"))
pause <- 20

# results are collected in an environment, so that record() can add to them
bench <- new.env()
bench$results <- list()
record <- function(method, k, n_cores, seconds, iterations, loglik) {
  row <- data.frame(
    method = method, k = k, n_cores = n_cores, seconds = round(seconds, 1),
    iterations = iterations, loglik = loglik
  )
  print(row, row.names = FALSE)
  bench$results[[length(bench$results) + 1]] <- row
}
n_iter <- function(adm) length(grep("^Iteration", adm$log[[1]]))

for (k in 2:4) {
  # gt_openadmixture(), default threshold, 1 and 4 threads
  for (n_cores in c(1, 4)) {
    Sys.sleep(pause)
    seconds <- system.time(
      adm <- gt_openadmixture(eur, k = k, seed = 1, n_cores = n_cores)
    )[["elapsed"]]
    record(
      "gt_openadmixture", k, n_cores, seconds, n_iter(adm), adm$loglik
    )
  }
  # gt_openadmixture() with ADMIXTURE's threshold (|change| < 1e-4), 1 thread
  Sys.sleep(pause)
  rtol <- 1e-4 / abs(adm$loglik)
  seconds <- system.time(
    adm <- gt_openadmixture(eur, k = k, seed = 1, n_cores = 1, rtol = rtol)
  )[["elapsed"]]
  record(
    "gt_openadmixture (matched threshold)", k, 1, seconds, n_iter(adm),
    adm$loglik
  )
  # ADMIXTURE, 1 thread, run in a temporary directory (it writes its output
  # files in the working directory)
  Sys.sleep(pause)
  out_dir <- tempfile("admixture_")
  dir.create(out_dir)
  old_wd <- setwd(out_dir)
  seconds <- system.time(
    adm_log <- system2(admixture_bin, c("-s", "1", bed_file, k),
      stdout = TRUE, stderr = TRUE
    )
  )[["elapsed"]]
  setwd(old_wd)
  unlink(out_dir, recursive = TRUE)
  # ADMIXTURE prints one line per iteration ("1 (EM)", "2 (QN/Block)", ...)
  # and the final "Loglikelihood: ..." line
  ll_line <- grep("^Loglikelihood:", adm_log, value = TRUE)
  record(
    "ADMIXTURE 1.3.0", k, 1, seconds,
    length(grep("^[0-9]+ \\(", adm_log)),
    as.numeric(sub(".*: ", "", ll_line))
  )
}

results <- do.call(rbind, bench$results)
cpu <- grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]
results$cpu <- trimws(sub(".*:", "", cpu))
results$date <- as.character(Sys.Date())
utils::write.csv(results, file.path("data-raw", "openadmixture", "benchmark_admixture.csv"),
  row.names = FALSE
)
print(results[, 1:6], row.names = FALSE)

# Benchmark of gt_openadmixture() on the HGDP data prepared by
# admixture_benchmark_data.R (in the same directory as this script), with
# the same settings as admixture_benchmark_cli.sh, which runs ADMIXTURE and
# fastmixture from the command line. The results are shown in the article
# admixture_benchmark.Rmd (in the parent directory of this script), which
# also shows the code of the chunks of this script (the sections that start
# with a "## ---- label" line, read with knitr::read_chunk()).
#
# Run from the directory of the article (vignettes/articles), with
# tidygenclust installed with R CMD INSTALL (devtools::load_all() compiles
# its C++ code without optimisation):
#
#   cd vignettes/articles
#   Rscript code/admixture_benchmark_openadmixture.R
#
# Run it on an otherwise idle machine with at least 8 cores. It saves the Q
# matrix, log and timing of each run in data/admixture_benchmark, which the
# article reads. `Rscript code/admixture_benchmark_openadmixture.R quick`
# runs everything on a small subset of the data, to check the script, and
# saves its results in data/admixture_benchmark_quick instead.

## ---- r-setup ----
start_time <- Sys.time()
quick <- identical(commandArgs(trailingOnly = TRUE), "quick")
if (quick) {
  bed_file <- "data/hgdp/hgdp_bench_quick.bed"
  out_dir <- "data/admixture_benchmark_quick"
} else {
  bed_file <- "data/hgdp/hgdp_bench.bed"
  out_dir <- "data/admixture_benchmark"
}
dir.create(file.path(out_dir, "Q"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), showWarnings = FALSE)

## ---- r-settings ----
# the same settings as in admixture_benchmark_cli.sh
k <- 5
seed <- 1
n_folds <- 5
cores <- c(1, 4, 8)

## ---- r-packages ----
library(tidygenclust)
# gt_openadmixture() runs in parallel with OpenMP (with n_cores threads), and
# only uses BLAS for small (k x k) systems: BLAS is kept to a single thread,
# so that a multithreaded BLAS does not use more cores than requested
bigparallelr::set_blas_ncores(1)

## ---- r-import ----
# the bed file is read into a gen_tibble (with a temporary backing file),
# which is timed on its own: this is done once, and the gen_tibble can then
# be used for any number of analyses
start <- proc.time()[["elapsed"]]
hgdp <- tidypopgen::gen_tibble(
  bed_file,
  backingfile = tempfile("hgdp_bench"),
  quiet = TRUE
)
import_seconds <- proc.time()[["elapsed"]] - start

## ---- r-runs ----
# the Q matrices are written as ADMIXTURE and fastmixture write theirs (one
# row per individual, with 6 decimals)
write_q <- function(q, file) {
  q_text <- matrix(sprintf("%.6f", q), nrow = nrow(q))
  write.table(q_text, file, quote = FALSE, row.names = FALSE, col.names = FALSE)
}
# each combination of number of cores and cross-validation (or not) is run
# once; its Q matrix and log are saved, and its log-likelihood,
# cross-validation error and number of iterations are recorded
timings <- list()
for (n_cores in cores) {
  for (crossval in c(FALSE, TRUE)) {
    run <- paste0("openadmixture_t", n_cores, if (crossval) "_cv" else "_nocv")
    message("Running ", run)
    start <- proc.time()[["elapsed"]]
    fit <- gt_openadmixture(
      hgdp,
      k = k,
      seed = seed,
      crossval = crossval,
      n_folds = n_folds,
      n_cores = n_cores
    )
    seconds <- proc.time()[["elapsed"]] - start
    write_q(fit$Q[[1]], file.path(out_dir, "Q", paste0(run, ".Q")))
    writeLines(fit$log[[1]], file.path(out_dir, "logs", paste0(run, ".log")))
    converged_line <- grep("^Converged in", fit$log[[1]], value = TRUE)[1]
    timings[[run]] <- data.frame(
      tool = "openadmixture",
      n_cores = n_cores,
      crossval = crossval,
      seconds = seconds,
      loglik = fit$loglik,
      cv_error = if (crossval) fit$cv else NA,
      iterations = as.integer(
        sub("^Converged in ([0-9]+) iterations.*", "\\1", converged_line)
      )
    )
  }
}
timings <- do.call(rbind, timings)

## ---- r-individuals ----
# the population and region (one of seven, with "East" misspelt as "Est" in
# the metadata) of each individual, in the order of the rows of the Q
# matrices (that of the fam file, as for ADMIXTURE and fastmixture)
populations <- read.delim("data/hgdp/hgdp650_id_pop.txt")
pop_rows <- match(hgdp$id, populations$Id)
individuals <- data.frame(
  id = hgdp$id,
  population = populations$population[pop_rows],
  region = sub("Est", "East", populations$Pop7Groups[pop_rows])
)

## ---- r-session ----
# the machine and the versions of the software used, and the time taken to
# read the data into a gen_tibble
cpu <- tryCatch(
  sub(
    "^.*:\\s*", "",
    grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1]
  ),
  error = function(e) NA_character_,
  warning = function(w) NA_character_
)
session <- data.frame(
  item = c(
    "date", "R", "platform", "BLAS", "tidygenclust", "tidypopgen", "CPU",
    "cores", "import_seconds"
  ),
  value = c(
    format(Sys.Date()),
    R.version.string,
    R.version$platform,
    sessionInfo()$BLAS,
    as.character(utils::packageVersion("tidygenclust")),
    as.character(utils::packageVersion("tidypopgen")),
    cpu,
    parallel::detectCores(),
    round(import_seconds, 2)
  )
)

## ---- r-save ----
write.csv(
  timings,
  file.path(out_dir, "timings_openadmixture.csv"),
  row.names = FALSE
)
write.csv(
  individuals,
  file.path(out_dir, "individuals.csv"),
  row.names = FALSE
)
write.csv(
  session,
  file.path(out_dir, "session_r.csv"),
  row.names = FALSE
)
print(timings)
message(
  "Results saved in ", out_dir, ", after ",
  format(round(difftime(Sys.time(), start_time, units = "mins"), 1))
)

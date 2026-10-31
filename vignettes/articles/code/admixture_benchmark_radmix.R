# Benchmark of gt_radmix() in an R session, on the HGDP data prepared by
# admixture_benchmark_data.R (in the same directory as this script), with
# the same settings as admixture_benchmark_cli.sh, which runs the programs
# from the command line (including gt_radmix(), with Rscript). Here, R is
# already running and tidygenclust is loaded, so the time of each run does
# not include starting them. The results are shown in the article
# admixture_benchmark.Rmd (in the parent directory of this script), which
# also shows the code of the chunks of this script (the sections that start
# with a "## ---- label" line, read with knitr::read_chunk()).
#
# Run from the directory of the article (vignettes/articles), with
# tidygenclust installed with R CMD INSTALL (devtools::load_all() compiles
# its C++ code without optimisation):
#
#   cd vignettes/articles
#   Rscript code/admixture_benchmark_radmix.R
#
# Run it on an otherwise idle machine with at least as many cores as the
# largest number in cores (below). It saves the Q matrix, log and timing of
# each run in data/admixture_benchmark, which the article reads (the
# gen_tibble and gt_admix object of each run are left in data/hgdp/runs).
# `Rscript code/admixture_benchmark_radmix.R quick` runs everything on a
# small subset of the data, to check the script, and saves its results in
# data/admixture_benchmark_quick instead.

## ---- r-setup ----
start_time <- Sys.time()
quick <- identical(commandArgs(trailingOnly = TRUE), "quick")
if (quick) {
  bed_file <- "data/hgdp/hgdp_bench_quick.bed"
  out_dir <- "data/admixture_benchmark_quick"
  runs_dir <- "data/hgdp/runs_quick"
} else {
  bed_file <- "data/hgdp/hgdp_bench.bed"
  out_dir <- "data/admixture_benchmark"
  runs_dir <- "data/hgdp/runs"
}
dir.create(file.path(out_dir, "Q"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), showWarnings = FALSE)

## ---- r-settings ----
# the same settings as in admixture_benchmark_cli.sh
k <- 5
seed <- 1
n_folds <- 5
cores <- c(1, 2, 4)

## ---- r-packages ----
library(tidygenclust)
# gt_radmix() runs in parallel with OpenMP (with n_cores threads), and
# only uses BLAS for small (k x k) systems: BLAS is kept to a single thread,
# so that a multithreaded BLAS does not use more cores than requested
bigparallelr::set_blas_ncores(1)

## ---- r-runs ----
# the Q matrices are written as ADMIXTURE and fastmixture write theirs (one
# row per individual, with 6 decimals)
write_q <- function(q, file) {
  q_text <- matrix(sprintf("%.6f", q), nrow = nrow(q))
  write.table(q_text, file, quote = FALSE, row.names = FALSE, col.names = FALSE)
}
# each combination of number of cores and cross-validation (or not) is run
# once, in its own directory. As ADMIXTURE and fastmixture read the bed file
# and write their estimates, the time of each run includes reading the bed
# file into a gen_tibble (with its backing file in the directory of the
# run), and saving the gt_admix object (with Q and P) to an RDS file. Its Q
# matrix and log are then saved with the results, and its log-likelihood,
# cross-validation error and number of iterations are recorded
timings <- list()
for (n_cores in cores) {
  for (crossval in c(FALSE, TRUE)) {
    run <- paste0("radmix_t", n_cores, if (crossval) "_cv" else "_nocv")
    run_dir <- file.path(runs_dir, run)
    unlink(run_dir, recursive = TRUE)
    dir.create(run_dir, recursive = TRUE)
    message("Running ", run)
    start <- proc.time()[["elapsed"]]
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
    seconds <- proc.time()[["elapsed"]] - start
    write_q(fit$Q[[1]], file.path(out_dir, "Q", paste0(run, ".Q")))
    writeLines(fit$log[[1]], file.path(out_dir, "logs", paste0(run, ".log")))
    converged_line <- grep("^Converged in", fit$log[[1]], value = TRUE)[1]
    timings[[run]] <- data.frame(
      tool = "radmix",
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
# matrices (that of the fam file, as for ADMIXTURE and fastmixture, and of
# the gen_tibble of the last run)
populations <- read.delim("data/hgdp/hgdp650_id_pop.txt")
pop_rows <- match(hgdp$id, populations$Id)
individuals <- data.frame(
  id = hgdp$id,
  population = populations$population[pop_rows],
  region = sub("Est", "East", populations$Pop7Groups[pop_rows])
)

## ---- r-session ----
# the machine and the versions of the software used
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
    "cores"
  ),
  value = c(
    format(Sys.Date()),
    R.version.string,
    R.version$platform,
    sessionInfo()$BLAS,
    as.character(utils::packageVersion("tidygenclust")),
    as.character(utils::packageVersion("tidypopgen")),
    cpu,
    parallel::detectCores()
  )
)

## ---- r-save ----
write.csv(
  timings,
  file.path(out_dir, "timings_radmix.csv"),
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

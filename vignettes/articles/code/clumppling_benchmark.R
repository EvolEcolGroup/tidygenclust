# Benchmark of gt_clumppling() against the python implementation of
# clumppling (version 2.0.10), on the example datasets of clumppling. The
# results are shown in the article clumppling_benchmark.Rmd (in the parent
# directory of this script), which also shows the code of the chunks of this
# script (the sections that start with a "## ---- label" line, read with
# knitr::read_chunk()).
#
# The python implementation is run through reticulate, with the wrapper used
# to create the reference results of the tests
# (data-raw/clumppling/py_clumppling.R), in the conda environment created by
# data-raw/clumppling/setup_clumppling_env.sh (clumppling_reference; set
# CLUMPPLING_REFERENCE_ENV to use another environment with the same
# versions of the python packages).
#
# Each implementation is timed in three settings: with the Q matrices
# already in memory, starting from Q files, and from the command line (in a
# new process, which also starts R or python and loads the packages).
#
# Run from the directory of the article (vignettes/articles), with
# tidygenclust installed with R CMD INSTALL (devtools::load_all() compiles
# its C++ code without optimisation), on an otherwise idle machine:
#
#   cd vignettes/articles
#   Rscript code/clumppling_benchmark.R
#
# It takes a few minutes, and saves its results in
# data/clumppling_benchmark, which the article reads.

## ---- clump-setup ----
library(tidygenclust)
start_time <- Sys.time()
out_dir <- "data/clumppling_benchmark"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# each implementation is timed n_repeats times on each dataset and in each
# setting, after a first run (not timed) that loads its code and gives the
# results compared
n_repeats <- 5
seed <- 1L

## ---- clump-python ----
# py_clumppling() runs the python implementation on a gt_admix object, with
# the same arguments as gt_clumppling(); the wrapper is sourced from the
# root of the package, where it finds its python module
withr::with_dir(
  "../..",
  source("data-raw/clumppling/py_clumppling.R")
)
# the same python module, returning python objects instead of converting
# them to R objects, to time the python code alone
py_module <- reticulate::import_from_path(
  "py_rclumppling",
  path = "../../data-raw/clumppling",
  convert = FALSE
)
# the executables of python (in the conda environment) and R, for the runs
# from the command line
python <- reticulate::py_exe()
rscript <- file.path(R.home("bin"), "Rscript")

## ---- clump-data ----
# the example datasets of clumppling, included in tidygenclust: 200 runs of
# ADMIXTURE (50 for each K from 2 to 5) on the Cape Verde data of Verdu et
# al. (2017), with 399 individuals (capeverde), and 60 runs of STRUCTURE (20
# for each of K = 17, 19 and 21) on the chicken data of Rosenberg et al.
# (2001), with 600 chickens
capeverde <- readRDS(
  system.file("extdata", "capeverde_admix.rds", package = "tidygenclust")
)
chicken <- readRDS(
  system.file("extdata", "chicken_admix.rds", package = "tidygenclust")
)
# the first n runs of each K of a gt_admix object
first_runs <- function(x, n) {
  keep <- sort(unlist(lapply(split(seq_along(x$k), x$k), utils::head, n)))
  x$k <- x$k[keep]
  x$Q <- x$Q[keep]
  x
}
# capeverde with 5, 10, 20 and 50 runs for each K, and chicken
datasets <- list(
  capeverde_5 = first_runs(capeverde, 5),
  capeverde_10 = first_runs(capeverde, 10),
  capeverde_20 = first_runs(capeverde, 20),
  capeverde_50 = capeverde,
  chicken_20 = chicken
)
# each dataset as a directory of Q files (one for each run, with one row for
# each individual and one column for each cluster, as written by ADMIXTURE),
# and as a python list of numpy arrays
write_q_files <- function(x, q_dir) {
  dir.create(q_dir)
  for (i in seq_along(x$Q)) {
    utils::write.table(
      as.matrix(x$Q[[i]]),
      file.path(q_dir, sprintf("run%03d.Q", i)),
      row.names = FALSE,
      col.names = FALSE
    )
  }
  q_dir
}
to_numpy <- function(x) {
  reticulate::r_to_py(lapply(x$Q, function(q) {
    reticulate::np_array(
      matrix(as.numeric(q), nrow = nrow(q)),
      dtype = "float64",
      order = "C"
    )
  }))
}
inputs <- lapply(names(datasets), function(name) {
  list(
    admix = datasets[[name]],
    q_dir = write_q_files(datasets[[name]], file.path(tempdir(), name)),
    q_py = to_numpy(datasets[[name]])
  )
})
names(inputs) <- names(datasets)

## ---- clump-compare ----
# compare the results of the two implementations: their modes (the runs of
# each mode, overall and for each K), their representatives, and, if the
# modes are the same, the largest differences between their aligned modes
# and their costs across K
compare_results <- function(r, py) {
  # the modes (named e.g. K2M1) of the values of K that match the pattern
  modes_of <- function(res, pattern) {
    modes <- res$mode_replicates[grep(pattern, names(res$mode_replicates))]
    lapply(modes, as.integer)
  }
  same_modes <- identical(modes_of(r, ""), modes_of(py, ""))
  same_modes_k <- Filter(function(k) {
    pattern <- paste0("^K", k, "M")
    identical(modes_of(r, pattern), modes_of(py, pattern))
  }, r$K_range)
  # (the elements of both lists are matched by name)
  max_diff <- function(a, b) {
    b <- b[names(a)]
    max(abs(unlist(lapply(a, as.matrix)) - unlist(lapply(b, as.matrix))))
  }
  data.frame(
    n_modes_r = length(r$mode_replicates),
    n_modes_python = length(py$mode_replicates),
    same_modes = same_modes,
    same_modes_k = paste(same_modes_k, collapse = ", "),
    same_representatives = identical(
      r$mode_stats$Representative,
      py$mode_stats$Representative
    ),
    max_diff_aligned = if (same_modes) {
      max_diff(r$aligned_modes, py$aligned_modes)
    } else {
      NA
    },
    max_diff_cost = if (same_modes) {
      max_diff(r$cost_acrossK, py$cost_acrossK)
    } else {
      NA
    }
  )
}

## ---- clump-settings ----
# run a command, stopping if it fails
run_command <- function(command, args) {
  status <- system2(command, args, stdout = FALSE, stderr = FALSE)
  if (status != 0) {
    stop("Failed: ", command, " ", paste(args, collapse = " "))
  }
}
# the timed runs of each implementation, in each setting, on an input (a
# dataset as a gt_admix object, a directory of Q files and numpy arrays)
run_r <- list(
  # the Q matrices are in a gt_admix object
  memory = function(input) {
    gt_clumppling(input$admix, seed = seed)
  },
  # the Q files are read with tidypopgen::read_q_files()
  files = function(input) {
    gt_clumppling(tidypopgen::read_q_files(input$q_dir), seed = seed)
  },
  # a new R process loads tidygenclust, reads the Q files, and saves the
  # results in an rds file
  command_line = function(input) {
    out <- tempfile(fileext = ".rds")
    on.exit(unlink(out))
    run_command(rscript, c("-e", shQuote(sprintf(
      paste(
        "library(tidygenclust);",
        "res <- gt_clumppling(tidypopgen::read_q_files('%s'), seed = %d);",
        "saveRDS(res, '%s')"
      ),
      input$q_dir, seed, out
    ))))
  }
)
run_python <- list(
  # the Q matrices are numpy arrays (converted from R before the timing),
  # and the results are python objects (not converted to R)
  memory = function(input) {
    py_module$clumppling_run(Q_list = input$q_py, random_seed = seed)
  },
  # the Q files are read by clumppling's own parser, which writes the
  # processed matrices to a directory and reads them back
  files = function(input) {
    processed_dir <- tempfile()
    on.exit(unlink(processed_dir, recursive = TRUE))
    py_module$clumppling_run(
      input_dir = input$q_dir,
      processed_dir = processed_dir,
      random_seed = seed
    )
  },
  # the command-line tool of clumppling, without figures, which reads the Q
  # files and writes its results to a directory
  command_line = function(input) {
    out <- tempfile()
    on.exit(unlink(out, recursive = TRUE))
    run_command(python, c(
      "-m", "clumppling", "-i", shQuote(input$q_dir), "-o", shQuote(out),
      "-f", "admixture", "-v", "False", "--random_seed", seed
    ))
  }
)

## ---- clump-runs ----
timings <- list()
agreement <- list()
for (name in names(inputs)) {
  input <- inputs[[name]]
  x <- input$admix
  message("Running ", name)
  # a first run of each implementation in memory, compared, and a first run
  # from the Q files, neither timed
  res_r <- suppressMessages(gt_clumppling(x, seed = seed))
  res_py <- py_clumppling(x, seed = seed)
  agreement[[name]] <- data.frame(
    dataset = name,
    n_runs = length(x$Q),
    n_individuals = nrow(x$Q[[1]]),
    k_values = paste(sort(unique(x$k)), collapse = ", "),
    compare_results(res_r, res_py)
  )
  invisible(suppressMessages(run_r$files(input)))
  invisible(run_python$files(input))
  # then n_repeats timed runs of each, in each setting
  for (setting in names(run_r)) {
    for (i in seq_len(n_repeats)) {
      seconds_r <- system.time(
        suppressMessages(run_r[[setting]](input))
      )[["elapsed"]]
      seconds_py <- system.time(run_python[[setting]](input))[["elapsed"]]
      timings[[paste(name, setting, i)]] <- data.frame(
        dataset = name,
        setting = setting,
        repeat_run = i,
        implementation = c("r", "python"),
        seconds = c(seconds_r, seconds_py)
      )
    }
  }
}
timings <- do.call(rbind, timings)
agreement <- do.call(rbind, agreement)

## ---- clump-startup ----
# the time taken to start R and load tidygenclust, and to start python and
# load the code of the command-line tool of clumppling
startup <- do.call(rbind, lapply(seq_len(n_repeats), function(i) {
  data.frame(
    repeat_run = i,
    implementation = c("r", "python"),
    seconds = c(
      system.time(
        run_command(rscript, c("-e", shQuote("library(tidygenclust)")))
      )[["elapsed"]],
      system.time(
        run_command(python, c("-c", shQuote("import clumppling.__main__")))
      )[["elapsed"]]
    )
  )
}))

## ---- clump-session ----
# the machine and the versions of the software used
cpu <- tryCatch(
  sub(
    "^.*:\\s*", "",
    grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1]
  ),
  error = function(e) NA_character_,
  warning = function(w) NA_character_
)
py_versions <- py_clumppling_versions()
session <- data.frame(
  item = c(
    "date", "R", "tidygenclust", "CPU",
    names(py_versions)
  ),
  value = c(
    format(Sys.Date()),
    R.version.string,
    as.character(utils::packageVersion("tidygenclust")),
    cpu,
    py_versions
  )
)

## ---- clump-save ----
write.csv(timings, file.path(out_dir, "timings.csv"), row.names = FALSE)
write.csv(agreement, file.path(out_dir, "agreement.csv"), row.names = FALSE)
write.csv(startup, file.path(out_dir, "startup.csv"), row.names = FALSE)
write.csv(session, file.path(out_dir, "session.csv"), row.names = FALSE)
print(agreement)
message(
  "Results saved in ", out_dir, ", after ",
  format(round(difftime(Sys.time(), start_time, units = "mins"), 1))
)

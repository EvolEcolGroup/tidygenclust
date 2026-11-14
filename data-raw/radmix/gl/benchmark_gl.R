# The cost of genotype likelihoods in gt_radmix(), and the path of genotypes
# before and after their introduction
#
# On simulated data (by default 200 individuals and 50,000 loci from K = 3
# ancestral populations; the genotypes, and genotype likelihoods from reads
# at a mean depth of 4), times each pass over the data (log-likelihood, EM
# step, block updates of Q and of P), as the median of `reps` runs after a
# warm-up, on a gen_tibble (bytes) and, if the build reads them, on the
# likelihoods, given by the example source of tidygenclust
# (src/radmix_gl_source.cpp, which keeps them in memory as doubles), for
# several numbers of threads; then times a full fit of each (K = 3, seed 1).
# A package that stores likelihoods gives them through its own source,
# whose cost depends on its storage.
#
# Run it with an optimised build of tidygenclust installed in a library of
# its own (never with devtools::load_all(), which compiles without
# optimisation), from the root of the package, e.g. for the branches dev and
# gl, each in its own R process:
#   TIDYGENCLUST_LIB=<lib_dev> Rscript data-raw/radmix/gl/benchmark_gl.R dev
#   TIDYGENCLUST_LIB=<lib_gl> Rscript data-raw/radmix/gl/benchmark_gl.R gl
# Each run appends its timings to data-raw/radmix/gl/results/benchmark_gl.csv,
# with the label given as argument and the load average of the machine.
# With NGSADMIX set to the path of NGSadmix, the likelihoods are also written
# to a beagle file, and NGSadmix (with its default convergence criteria, but
# no filter on the minor allele frequency, `-minMaf 0`, so that both programs
# fit the same loci) is timed on them, for the same K and the largest number
# of threads.
# The sizes can be changed with BENCH_N and BENCH_M, and the numbers of
# threads with BENCH_CORES (e.g. "1,4").
# Time on an otherwise idle machine (or count instructions and cycles with
# perf_gl.R).

args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args) >= 1) args[1] else "build"
n_indiv <- as.integer(Sys.getenv("BENCH_N", "200"))
n_loci <- as.integer(Sys.getenv("BENCH_M", "50000"))
cores <- as.integer(strsplit(Sys.getenv("BENCH_CORES", "1,4"), ",")[[1]])
reps <- 5
k <- 3

lib <- Sys.getenv("TIDYGENCLUST_LIB")
if (nzchar(lib)) {
  library(tidygenclust, lib.loc = lib)
} else {
  library(tidygenclust)
}
message("tidygenclust from ", find.package("tidygenclust"), " (", label, ")")
ns <- asNamespace("tidygenclust")

# the simulation helpers of the tests: radmix_simulate_admixed(),
# radmix_test_gt() and radmix_simulate_gl()
source(file.path("tests", "testthat", "helper-radmix-data.R"))
sim <- radmix_simulate_admixed(
  n = n_indiv, m = n_loci, k = k, missing = 0.01, seed = 1
)
datasets <- list(
  genotypes = ns$radmix_prepare(radmix_test_gt(sim$genotypes))
)
# The likelihoods, given by the example source of tidygenclust (see
# radmix_test_gl() in the helpers of the tests, which calls an internal
# function of the package). A build without genotype likelihoods (dev)
# does not have it.
gl_ok <- "gt_radmix_data" %in% getNamespaceExports(ns)
no_data_ll <- list(genotypes = 0)
if (gl_ok) {
  radmix_gl_dense_source_cpp <- ns$radmix_gl_dense_source_cpp
  gl_dat <- radmix_test_gl(radmix_simulate_gl(sim$genotypes, depth = 4))
  datasets$likelihoods <- gl_dat
  # the log-likelihoods of fits on likelihoods are reported in NGSadmix's
  # convention, which counts log(1/3) for each individual without data at a
  # locus
  no_data_ll$likelihoods <- log(1 / 3) * sum(!gl_dat$obs)
}
rm(sim)
invisible(gc())

# starting values: a few EM steps from a random start, as in a fit
set.seed(42)
q <- matrix(stats::runif(k * n_indiv), k, n_indiv)
q <- sweep(q, 2, colSums(q), "/")
p <- matrix(stats::runif(k * n_loci, 0.05, 0.95), k, n_loci)

passes <- list(
  loglik = ns$radmix_loglik_cpp,
  em = ns$radmix_em_cpp,
  update_q = ns$radmix_update_q_cpp,
  update_p = ns$radmix_update_p_cpp
)
load_avg <- function() as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1])

rows <- list()
for (kind in names(datasets)) {
  dat <- datasets[[kind]]
  for (nc in cores) {
    for (pass in names(passes)) {
      f <- passes[[pass]]
      invisible(f(dat, q, p, nc))
      times <- vapply(seq_len(reps), function(r) {
        system.time(f(dat, q, p, nc))[["elapsed"]]
      }, numeric(1))
      rows[[length(rows) + 1]] <- data.frame(
        label = label, data = kind, what = pass, n_cores = nc,
        seconds = stats::median(times), min_seconds = min(times),
        iterations = NA_integer_, loglik = NA_real_, load = load_avg()
      )
      message(kind, ", ", pass, ", ", nc, " cores: ", stats::median(times), " s")
    }
    # a full fit from a random start
    set.seed(1)
    t_fit <- system.time(fit <- ns$radmix_fit(dat, k = k, n_cores = nc))
    rows[[length(rows) + 1]] <- data.frame(
      label = label, data = kind, what = "fit", n_cores = nc,
      seconds = t_fit[["elapsed"]], min_seconds = t_fit[["elapsed"]],
      iterations = fit$n_iter, loglik = fit$loglik + no_data_ll[[kind]],
      load = load_avg()
    )
    message(
      kind, ", fit, ", nc, " cores: ", t_fit[["elapsed"]], " s, ",
      fit$n_iter, " iterations"
    )
  }
}
# NGSadmix on the same likelihoods (normalised to sum to 1, with 8
# significant digits), written in chunks of loci to limit memory
ngsadmix <- Sys.getenv("NGSADMIX")
if (gl_ok && nzchar(ngsadmix)) {
  beagle_file <- tempfile("bench_", fileext = ".beagle.gz")
  con <- gzfile(beagle_file, "w")
  writeLines(paste(
    c("marker", "allele1", "allele2", rep(paste0("Ind", seq_len(n_indiv) - 1),
      each = 3
    )),
    collapse = "\t"
  ), con)
  for (j0 in seq(1, n_loci, by = 5000)) {
    js <- j0:min(n_loci, j0 + 4999)
    e1 <- gl_dat$e1[, js, drop = FALSE]
    e2 <- gl_dat$e2[, js, drop = FALSE]
    total <- 1 + e1 + e2
    lik <- array(c(1 / total, e1 / total, e2 / total), c(dim(e1), 3))
    values <- matrix(
      sprintf("%.8g", aperm(lik, c(3, 1, 2))), length(js),
      byrow = TRUE
    )
    writeLines(paste(paste0("1_", js), 0, 3,
      apply(values, 1, paste, collapse = "\t"),
      sep = "\t"
    ), con)
  }
  close(con)
  # (NGSadmix is slow: it is only timed with the largest number of threads)
  for (nc in max(cores)) {
    prefix <- tempfile("ngsadmix_")
    t_ngs <- system.time(system2(ngsadmix,
      c(
        "-likes", beagle_file, "-K", k, "-minMaf", 0, "-seed", 1,
        "-P", nc, "-o", prefix
      ),
      stdout = FALSE, stderr = FALSE
    ))[["elapsed"]]
    like_line <- grep("best like=", readLines(paste0(prefix, ".log")),
      value = TRUE
    )
    rows[[length(rows) + 1]] <- data.frame(
      label = label, data = "likelihoods (NGSadmix)", what = "fit",
      n_cores = nc, seconds = t_ngs, min_seconds = t_ngs,
      iterations = as.integer(sub(
        ".* after ([0-9]+) iterations.*", "\\1", like_line
      )),
      loglik = as.numeric(sub(".*best like=([-0-9.e]+).*", "\\1", like_line)),
      load = load_avg()
    )
    message(
      "NGSadmix, fit, ", nc, " cores: ", t_ngs, " s, ", like_line
    )
    unlink(paste0(prefix, c(".log", ".qopt", ".fopt.gz")))
  }
  unlink(beagle_file)
}

results <- do.call(rbind, rows)
results$n_indiv <- n_indiv
results$n_loci <- n_loci
out_dir <- file.path("data-raw", "radmix", "gl", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, "benchmark_gl.csv")
# (appended to the timings of earlier runs)
all_results <- if (file.exists(out_file)) {
  rbind(utils::read.csv(out_file), results)
} else {
  results
}
utils::write.csv(all_results, out_file, row.names = FALSE)
print(results)

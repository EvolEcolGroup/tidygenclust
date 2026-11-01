# Checks that fits of a build on different numbers of cores are identical
# (the kernels are written to be thread-invariant: the tests require
# identical results on 1 and 2 cores), from the fits_<label>_<n_cores>_*.rds
# files that bench_fits.R saves.
#
#   Rscript check_thread_invariance.R <out dir> <label> [rtol, default 1e-07]
args <- commandArgs(TRUE)
out <- args[1]
label <- args[2]
rtol <- if (length(args) >= 3) args[3] else "1e-07"
files <- list.files(out, paste0("^fits_", label, "_[0-9]+_", rtol, "\\.rds$"))
n_cores <- as.integer(sub(paste0("^fits_", label, "_([0-9]+)_.*"), "\\1", files))
files <- files[order(n_cores)]
n_cores <- sort(n_cores)
ref <- readRDS(file.path(out, files[1]))
cat("reference:", files[1], "\n")
for (i in seq_along(files)[-1]) {
  f <- readRDS(file.path(out, files[i]))
  for (n in intersect(names(ref), names(f))) {
    cat(sprintf(
      "%3d cores vs %d, %-14s Q %.3g, P %.3g, loglik %.3g\n", n_cores[i],
      n_cores[1], n, max(abs(f[[n]]$q - ref[[n]]$q)),
      max(abs(f[[n]]$p - ref[[n]]$p)), f[[n]]$loglik - ref[[n]]$loglik
    ))
  }
}

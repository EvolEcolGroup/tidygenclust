# Run the per-file timing of run_files.R inside covr's instrumented session.
# Arguments: label, optimisation ("O0" = covr's defaults, "O1" = the flags of
# PR #37: -O1 plus -D_GLIBCXX_ASSERTIONS for C++).
W <- "/media/andrea/annex/tmp/claude/tidygenclust/covrtime"
.libPaths(c(file.path(W, "rlib"), .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
label <- args[1]
opt <- args[2]
loadNamespace("covr")
if (opt == "O1") {
  flags <- sub("-O0", "-O1", getOption("covr.flags"), fixed = TRUE)
  cxx <- startsWith(names(flags), "CXX")
  flags[cxx] <- paste(flags[cxx], "-D_GLIBCXX_ASSERTIONS")
  options(covr.flags = flags)
}
print(getOption("covr.flags"))
# a fresh copy of the sources (covr compiles in the source directory)
build <- file.path(W, paste0("build_covr_", label))
unlink(build, recursive = TRUE)
dir.create(build)
file.copy(list.files(file.path(W, "src"), full.names = TRUE, all.files = TRUE,
                     no.. = TRUE), build, recursive = TRUE)
Sys.setenv(
  TGC_LIB = "", TGC_TESTS = Sys.getenv("TGC_SMOKE_TESTS", file.path(W, "src", "tests", "testthat")),
  TGC_OUT = file.path(W, "results", paste0(label, ".csv"))
)
t0 <- proc.time()
cov <- covr::package_coverage(build,
  type = "none",
  code = sprintf('source("%s")', file.path(W, "run_files.R")),
  quiet = TRUE,
  install_path = file.path(W, paste0("covr_install_", label))
)
t <- proc.time() - t0
cat(sprintf("covr total: elapsed %.1f s\n", t[["elapsed"]]))
cat(sprintf("coverage %.2f%%\n", covr::percent_coverage(cov)))
# the copies are only needed for this run
unlink(build, recursive = TRUE)
unlink(file.path(W, paste0("covr_install_", label)), recursive = TRUE)

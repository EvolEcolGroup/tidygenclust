# Time each test file of tidygenclust in turn, in one R process.
# Environment variables:
#   TGC_LIB    library holding the build to test ("" when covr has already
#              loaded its instrumented build)
#   TGC_TESTS  directory of the test files
#   TGC_OUT    CSV file for the results
lib <- Sys.getenv("TGC_LIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
test_dir <- Sys.getenv("TGC_TESTS")
out <- Sys.getenv("TGC_OUT")
suppressMessages(library(tidygenclust))
cat("tidygenclust from", find.package("tidygenclust"), "\n")
files <- sort(list.files(test_dir, "^test.*\\.R$"), method = "radix")
res <- lapply(files, function(f) {
  invisible(gc())
  t0 <- proc.time()
  r <- testthat::test_file(file.path(test_dir, f),
    package = "tidygenclust", load_package = "installed",
    reporter = "silent", stop_on_failure = FALSE
  )
  t <- proc.time() - t0
  df <- as.data.frame(r)
  row <- data.frame(
    file = f, elapsed = t[["elapsed"]],
    cpu = t[["user.self"]] + t[["sys.self"]] + t[["user.child"]] +
      t[["sys.child"]],
    n = sum(df$nb), failed = sum(df$failed), skipped = sum(df$skipped),
    error = any(df$error)
  )
  cat(sprintf("%-32s %8.1f s  cpu %8.1f s  n %d\n", f, row$elapsed, row$cpu,
              row$n))
  row
})
res <- do.call(rbind, res)
write.csv(res, out, row.names = FALSE)
cat(sprintf("TOTAL elapsed %.1f s, cpu %.1f s\n", sum(res$elapsed),
            sum(res$cpu)))

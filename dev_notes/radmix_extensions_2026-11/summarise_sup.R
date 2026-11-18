# Summary of run_sup_timing.sh: median (min-max) time of each build, case
# and number of threads, the ratio to dev, and the agreement of the fits.
args <- commandArgs(trailingOnly = TRUE)
res <- utils::read.csv(args[1])
rds_dir <- args[2]
res$build <- sub("^lib_", "", res$build)
key <- unique(res[, c("case", "n_cores")])
rows <- list()
for (r in seq_len(nrow(key))) {
  sub <- res[res$case == key$case[r] & res$n_cores == key$n_cores[r], ]
  dev_med <- stats::median(sub$time_s[sub$build == "dev"])
  for (b in unique(sub$build)) {
    t <- sub$time_s[sub$build == b]
    rows[[length(rows) + 1]] <- data.frame(
      case = key$case[r], n_cores = key$n_cores[r], build = b,
      n = length(t), median_s = stats::median(t), min_s = min(t),
      max_s = max(t), ratio = stats::median(t) / dev_med,
      iterations = paste(unique(sub$iterations[sub$build == b]),
                         collapse = "/"),
      rejected = paste(unique(sub$rejected[sub$build == b]), collapse = "/")
    )
  }
}
tab <- do.call(rbind, rows)
print(tab, digits = 3, row.names = FALSE)
# agreement with dev (Q, P, log-likelihood)
cat("\nAgreement with dev (largest absolute differences):\n")
for (f in list.files(rds_dir, pattern = "_lib_dev\\.rds$")) {
  d <- readRDS(file.path(rds_dir, f))
  for (b in c("new", "allB")) {
    g <- file.path(rds_dir, sub("_lib_dev", paste0("_lib_", b), f))
    if (!file.exists(g)) next
    n <- readRDS(g)
    cat(sprintf("%-28s %-5s Q %.2e  P %.2e  loglik %.2e (relative %.1e)\n",
                f, b, max(abs(n$Q - d$Q)), max(abs(n$P - d$P)),
                n$loglik - d$loglik, (n$loglik - d$loglik) / abs(d$loglik)))
  }
}

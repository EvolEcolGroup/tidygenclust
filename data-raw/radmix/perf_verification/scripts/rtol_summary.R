# The convergence thresholds compared (as in section 2.2 of HANDOVER.md):
# fits at rtol 1e-7 and 1e-8 (from final_compare.sh) and at ADMIXTURE's
# absolute threshold (1e-4 = rtol x |loglik|, from bench_fits.R with that
# rtol), seed 1, 1 core. For each, the time, the iterations, how far the
# log-likelihood is below that of the fit at ADMIXTURE's threshold, and
# the largest difference of Q from that fit.
#
#   Rscript rtol_summary.R <final_compare dir> <abs threshold dir> [labels]
args <- commandArgs(TRUE)
d_final <- args[1]
d_abs <- args[2]
labels <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else c("base", "comb")
f <- rbind(
  read.csv(file.path(d_final, "fits.csv")),
  read.csv(file.path(d_abs, "fits.csv"))
)
f <- subset(f, n_cores == 1 & !crossval & label %in% labels)
fm <- aggregate(cbind(seconds, iterations, loglik) ~ label + data + k + rtol,
  f, median)
read_q <- function(dir, label, rtol) {
  file <- file.path(dir, paste0("fits_", label, "_1_", rtol, ".rds"))
  if (file.exists(file)) readRDS(file) else list()
}
rows <- list()
for (lab in labels) {
  q_abs <- c(read_q(d_abs, lab, 6.42e-12), read_q(d_abs, lab, 2.58e-11))
  q_7 <- read_q(d_final, lab, 1e-07)
  q_8 <- read_q(d_final, lab, 1e-08)
  for (case in names(q_abs)) {
    dk <- strsplit(case, " ")[[1]]
    r <- subset(fm, label == lab & data == dk[1] & k == as.integer(dk[2]))
    ref <- r[r$rtol < 1e-9, ]
    t7 <- r$seconds[r$rtol == 1e-7]
    for (i in order(-r$rtol)) {
      q <- switch(as.character(r$rtol[i]),
        "1e-07" = q_7[[case]]$q, "1e-08" = q_8[[case]]$q, q_abs[[case]]$q
      )
      rows[[length(rows) + 1]] <- data.frame(
        label = lab, case = case,
        threshold = if (r$rtol[i] < 1e-9) "abs 1e-4" else format(r$rtol[i]),
        seconds = r$seconds[i], time_vs_1e7 = round(r$seconds[i] / t7, 2),
        iterations = r$iterations[i],
        below_abs = signif(ref$loglik - r$loglik[i], 3),
        max_dq_from_abs = signif(max(abs(q - q_abs[[case]]$q)), 3)
      )
    }
  }
}
print(do.call(rbind, rows), row.names = FALSE)

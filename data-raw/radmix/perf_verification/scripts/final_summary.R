# Summary of an output directory of final_compare.sh (or of any directory
# with kernels.csv and fits.csv): medians, and ratios to a reference build.
#
#   Rscript final_summary.R <out dir> [reference label, default base]
args <- commandArgs(TRUE)
out <- args[1]
ref <- if (length(args) >= 2) args[2] else "base"

kern <- read.csv(file.path(out, "kernels.csv"))
km <- aggregate(seconds ~ label + data + k + kernel + n_cores, kern, median)
key <- function(d, cols) do.call(paste, d[cols])
kcols <- c("data", "k", "kernel", "n_cores")
base <- subset(km, label == ref)
km$ref_seconds <- base$seconds[match(key(km, kcols), key(base, kcols))]
km$ratio <- round(km$seconds / km$ref_seconds, 3)
w <- reshape(km[, c(kcols, "label", "ratio")],
  idvar = kcols, timevar = "label", direction = "wide"
)
w$ref_seconds <- base$seconds[match(key(w, kcols), key(base, kcols))]
cat("== kernels (median seconds of", ref, "and the ratios to it)\n")
print(w[order(w$n_cores, w$data, w$k, w$kernel), ], row.names = FALSE)

f <- read.csv(file.path(out, "fits.csv"))
fm <- aggregate(cbind(seconds, iterations, loglik, cv) ~ label + data + k +
  n_cores + crossval + rtol, f, median, na.action = na.pass)
fcols <- c("data", "k", "n_cores", "crossval")
base <- subset(fm, label == ref & rtol == 1e-7)
m <- match(key(fm, fcols), key(base, fcols))
fm$ratio_to_ref_1e7 <- round(fm$seconds / base$seconds[m], 3)
fm$loglik_minus_ref <- signif(fm$loglik - base$loglik[m], 3)
cat(
  "\n== fits, seed 1 (median seconds; ratio to", ref,
  "at the same cores and CV, rtol 1e-7)\n"
)
print(fm[
  order(fm$crossval, fm$n_cores, fm$rtol, fm$data, fm$k, fm$label),
  c(
    "data", "k", "n_cores", "crossval", "rtol", "label", "seconds",
    "iterations", "ratio_to_ref_1e7", "loglik_minus_ref", "cv"
  )
], row.names = FALSE)

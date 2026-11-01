# Compares two builds (labels A and B) in an output directory of compare.sh:
# median times and their ratios (B / A), and the largest differences between
# their outputs (0 = bit-identical).
#
#   Rscript summarise.R <out dir> <label A> <label B>
args <- commandArgs(TRUE)
out <- args[1]
a_lab <- args[2]
b_lab <- args[3]

kern <- read.csv(file.path(out, "kernels.csv"))
kern <- subset(kern, label %in% c(a_lab, b_lab) & n_cores == 1)
km <- aggregate(seconds ~ label + data + k + kernel, kern, median)
w <- reshape(km,
  idvar = c("data", "k", "kernel"), timevar = "label",
  direction = "wide"
)
w$ratio <- w[[paste0("seconds.", b_lab)]] / w[[paste0("seconds.", a_lab)]]
cat("== kernels, 1 core (median seconds; ratio = B / A)\n")
print(w[order(w$data, w$k, w$kernel), ], digits = 3, row.names = FALSE)

ka <- readRDS(file.path(out, paste0("kernels_", a_lab, ".rds")))
kb <- readRDS(file.path(out, paste0("kernels_", b_lab, ".rds")))
cat("\n== kernel outputs, largest absolute difference (0 = bit-identical)\n")
for (n in names(ka)) {
  cat(sprintf(
    "  %-26s %.3g\n", n,
    max(abs(unlist(ka[[n]]) - unlist(kb[[n]])))
  ))
}

if (file.exists(file.path(out, "fits.csv"))) {
  f <- read.csv(file.path(out, "fits.csv"))
  f <- subset(f, label %in% c(a_lab, b_lab))
  fm <- aggregate(
    cbind(seconds, iterations, loglik) ~ label + data + k + n_cores +
      crossval + rtol,
    f, median
  )
  w <- reshape(fm,
    idvar = c("data", "k", "n_cores", "crossval", "rtol"),
    timevar = "label", direction = "wide"
  )
  w$ratio <- w[[paste0("seconds.", b_lab)]] / w[[paste0("seconds.", a_lab)]]
  cat("\n== fits (median seconds; ratio = B / A)\n")
  print(w, digits = 10, row.names = FALSE)
  for (rds in list.files(out, paste0("^fits_", a_lab, "_[0-9].*rds$"))) {
    fb <- file.path(out, sub(
      paste0("^fits_", a_lab),
      paste0("fits_", b_lab), rds
    ))
    if (!file.exists(fb)) next
    fa <- readRDS(file.path(out, rds))
    fb <- readRDS(fb)
    cat("\n== fit outputs", rds, ": largest differences in Q, P and loglik\n")
    for (n in intersect(names(fa), names(fb))) {
      cat(sprintf(
        "  %-14s Q %.3g, P %.3g, loglik %.3g\n", n,
        max(abs(fa[[n]]$q - fb[[n]]$q)), max(abs(fa[[n]]$p - fb[[n]]$p)),
        fa[[n]]$loglik - fb[[n]]$loglik
      ))
    }
  }
}

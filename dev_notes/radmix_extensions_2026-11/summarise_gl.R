# Summary of run_gl_timing.sh: median over rounds of each pass and of the
# fit, by combination and threads, with the ratio to ext_base.
args <- commandArgs(trailingOnly = TRUE)
r <- utils::read.csv(args[1])
r$combo <- ifelse(r$source == "dense", paste0(r$label, "_dense"), r$label)
base <- if (length(args) >= 2) args[2] else "ext_base"
agg <- stats::aggregate(seconds ~ combo + what + n_cores, data = r,
                        FUN = stats::median)
rng <- stats::aggregate(seconds ~ combo + what + n_cores, data = r,
                        FUN = function(x) diff(range(x)) / stats::median(x))
agg$spread <- rng$seconds
ref <- agg[agg$combo == base, c("what", "n_cores", "seconds")]
names(ref)[3] <- "ref"
agg <- merge(agg, ref)
agg$ratio <- agg$seconds / agg$ref
for (nc in sort(unique(agg$n_cores))) {
  cat("\n== threads:", nc, "(ratio to", base, "; seconds of", base, ")\n")
  sub <- agg[agg$n_cores == nc, ]
  w <- stats::reshape(sub[, c("combo", "what", "ratio")], idvar = "combo",
                      timevar = "what", direction = "wide")
  names(w) <- sub("ratio.", "", names(w))
  print(w, digits = 3, row.names = FALSE)
  print(unique(sub[sub$combo == base, c("what", "ref")]), digits = 3,
        row.names = FALSE)
  cat("max spread (max-min)/median:", signif(max(sub$spread), 2), "\n")
}

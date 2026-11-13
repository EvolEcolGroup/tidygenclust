# Part B summary: passes, iterations and log-likelihoods of each start or
# warm-up, from the .rds files of warmup_compare.R.
# Usage: Rscript warmup_summary.R <warmup_dir> [<out_csv>]
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
out_csv <- if (length(args) >= 2) args[2] else file.path(dir, "warmup.csv")
files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
get <- function(r, nm) if (is.null(r[[nm]])) NA else r[[nm]]
tab <- do.call(rbind, lapply(files, function(f) {
  r <- readRDS(f)
  data.frame(
    dataset = r$dataset, k = r$k, seed = r$seed, method = r$method,
    passes = r$passes, warm_passes = r$warm_passes,
    fit_passes = r$fit_passes, fit_iter = r$n_iter, loglik = r$loglik,
    converged = r$converged, time = r$time,
    time_svd = get(r, "time_svd"), epochs = get(r, "epochs"),
    sub_iter = get(r, "sub_iter"), load1 = as.numeric(strsplit(r$load,
      " ")[[1]][1])
  )
}))
# the best log-likelihood found for each dataset and K (any method, any
# seed), and how far each run ended below it
best <- stats::aggregate(loglik ~ dataset + k, tab, max)
names(best)[3] <- "best"
tab <- merge(tab, best)
tab$below_best <- tab$best - tab$loglik
tab <- tab[order(tab$dataset, tab$k, tab$method, tab$seed), ]
utils::write.csv(tab, out_csv, row.names = FALSE)

fmt_tab <- function(d) {
  if (requireNamespace("knitr", quietly = TRUE)) {
    print(knitr::kable(d, row.names = FALSE, digits = 2))
  } else {
    print(d, row.names = FALSE)
  }
  cat("\n")
}
cat("## Per run\n\n")
fmt_tab(tab[, c("dataset", "k", "seed", "method", "passes", "warm_passes",
  "fit_iter", "below_best", "time", "load1")])
cat("## Per method (medians over seeds; passes relative to random, same seeds)\n\n")
s <- do.call(rbind, lapply(split(tab, list(tab$dataset, tab$k, tab$method),
  drop = TRUE), function(d) {
  ref <- tab[tab$dataset == d$dataset[1] & tab$k == d$k[1] &
    tab$method == "random" & tab$seed %in% d$seed, ]
  data.frame(
    dataset = d$dataset[1], k = d$k[1], method = d$method[1],
    n_seeds = nrow(d), passes = stats::median(d$passes),
    rel_passes = if (nrow(ref)) stats::median(d$passes) /
      stats::median(ref$passes) else NA,
    warm_passes = stats::median(d$warm_passes),
    fit_iter = stats::median(d$fit_iter),
    below_best = paste(round(sort(d$below_best), 1), collapse = ", "),
    time = stats::median(d$time)
  )
}))
fmt_tab(s)

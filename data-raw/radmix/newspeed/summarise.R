# Summary of the results of run_bench.sh (bench_fits.R, bench_kernels.R,
# bench_cv.R) as markdown tables. Fits: for each case (data, K, threads,
# iteration limit), the median time of each build, its ratio to the
# reference build, the iterations and passes over the genotypes, and whether
# the results agree with the reference (iterations, log-likelihoods) and
# across numbers of threads. Kernels: the median time of each pass, and the
# cost of an iteration made of them. CV: the cost of fold fits relative to
# plain fits.
#
#   Rscript data-raw/radmix/newspeed/summarise.R <out dir> [reference label]
#
# The reference is the first label of the results unless given.
args <- commandArgs(TRUE)
out <- args[1]
ref_arg <- if (length(args) >= 2) args[2] else NULL
read_res <- function(name) {
  f <- file.path(out, name)
  if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
}
fits <- read_res("fits.csv")
kern <- read_res("kernels.csv")
cvres <- read_res("cv.csv")

# Kernels: medians over the runs of the median of each run
if (!is.null(kern)) {
  ref <- if (is.null(ref_arg)) kern$label[1] else ref_arg
  cat("# Passes over the genotypes:", out, "\n\n")
  cat("Median seconds of one pass (medians over runs of medians); ",
    "iteration = the passes of an accepted iteration: before, 2 update_q + ",
    "2 update_p + loglik; fused, update_q + update_q_loglik + 2 update_p. ",
    "Ratio: to the iteration of the reference build `", ref, "` (before).\n\n",
    sep = ""
  )
  agg <- stats::aggregate(seconds ~ label + data + k + n_cores + view + kernel,
    data = kern, FUN = stats::median
  )
  wide <- stats::reshape(agg,
    idvar = c("label", "data", "k", "n_cores", "view"),
    timevar = "kernel", direction = "wide"
  )
  names(wide) <- sub("^seconds\\.", "", names(wide))
  if (is.null(wide$update_q_loglik)) wide$update_q_loglik <- NA
  wide$iter_before <- 2 * wide$update_q + 2 * wide$update_p + wide$loglik
  wide$iter_fused <- wide$update_q + wide$update_q_loglik + 2 * wide$update_p
  cat("| data | K | threads | view | build | loglik | update_q |",
    "update_p | update_q_loglik | iteration before | iteration fused |",
    "ratio fused |\n")
  cat("|---|---|---|---|---|---|---|---|---|---|---|---|\n")
  for (r in seq_len(nrow(wide))) {
    w <- wide[r, ]
    ref_iter <- wide$iter_before[wide$label == ref & wide$data == w$data &
      wide$k == w$k & wide$n_cores == w$n_cores & wide$view == w$view]
    cat(sprintf(
      "| %s | %d | %d | %s | %s | %.3f | %.3f | %.3f | %s | %.3f | %s | %s |\n",
      w$data, w$k, w$n_cores, w$view, w$label, w$loglik, w$update_q,
      w$update_p,
      if (is.na(w$update_q_loglik)) "" else sprintf("%.3f", w$update_q_loglik),
      w$iter_before,
      if (is.na(w$iter_fused)) "" else sprintf("%.3f", w$iter_fused),
      if (is.na(w$iter_fused) || !length(ref_iter)) "" else
        sprintf("%.3f", w$iter_fused / ref_iter[1])
    ))
  }
  cat("\n")
}

# CV: fold fits against plain fits of the same number of iterations
if (!is.null(cvres)) {
  cat("# Cross-validation fold fits:", out, "\n\n")
  cat("Median seconds of fits of a fixed number of iterations (each run ",
    "reports a median), plain view and training views; ratio: to the plain ",
    "view of the same run.\n\n",
    sep = ""
  )
  agg <- stats::aggregate(cbind(seconds, ratio_to_none) ~
    label + data + k + n_cores + max_iter + view + mask_mb,
  data = cvres, FUN = stats::median
  )
  cat("| data | K | threads | iterations | build | view | median s |",
    "ratio to plain | mask MB |\n")
  cat("|---|---|---|---|---|---|---|---|---|\n")
  for (r in seq_len(nrow(agg))) {
    a <- agg[r, ]
    cat(sprintf(
      "| %s | %d | %d | %d | %s | %s | %.2f | %.3f | %.1f |\n",
      a$data, a$k, a$n_cores, a$max_iter, a$label, a$view, a$seconds,
      a$ratio_to_none, a$mask_mb
    ))
  }
  cat("\n")
}

if (is.null(fits)) quit(save = "no")
labels <- unique(fits$label)
ref <- if (is.null(ref_arg)) labels[1] else ref_arg
fits$loglik <- as.numeric(fits$loglik)
keys <- c("data", "k", "n_cores", "max_iter", "rtol")
cases <- unique(fits[keys])
cases <- cases[do.call(order, cases), , drop = FALSE]

cat("# Fits:", out, "\n\n")
cat("Reference build: `", ref, "`. Times are medians over the repeats ",
  "(n) and seeds; ratio = median of the build / median of the reference ",
  "(below 1 is faster). Load: the range of the 1-minute load average ",
  "before the fits. Passes: over the genotypes in the main loop (the ",
  "start included, EM priming not), counted from the log, before / with ",
  "the fused pass.\n\n",
  sep = ""
)
cat("| data | K | threads | max_iter | build | n | median s | ratio |",
  "s / iteration | iterations | rejected | passes before / fused |",
  "load | max rel. loglik diff vs ref |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (r in seq_len(nrow(cases))) {
  cs <- cases[r, ]
  sel <- Reduce(`&`, lapply(keys, function(k) fits[[k]] == cs[[k]]))
  sub <- fits[sel, ]
  ref_sub <- sub[sub$label == ref, ]
  ref_med <- stats::median(ref_sub$seconds)
  for (lab in unique(sub$label)) {
    s <- sub[sub$label == lab, ]
    med <- stats::median(s$seconds)
    # the log-likelihoods of the same seeds, against the reference
    m <- merge(
      unique(s[c("seed", "loglik", "iterations")]),
      unique(ref_sub[c("seed", "loglik", "iterations")]),
      by = "seed", suffixes = c("", "_ref")
    )
    lldiff <- if (nrow(m)) {
      signif(max(abs(m$loglik - m$loglik_ref) / abs(m$loglik_ref)), 2)
    } else {
      NA
    }
    same_iter <- nrow(m) == 0 || all(m$iterations == m$iterations_ref)
    cat(sprintf(
      "| %s | %d | %d | %d | %s | %d | %.2f | %s | %.3f | %s%s | %s | %s / %s | %.1f-%.1f | %s |\n",
      cs$data, cs$k, cs$n_cores, cs$max_iter, lab, nrow(s), med,
      if (length(ref_med) && !is.na(ref_med)) sprintf("%.3f", med / ref_med) else "",
      stats::median(s$seconds_per_iter),
      paste(sort(unique(s$iterations)), collapse = ","),
      if (same_iter) "" else " (differ)",
      paste(sort(unique(s$rejected)), collapse = ","),
      paste(sort(unique(s$passes_before)), collapse = ","),
      paste(sort(unique(s$passes_fused)), collapse = ","),
      min(s$load_before), max(s$load_before), format(lldiff)
    ))
  }
}

# Thread invariance: for each build, case and seed, the log-likelihood
# histories and Q of the saved fits must be identical for every number of
# threads
cat("\n## Thread invariance (identical log-likelihood histories and Q)\n\n")
rds <- list.files(file.path(out, "fits"), pattern = "\\.rds$")
if (length(rds)) {
  # <label>_<data>_k<K>_s<seed>_c<n_cores>_i<max_iter>.rds
  parts <- regmatches(rds, regexec(
    "^(.*)_(eur|hgdp_quick|hgdp_full|hgdp)_k([0-9]+)_s([0-9]+)_c([0-9]+)_i([0-9]+)\\.rds$",
    rds
  ))
  info <- do.call(rbind, lapply(parts, function(p) {
    data.frame(
      file = p[1], label = p[2], data = p[3], k = p[4], seed = p[5],
      n_cores = p[6], max_iter = p[7]
    )
  }))
  groups <- split(info, info[c("label", "data", "k", "seed", "max_iter")],
    drop = TRUE
  )
  for (g in groups) {
    if (length(unique(g$n_cores)) < 2) next
    res <- lapply(g$file, function(f) readRDS(file.path(out, "fits", f)))
    same <- all(vapply(res[-1], function(z) {
      identical(z$loglik_history, res[[1]]$loglik_history) &&
        identical(z$q, res[[1]]$q)
    }, logical(1)))
    cat(sprintf(
      "- %s %s K = %s seed %s max_iter %s, threads %s: %s\n",
      g$label[1], g$data[1], g$k[1], g$seed[1], g$max_iter[1],
      paste(sort(as.integer(g$n_cores)), collapse = ","),
      if (same) "identical" else "DIFFERENT"
    ))
  }
}

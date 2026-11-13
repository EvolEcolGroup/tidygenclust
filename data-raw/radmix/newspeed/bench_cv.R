# The cost of the passes of a cross-validation fold fit, with one build of
# tidygenclust: fits of a fixed number of iterations (rtol = 0, so that every
# build does the same work), warm-started from the same estimates (a random
# start and 5 EM steps, seed 1, both blocks given, so no EM priming), on the
# plain view ("none"), and on the training view of fold 1 of 5, with the
# hash computed for every genotype ("train_hash", cv element mask = FALSE)
# and with the precomputed bit mask ("train_mask", the default in builds that
# have it; builds without it ignore the element and hash). The time of a
# fold fit includes building the mask.
#
#   Rscript data-raw/radmix/newspeed/bench_cv.R <lib> <label> <out dir> \
#     <cases> [n_cores] [max_iter] [reps]
#
# cases: dataset:K pairs, e.g. "hgdp:5"; max_iter default 4, reps default 3
# (the median is reported). Appends to <out dir>/cv.csv, with the size of the
# bit mask (one bit per genotype, rounded up to 64 individuals per locus).
args <- commandArgs(TRUE)
lib <- args[1]
label <- args[2]
out <- args[3]
cases <- args[4]
n_cores <- if (length(args) >= 5) as.integer(args[5]) else 1L
max_iter <- if (length(args) >= 6) as.integer(args[6]) else 4L
reps <- if (length(args) >= 7) as.integer(args[7]) else 3L
source("data-raw/radmix/newspeed/common.R")
ns <- tgc_load(lib)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

for (cs in tgc_cases(cases)) {
  x <- tgc_read(cs$data)
  dat <- ns$radmix_prepare(x)
  k <- cs$k
  set.seed(1)
  q <- matrix(stats::runif(k * dat$n_indiv), k)
  q <- sweep(q, 2, colSums(q), "/")
  p <- matrix(stats::runif(k * dat$n_loci), k)
  for (i in 1:5) {
    em <- ns$radmix_em_cpp(dat, q, p, n_cores)
    q <- em$q
    p <- em$p
  }
  # the user layout of starting values (individuals x K, loci x K)
  init_q <- t(ns$radmix_project_cpp(q, TRUE))
  init_p <- t(ns$radmix_project_cpp(p, FALSE))
  views <- list(
    none = dat,
    train_hash = ns$radmix_set_cv(dat, list(
      mode = "train", fold = 1, n_folds = 5, seed = 1, mask = FALSE
    )),
    train_mask = ns$radmix_set_cv(dat, list(
      mode = "train", fold = 1, n_folds = 5, seed = 1, mask = TRUE
    ))
  )
  mask_mb <- ceiling(dat$n_indiv / 64) * 8 * dat$n_loci / 2^20
  results <- list()
  for (r in seq_len(reps)) {
    # the views alternate within each repeat
    for (view in names(views)) {
      t <- system.time(fit <- ns$radmix_fit(views[[view]],
        k = k, init_q = init_q, init_p = init_p, max_iter = max_iter,
        rtol = 0, n_cores = n_cores
      ))[["elapsed"]]
      results[[length(results) + 1]] <- data.frame(
        view = view, seconds = t, loglik = fit$loglik
      )
    }
  }
  res <- do.call(rbind, results)
  for (view in names(views)) {
    s <- res[res$view == view, ]
    tgc_append(data.frame(
      label = label, data = cs$data, k = k, n_cores = n_cores,
      max_iter = max_iter, view = view, reps = reps,
      seconds = stats::median(s$seconds), min_seconds = min(s$seconds),
      ratio_to_none = stats::median(s$seconds) /
        stats::median(res$seconds[res$view == "none"]),
      loglik = sprintf("%.10f", s$loglik[1]), mask_mb = mask_mb,
      load = tgc_load_avg(), time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ), file.path(out, "cv.csv"))
  }
  # the two training views must give identical fits
  stopifnot(identical(
    res$loglik[res$view == "train_hash"], res$loglik[res$view == "train_mask"]
  ))
}

# SVD starts against random starts for gt_radmix(), on one dataset: time
# (SVD, start, fit), iterations and log-likelihood of each run, and how
# often each kind of start reaches the best log-likelihood found by any
# run at that K.
#
#   Rscript svd_vs_random.R <lib> <out dir> <dataset> <ks> <seeds> <n_cores> \
#     [rtol] [strategies]
#
# e.g. Rscript svd_vs_random.R $TGC_WORK/libs/comb out hgdp 3,5,7,10 1:5 8
# strategies (comma separated, default "random,svd,svd_call"):
# - random: gt_radmix(x, k, seed, n_cores, rtol), the default start;
# - svd: the steps of gt_radmix(init = "svd") timed one by one, with the
#   internal functions: the SVD once (max(ks) + 5 components), then for each
#   run set.seed(seed), the start (best of 3 ALS solutions) and the fit;
# - svd_call: one call of gt_radmix(x, k = ks, seed = 1, init = "svd"), timed
#   as a whole (the SVD shared by all the K), as a user would run it (its
#   time is in the column call_s of each of its rows).
# Missing genotypes are imputed (gt_impute_simple(method = "mean0")) before
# all the runs; the fits ignore the imputed values (use_imputed = FALSE).
# Rows are appended to <out>/svd_vs_random.csv after each run; Q of each
# run is saved in <out>/svd_vs_random_q_<dataset>.rds; a summary per K and
# strategy is written to <out>/svd_vs_random_summary_<dataset>.txt.
args <- commandArgs(TRUE)
lib <- args[1]
out <- args[2]
dataset <- args[3]
source("data-raw/radmix/perf_verification/scripts/common.R")
ks <- tgc_ints(args[4])
seeds <- eval(parse(text = args[5]))
n_cores <- as.integer(args[6])
rtol <- if (length(args) >= 7) as.numeric(args[7]) else 1e-7
strategies <- strsplit(
  if (length(args) >= 8) args[8] else "random,svd,svd_call", ","
)[[1]]
ns <- tgc_load(lib)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
csv <- file.path(out, "svd_vs_random.csv")
q_file <- file.path(out, paste0("svd_vs_random_q_", dataset, ".rds"))
qs <- if (file.exists(q_file)) readRDS(q_file) else list()
blas_threads <- if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_get_num_procs()
} else {
  NA
}
cat("BLAS:", extSoftVersion()[["BLAS"]], "threads:", blas_threads, "\n")

x <- tgc_read(dataset)
t_impute <- 0
if (any(tidypopgen::loci_missingness(x, n_cores = n_cores) > 0)) {
  t_impute <- system.time(
    x <- tidypopgen::gt_impute_simple(x, method = "mean0", n_cores = n_cores)
  )[["elapsed"]]
}
data <- ns$radmix_prepare(x)
cat(dataset, ":", data$n_indiv, "individuals,", data$n_loci, "loci\n")

record <- function(strategy, k, seed, svd_s, start_s, fit_s, fit,
                   call_s = NA) {
  qs[[paste(strategy, k, seed, rtol, sep = "_")]] <<- unclass(fit$q)
  saveRDS(qs, q_file)
  row <- data.frame(
    dataset = dataset, n_indiv = data$n_indiv, n_loci = data$n_loci,
    strategy = strategy, k = k, seed = seed, n_cores = n_cores, rtol = rtol,
    impute_s = t_impute, svd_s = svd_s, start_s = start_s, fit_s = fit_s,
    call_s = call_s,
    iterations = fit$n_iter, loglik = fit$loglik, blas_threads = blas_threads
  )
  tgc_append(row, csv)
  print(row[, c("strategy", "k", "seed", "svd_s", "start_s", "fit_s",
                "iterations", "loglik")], row.names = FALSE)
}

if ("svd" %in% strategies) {
  t_svd <- system.time(
    factors <- ns$radmix_svd_factors(x, rank = max(ks) + 5, n_cores = n_cores)
  )[["elapsed"]]
  cat("SVD:", t_svd, "s\n")
}
for (k in ks) {
  for (seed in seeds) {
    if ("random" %in% strategies) {
      Sys.sleep(2)
      t <- system.time(a <- gt_radmix(x,
        k = k, seed = seed,
        n_cores = n_cores, rtol = rtol
      ))[["elapsed"]]
      record("random", k, seed, 0, 0, t, list(
        q = a$Q[[1]], n_iter = tgc_iterations(a), loglik = a$loglik
      ))
    }
    if ("svd" %in% strategies) {
      Sys.sleep(2)
      set.seed(seed)
      t_start <- system.time(st <- ns$radmix_svd_start(data,
        factors = factors, k = k, n_cores = n_cores
      ))[["elapsed"]]
      t_fit <- system.time(fit <- ns$radmix_fit(data,
        k = k,
        init_q = st$q, init_p = st$p, rtol = rtol, n_cores = n_cores
      ))[["elapsed"]]
      record("svd", k, seed, t_svd, t_start, t_fit, fit)
    }
  }
}
if ("svd_call" %in% strategies) {
  Sys.sleep(2)
  t <- system.time(a <- gt_radmix(x,
    k = ks, seed = 1, init = "svd",
    n_cores = n_cores, rtol = rtol
  ))[["elapsed"]]
  for (i in seq_along(ks)) {
    record(
      "svd_call", ks[i], 1, NA, NA, NA,
      list(
        q = a$Q[[i]], n_iter = tgc_iterations(list(log = a$log[i])),
        loglik = a$loglik[i]
      ),
      call_s = t
    )
  }
  cat("svd_call: one call for K =", paste(ks, collapse = ","), "took", t, "s\n")
}

# summary of all the runs of this dataset in the csv (also from earlier
# calls), per K, rtol and strategy: gaps to the best log-likelihood of any
# run at that K, and the largest difference in Q from the best run (after
# matching the clusters)
res <- read.csv(csv)
res <- res[res$dataset == dataset & res$n_cores == n_cores, ]
qs <- readRDS(q_file)
sink(file.path(out, paste0("svd_vs_random_summary_", dataset, ".txt")))
cat("dataset", dataset, "n_cores", n_cores, "\n")
for (kk in sort(unique(res$k))) {
  for (rt in sort(unique(res$rtol))) {
    r <- subset(res, k == kk & rtol == rt)
    best <- r[which.max(r$loglik), ]
    q_best <- qs[[paste(best$strategy, kk, best$seed, rt, sep = "_")]]
    r$gap <- best$loglik - r$loglik
    r$max_dq <- vapply(seq_len(nrow(r)), function(i) {
      q <- qs[[paste(r$strategy[i], kk, r$seed[i], rt, sep = "_")]]
      if (is.null(q) || is.null(q_best)) return(NA_real_)
      max(abs(tgc_match_clusters(q, q_best) - q_best))
    }, 0)
    s <- do.call(rbind, lapply(split(r, r$strategy), function(d) {
      data.frame(
        k = kk, rtol = rt, strategy = d$strategy[1], runs = nrow(d),
        svd_s = round(max(d$svd_s), 1),
        start_s = round(stats::median(d$start_s), 1),
        fit_s = round(stats::median(d$fit_s), 1),
        call_s = round(max(d$call_s), 1),
        iter = stats::median(d$iterations),
        at_best = sum(d$gap < 0.1), near_5 = sum(d$gap < 5),
        modes = length(unique(round(d$loglik))),
        worst_gap = signif(max(d$gap), 3),
        max_dq_from_best = signif(max(d$max_dq), 3)
      )
    }))
    print(s, row.names = FALSE)
    cat("\n")
  }
}
sink()
cat(readLines(file.path(out, paste0(
  "svd_vs_random_summary_", dataset, ".txt"
))), sep = "\n")

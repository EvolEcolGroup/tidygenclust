# Shared helpers of the algorithmic experiments on gt_radmix() (fewer passes
# over the genotypes): data loading, pass counting, an R replica of the
# fitting loop that records each iteration (for the cross-validation
# policies), and the warm-ups. Sourced by the scripts in this folder; they
# need tidygenclust loaded first (library(tidygenclust, lib.loc = ...)), and
# use its internals through `:::`.
#
# Units. A "pass" is one pass over all the genotypes of the data used (one
# EM step, one Q or P block update, one log-likelihood, one deviance). A pass
# over a subset of the loci counts as its fraction of a full pass ("full-pass
# equivalents"). Passes with a cross-validation fold hidden cost ~1.5x a
# plain pass in time (dev_notes/radmix_performance.md); they are counted as
# passes here, and the time ratio is noted where it matters.

tg <- asNamespace("tidygenclust")

# Paths -------------------------------------------------------------------

algo_paths <- function(scratch = Sys.getenv(
                         "ALGO_SCRATCH",
                         "/media/andrea/annex/tmp/claude/tidygenclust/algo"
                       )) {
  hgdp_dir <- Sys.getenv(
    "HGDP_DIR",
    "/media/andrea/annex/git/tidygenclust/vignettes/articles/data/hgdp"
  )
  list(
    scratch = scratch,
    data = file.path(scratch, "data"),
    out = file.path(scratch, "out"),
    bed = c(
      hgdp = file.path(hgdp_dir, "hgdp_bench.bed"),
      hgdp_full = file.path(hgdp_dir, "hgdp_full.bed"),
      hgdp_quick = file.path(hgdp_dir, "hgdp_bench_quick.bed"),
      eur = Sys.getenv(
        "EUR_BED",
        file.path(tools::R_user_dir("tidygenclust", "cache"), "EUR_subset.bed")
      )
    )
  )
}

# A gen_tibble of one of the benchmark datasets, created once and then
# reloaded from its .gt file (backing files in the scratch data folder).
# "<name>_imputed" is a separate copy with the missing genotypes imputed
# (gt_impute_simple(), the most common genotype), for the SVD start, which
# needs them; imputing modifies the backing file, so the copy keeps the
# plain one untouched for the jobs reading it.
load_gt <- function(name, paths = algo_paths()) {
  dir.create(paths$data, showWarnings = FALSE, recursive = TRUE)
  gt_file <- file.path(paths$data, paste0(name, ".gt"))
  if (file.exists(gt_file)) {
    return(tidypopgen::gt_load(gt_file))
  }
  imputed <- grepl("_imputed$", name)
  x <- tidypopgen::gen_tibble(paths$bed[[sub("_imputed$", "", name)]],
    backingfile = file.path(paths$data, name), quiet = TRUE
  )
  if (imputed) x <- tidypopgen::gt_impute_simple(x, method = "mode")
  tidypopgen::gt_save(x, quiet = TRUE)
  x
}

# Pass counting -------------------------------------------------------------

# Passes over the genotypes of a radmix_fit() run, from its log (verbose =
# FALSE): the EM steps, the initial log-likelihood, and per iteration the
# block updates (2 per free block), the log-likelihood of the quasi-Newton
# point (unless the step failed, QN = -Inf) and that of x2 (when the point
# was rejected, "basic" is computed)
fit_passes <- function(fit, n_em_iter = 5, em_ran = TRUE, fix_q = FALSE,
                       fix_p = FALSE) {
  it <- grep("^Iteration [0-9]+:", fit$log, value = TRUE)
  qn <- !grepl("QN = -Inf", it, fixed = TRUE)
  basic <- !grepl("basic = not computed", it, fixed = TRUE)
  n_upd <- 2 * (!fix_q) + 2 * (!fix_p)
  (if (em_ran) n_em_iter else 0) + 1 + length(it) * n_upd + sum(qn) +
    sum(basic)
}

# The random start of radmix_fit_cpp(), drawn in the same order from R's
# generator (q, rescaled so that each individual sums to one, then p), in the
# internal layout (K x individuals, K x loci). The sums are accumulated in
# the same order as in C++, so the start is bit-identical.
random_start <- function(k, n_indiv, n_loci) {
  q <- matrix(stats::runif(k * n_indiv), k, n_indiv)
  s <- q[1, ]
  if (k > 1) for (kk in 2:k) s <- s + q[kk, ]
  q <- sweep(q, 2, s, "/")
  p <- matrix(stats::runif(k * n_loci), k, n_loci)
  list(q = q, p = p)
}

# An R replica of the main loop of radmix_fit_cpp() (src/radmix_fit.cpp),
# from (q, p) in the internal layout, calling the same kernels in the same
# order, so that it gives bit-identical iterates (checked against
# radmix_fit() by check_replica.R). `on_iter(iter, q, p, ll, reldiff,
# passes)` is called at the start (iter = 0) and after each iteration, e.g.
# to score the held-out genotypes along the trajectory of a
# cross-validation refit (`passes` counts the passes of the fit so far).
# Returns the estimates (internal layout), the history and the passes.
fit_replica <- function(data, q, p, max_iter = 1000, rtol = 1e-8, n_qn = 3,
                        n_cores = 1, on_iter = NULL) {
  k <- nrow(q)
  ll_new <- tg$radmix_loglik_cpp(data, q, p, n_cores)
  passes <- 1
  history <- ll_new
  if (!is.null(on_iter)) on_iter(0, q, p, ll_new, NA_real_, passes)
  n_par <- length(q) + length(p)
  u_mat <- matrix(0, n_par, n_qn)
  v_mat <- matrix(0, n_par, n_qn)
  converged <- FALSE
  iter <- 0
  for (iter in seq_len(max_iter)) {
    ll_prev <- ll_new
    q1 <- tg$radmix_update_q_cpp(data, q, p, n_cores)
    p1 <- tg$radmix_update_p_cpp(data, q1, p, n_cores)
    q2 <- tg$radmix_update_q_cpp(data, q1, p1, n_cores)
    p2 <- tg$radmix_update_p_cpp(data, q2, p1, n_cores)
    passes <- passes + 4
    # the secant pair in column (iter - 1) %% n_qn, as SecantPairs::store()
    col <- (iter - 1) %% n_qn + 1
    u_mat[, col] <- c(q1 - q, p1 - p)
    v_mat[, col] <- c(q2 - q1, p2 - p1)
    n_used <- min(iter, n_qn)
    pt <- tg$radmix_qn_step_cpp(
      u_mat[, seq_len(n_used), drop = FALSE],
      v_mat[, seq_len(n_used), drop = FALSE], q, p, q1, p1, n_cores
    )
    ll_qn <- -Inf
    if (!is.null(pt)) {
      ll_qn <- tg$radmix_loglik_cpp(data, pt$q, pt$p, n_cores)
      passes <- passes + 1
    }
    if (is.finite(ll_qn) && ll_prev < ll_qn) {
      q <- pt$q
      p <- pt$p
      ll_new <- ll_qn
    } else {
      ll_new <- tg$radmix_loglik_cpp(data, q2, p2, n_cores)
      passes <- passes + 1
      q <- q2
      p <- p2
    }
    history <- c(history, ll_new)
    reldiff <- abs((ll_new - ll_prev) / ll_prev)
    if (!is.null(on_iter)) on_iter(iter, q, p, ll_new, reldiff, passes)
    if (!is.na(reldiff) && reldiff < rtol) {
      converged <- TRUE
      break
    }
  }
  list(
    q = q, p = p, loglik = ll_new, loglik_history = history, n_iter = iter,
    converged = converged, passes = passes
  )
}

# The data restricted to a subset of the loci (positions in data$cols): the
# kernels take the loci from `cols`, so this is a view, not a copy
subset_loci <- function(data, idx) {
  data$cols <- data$cols[idx]
  data$n_loci <- length(idx)
  data
}

# Warm-ups -------------------------------------------------------------------

# P for all the loci given Q (internal layout, K x individuals): one EM step
# from flat frequencies (which gives the Q-weighted allele frequencies,
# p_kj = sum_i q_ki g_ij / sum_i n_i q_ki) and `n_newton` Newton P updates.
# The EM step also returns an updated Q, which is discarded. Returns P (K x
# loci) and the number of full passes.
extend_p <- function(data, q, n_newton = 2, n_cores = 1) {
  k <- nrow(q)
  f <- matrix(0.5, k, data$n_loci)
  p <- tg$radmix_em_cpp(data, q, f, n_cores)$p
  for (i in seq_len(n_newton)) p <- tg$radmix_update_p_cpp(data, q, p, n_cores)
  list(p = p, passes = 1 + n_newton)
}

# Warm-up on a random subset of the loci: a random start and 5 EM steps on
# the subset, a fit to `rtol_sub` on the subset, P extended to all the loci
# given its Q (extend_p()). Returns (q, p) in the user layout (individuals x
# K, loci x K), ready for radmix_fit(init_q, init_p), and the passes in
# full-pass equivalents.
warmup_subset <- function(data, k, frac, rtol_sub = 1e-6, n_newton = 2,
                          n_cores = 1) {
  m <- data$n_loci
  idx <- sort(sample.int(m, max(1, round(frac * m))))
  sub <- subset_loci(data, idx)
  fit_s <- tg$radmix_fit(sub, k = k, rtol = rtol_sub, n_cores = n_cores)
  ext <- extend_p(data, t(fit_s$q), n_newton = n_newton, n_cores = n_cores)
  list(
    q = fit_s$q, p = t(ext$p),
    passes = fit_passes(fit_s) * length(idx) / m + ext$passes,
    sub_iter = fit_s$n_iter, sub_loglik = fit_s$loglik
  )
}

# admixer's mini-batch warm-up (Fitter::warmup in admixer's src/fit.hpp, as
# summarised in the speed-up plan): the loci are shuffled once and split
# into B = min(32, M / 500) batches; for each batch, a Newton P step on its
# loci, then a Newton Q step from its loci only (so Q is updated B times per
# epoch); after each epoch the full log-likelihood is computed, and B is
# halved whenever an epoch does not raise it (the estimates before that
# epoch are then kept), until B = 1, or `max_epochs`. With `rtol_halve` > 0,
# B is also halved when an epoch raises it by less than that relative
# amount (its estimates are then kept). From (q, p) in the
# internal layout; returns (q, p) in the user layout, the passes (an epoch
# is 2 passes of updates and 1 of log-likelihood) and the epochs.
warmup_minibatch <- function(data, q, p, max_batches = 32, min_batch = 500,
                             max_epochs = 100, rtol_halve = 0,
                             n_cores = 1) {
  m <- data$n_loci
  n_batches <- max(1, min(max_batches, floor(m / min_batch)))
  order_loci <- sample.int(m)
  ll <- tg$radmix_loglik_cpp(data, q, p, n_cores)
  passes <- 1
  epochs <- 0
  history <- ll
  while (n_batches > 1 && epochs < max_epochs) {
    epochs <- epochs + 1
    q_new <- q
    p_new <- p
    # batches of consecutive positions in the shuffled order, each sorted
    # (reading the loci in file order)
    batch <- cut(seq_len(m), n_batches, labels = FALSE)
    for (b in seq_len(n_batches)) {
      idx <- sort(order_loci[batch == b])
      sub <- subset_loci(data, idx)
      p_new[, idx] <- tg$radmix_update_p_cpp(sub, q_new, p_new[, idx,
        drop = FALSE
      ], n_cores)
      q_new <- tg$radmix_update_q_cpp(sub, q_new, p_new[, idx, drop = FALSE],
        n_cores)
    }
    ll_new <- tg$radmix_loglik_cpp(data, q_new, p_new, n_cores)
    passes <- passes + 3
    history <- c(history, ll_new)
    if (ll_new > ll) {
      gain <- (ll_new - ll) / abs(ll)
      q <- q_new
      p <- p_new
      ll <- ll_new
      if (gain < rtol_halve) n_batches <- n_batches %/% 2
    } else {
      n_batches <- n_batches %/% 2
    }
  }
  list(q = t(q), p = t(p), passes = passes, epochs = epochs, loglik = ll,
    history = history)
}

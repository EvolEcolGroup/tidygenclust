# Part A summary: the cross-validation error of each refit policy, read
# from the trajectories recorded by cv_trajectories.R.
#
# A policy stops each fold's refit at some iteration of the (deterministic)
# trajectory from the full-data estimates:
#   conv       the current policy: rtol 1e-8 (the end of the trajectory)
#   it<m>      at most m iterations (radmix_cv(max_iter = m)); it0 is no
#              refit at all (the full-data estimates score the fold)
#   rtol<r>    the first iteration whose relative change in log-likelihood
#              is below r (radmix_cv(rtol = r))
# Its cost is the passes of the refits up to that iteration (with the
# initial log-likelihood) plus one deviance pass per fold; passes with a
# fold hidden cost ~1.5x a plain pass in time.
#
# Usage: Rscript cv_policies.R <cv_dir> [<out_prefix>] [<admixture_fixture>]
# Writes <out_prefix>_errors.csv (one row per dataset, K, seed, policy) and
# prints the tables of the report. With the ADMIXTURE fixture of the tests
# (tests/testthat/fixtures/admixture_eur_reference.rds, CV errors of
# ADMIXTURE 1.3.0 on EUR_subset, K = 2-4, seed 1), the EUR errors of each
# policy are also compared with ADMIXTURE's (whose folds differ from ours).
args <- commandArgs(trailingOnly = TRUE)
cv_dir <- args[1]
out_prefix <- if (length(args) >= 2) args[2] else file.path(cv_dir, "policies")
fixture <- if (length(args) >= 3) args[3] else NULL
files <- list.files(cv_dir, pattern = "_k[0-9]+_s[0-9]+\\.rds$",
  full.names = TRUE)
policies <- c("conv", "it0", "it1", "it2", "it3", "it5", "it10",
  "rtol1e-7", "rtol1e-6", "rtol1e-5", "rtol1e-4")

# the iteration at which a policy stops one fold's trajectory
stop_iter <- function(tr, policy) {
  n_conv <- max(tr$iter)
  if (policy == "conv") return(n_conv)
  if (grepl("^it", policy)) {
    return(min(as.integer(sub("^it", "", policy)), n_conv))
  }
  r <- as.numeric(sub("^rtol", "", policy))
  hit <- tr$iter[!is.na(tr$reldiff) & tr$reldiff < r]
  if (length(hit) == 0) n_conv else min(hit)
}

rows <- list()
for (f in files) {
  res <- readRDS(f)
  tr_all <- res$traj
  for (policy in policies) {
    dev <- 0
    n <- 0
    passes <- 0
    iters <- integer(0)
    for (fold in sort(unique(tr_all$fold))) {
      tr <- tr_all[tr_all$fold == fold, ]
      it <- stop_iter(tr, policy)
      row <- tr[tr$iter == it, ]
      dev <- dev + row$dev
      n <- n + row$n
      passes <- passes + row$passes + 1
      iters <- c(iters, it)
    }
    rows[[length(rows) + 1]] <- data.frame(
      dataset = res$dataset, k = res$k, seed = res$seed, policy = policy,
      cv = dev / n, mean_iter = mean(iters), max_iter = max(iters),
      cv_passes = passes, fit_passes = res$fit$passes,
      fit_loglik = res$fit$loglik
    )
  }
}
tab <- do.call(rbind, rows)
conv <- tab[tab$policy == "conv", c("dataset", "k", "seed", "cv", "cv_passes")]
names(conv)[4:5] <- c("cv_conv", "passes_conv")
tab <- merge(tab, conv)
tab$delta <- tab$cv - tab$cv_conv
tab$rel_passes <- tab$cv_passes / tab$passes_conv
tab <- tab[order(tab$dataset, tab$seed, tab$k, match(tab$policy, policies)), ]
utils::write.csv(tab, paste0(out_prefix, "_errors.csv"), row.names = FALSE)

fmt_tab <- function(d) {
  if (requireNamespace("knitr", quietly = TRUE)) {
    print(knitr::kable(d, row.names = FALSE, digits = 6))
  } else {
    print(d, row.names = FALSE)
  }
  cat("\n")
}

for (ds in unique(tab$dataset)) {
  d <- tab[tab$dataset == ds, ]
  cat("\n## ", ds, "\n\n", sep = "")
  # CV error per K and policy, averaged over seeds
  cat("CV error per K (mean over seeds", paste(sort(unique(d$seed)),
    collapse = ","), ")\n\n")
  w <- tapply(d$cv, list(d$k, d$policy), mean)[, policies, drop = FALSE]
  fmt_tab(data.frame(k = rownames(w), round(w, 6), check.names = FALSE))
  cat("Difference from the converged policy (x 1e-4), mean over seeds\n\n")
  w <- tapply(d$delta * 1e4, list(d$k, d$policy), mean)[, policies,
    drop = FALSE]
  fmt_tab(data.frame(k = rownames(w), round(w, 3), check.names = FALSE))
  cat("K with the lowest CV error, per seed\n\n")
  best <- sapply(policies, function(pol) {
    sapply(sort(unique(d$seed)), function(s) {
      ds_ <- d[d$policy == pol & d$seed == s, ]
      ds_$k[which.min(ds_$cv)]
    })
  })
  best <- matrix(best, ncol = length(policies), dimnames = list(NULL, policies))
  fmt_tab(data.frame(seed = sort(unique(d$seed)), best, check.names = FALSE))
  cat("Cost: refit iterations per fold (mean over K, seeds and folds), CV",
    "passes relative to the converged policy and to the main fits (ratios of",
    "the totals over K and seeds)\n\n")
  cost <- data.frame(
    policy = policies,
    mean_iter = sapply(policies, function(p) mean(d$mean_iter[d$policy == p])),
    rel_passes = sapply(policies, function(p) {
      sum(d$cv_passes[d$policy == p]) / sum(d$passes_conv[d$policy == p])
    }),
    cv_over_fit = sapply(policies, function(p) {
      sum(d$cv_passes[d$policy == p]) / sum(d$fit_passes[d$policy == p])
    })
  )
  fmt_tab(cost)
}

if (!is.null(fixture) && "eur" %in% tab$dataset) {
  runs <- readRDS(fixture)$runs
  adm <- data.frame(
    k = vapply(runs, function(r) r$k, numeric(1)),
    admixture = vapply(runs, function(r) r$cv, numeric(1))
  )
  d <- tab[tab$dataset == "eur" & tab$k %in% adm$k, ]
  w <- tapply(d$cv, list(d$k, d$policy), mean)[, policies, drop = FALSE]
  cat("\n## EUR: CV error minus ADMIXTURE's (x 1e-4; ADMIXTURE prints 5",
    "decimals)\n\n")
  diff <- (w - adm$admixture[match(rownames(w), adm$k)]) * 1e4
  fmt_tab(data.frame(k = rownames(w), admixture = adm$admixture[match(
    rownames(w), adm$k
  )], round(diff, 2), check.names = FALSE))
}

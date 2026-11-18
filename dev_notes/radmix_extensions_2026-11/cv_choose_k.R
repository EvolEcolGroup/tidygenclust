# Does the cross-validation error of genotype likelihoods choose the K of
# the simulation? Simulated data (K = 3, admixed), genotypes and their
# likelihoods at depths 1, 2 and 4; CV errors for K = 1 to 5, 3 seeds.
# Usage: Rscript cv_choose_k.R <tidygenclust lib> <out csv>
args <- commandArgs(trailingOnly = TRUE)
library(tidygenclust, lib.loc = args[1])
ns <- asNamespace("tidygenclust")
for (nm in ls(ns, pattern = "^radmix_")) assign(nm, get(nm, ns))
source("/home/andrea/git/tidygenclust/tests/testthat/helper-radmix-data.R")
n_cores <- as.integer(Sys.getenv("N_CORES", "2"))
rows <- list()
for (seed in 1:3) {
  sim <- radmix_simulate_admixed(n = 200, m = 20000, k = 3, missing = 0.01,
                                 seed = seed)
  gt <- radmix_test_gt(sim$genotypes)
  datasets <- list(genotypes = gt)
  for (depth in c(1, 2, 4)) {
    set.seed(100 + seed)
    ll <- radmix_simulate_gl(sim$genotypes, depth = depth)
    x <- gt
    class(x) <- c(paste0("cvk_gl_", depth), class(x))
    registerS3method("gt_radmix_data", paste0("cvk_gl_", depth),
                     local({ d <- radmix_test_gl(ll); function(x, ...) d }),
                     envir = ns)
    datasets[[paste0("depth_", depth)]] <- x
  }
  for (nm in names(datasets)) {
    t <- system.time(adm <- gt_radmix(datasets[[nm]], k = 1:5, seed = seed,
                                      crossval = TRUE, n_cores = n_cores))
    rows[[length(rows) + 1]] <- data.frame(
      seed = seed, data = nm, k = adm$k, cv = adm$cv, loglik = adm$loglik,
      best = adm$k == adm$k[which.min(adm$cv)], time_s = t[["elapsed"]]
    )
    message(seed, " ", nm, ": best K = ", adm$k[which.min(adm$cv)],
            "; cv = ", paste(signif(adm$cv, 5), collapse = " "))
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, args[2], row.names = FALSE)

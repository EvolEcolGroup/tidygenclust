# Compare gt_radmix() with the ADMIXTURE binary on EUR_subset, K = 4
# Requires the EUR_subset data: run data-raw/radmix/download_eur_subset.R first.
# Run from the root of the package, with tidygenclust installed (timings
# from devtools::load_all() are misleading, as it compiles without
# optimisation).
library(tidygenclust)
x <- tidypopgen::gt_load("data-raw/radmix/data/eur.gt")
t_radmix <- system.time(fit_radmix <- gt_radmix(x, k = 4, seed = 1, n_cores = 4))
t_adm <- system.time(
  adm <- tidypopgen::gt_admixture(x,
    k = 4, seed = 1, n_cores = 4,
    conda_env = "ctidygenclust"
  )
)
cat(
  "gt_radmix: loglik", format(fit_radmix$loglik, digits = 12),
  " time", t_radmix[3], "s\n"
)
cat(
  "ADMIXTURE binary : loglik", format(adm$loglik, digits = 12),
  " time", t_adm[3], "s\n"
)
q1 <- unclass(fit_radmix$Q[[1]])
q2 <- as.matrix(adm$Q[[1]])
cc <- stats::cor(q1, q2)
perm <- apply(cc, 1, which.max)
cat("matched columns:", perm, "\n")
cat("max |Q difference|:", max(abs(q1 - q2[, perm])), "\n")
cat("mean |Q difference|:", mean(abs(q1 - q2[, perm])), "\n")

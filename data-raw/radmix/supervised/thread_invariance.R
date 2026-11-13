# Thread invariance of supervised gt_radmix() on HGDP pruned, K = 5: the
# fits with n_cores = 1 and n_cores = 8 must give identical() Q, P and
# log-likelihood, for both sets of references of hgdp_supervised.R.
#
# Usage, from the root of the package:
#   Rscript data-raw/radmix/supervised/thread_invariance.R <library> \
#     <scratch dir> [n_cores]
# (<scratch dir> holds the gen_tibble made by hgdp_supervised.R)
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
scratch <- args[2]
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 8L
library(tidygenclust, lib.loc = lib)

hgdp_dir <- Sys.getenv("TGC_HGDP_DIR", "vignettes/articles/data/hgdp")
hgdp <- tidypopgen::gt_load(file.path(scratch, "hgdp_bench.gt"))
lab <- utils::read.delim(file.path(hgdp_dir, "hgdp650_id_pop.txt"))
lab <- lab[match(hgdp$id, lab$Id), ]
cases <- list(
  pops = ifelse(lab$population %in% c("Yoruba", "French", "Han", "Pima", "Papuan"),
    lab$population, NA
  ),
  regions = ifelse(lab$Pop7Groups %in% c("Africa", "Europe", "Est_Asia", "America", "Oceania"),
    lab$Pop7Groups, NA
  )
)
res <- do.call(rbind, lapply(names(cases), function(case) {
  fits <- lapply(c(1L, n_cores), function(nc) {
    gt_radmix(hgdp, k = 5, supervised = cases[[case]], seed = 1, n_cores = nc)
  })
  data.frame(
    case = case, n_cores = n_cores,
    q_identical = identical(fits[[1]]$Q, fits[[2]]$Q),
    p_identical = identical(fits[[1]]$P, fits[[2]]$P),
    loglik_identical = identical(fits[[1]]$loglik, fits[[2]]$loglik),
    loglik = fits[[1]]$loglik
  )
}))
print(res, digits = 15)
utils::write.csv(res, file.path(scratch, "thread_invariance.csv"), row.names = FALSE)
stopifnot(all(res$q_identical), all(res$p_identical), all(res$loglik_identical))

# Cost of the passes over the genotypes with and without the references
# (HGDP pruned, K = 5, references as in hgdp_supervised.R): the median time
# of update_q() and of an EM step without a mask and with the masks of the
# references, next to update_p() and loglik(), which a supervised fit does
# not change. One thread.
#
# Usage, from the root of the package:
#   Rscript data-raw/radmix/supervised/time_kernels.R <library> \
#     <scratch dir> [n_reps]
# (<scratch dir> holds the gen_tibble made by hgdp_supervised.R)
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
scratch <- args[2]
n_reps <- if (length(args) >= 3) as.integer(args[3]) else 3L
library(tidygenclust, lib.loc = lib)

hgdp_dir <- "/media/andrea/annex/git/tidygenclust/vignettes/articles/data/hgdp"
hgdp <- tidypopgen::gt_load(file.path(scratch, "hgdp_bench.gt"))
lab <- utils::read.delim(file.path(hgdp_dir, "hgdp650_id_pop.txt"))
lab <- lab[match(hgdp$id, lab$Id), ]
ref_regions <- c("Africa", "Europe", "Est_Asia", "America", "Oceania")
ref_pops <- c("Yoruba", "French", "Han", "Pima", "Papuan")
masks <- list(
  none = NULL,
  regions = lab$Pop7Groups %in% ref_regions,
  pops = lab$population %in% ref_pops
)
dat <- tidygenclust:::radmix_prepare(hgdp)
k <- 5
set.seed(1)
q <- matrix(stats::runif(k * dat$n_indiv), k)
q <- sweep(q, 2, colSums(q), "/")
p <- matrix(stats::runif(k * dat$n_loci, 0.05, 0.95), k)

time_it <- function(f) {
  stats::median(vapply(seq_len(n_reps), function(r) {
    system.time(f())[["elapsed"]]
  }, numeric(1)))
}
res <- NULL
for (m in names(masks)) {
  fixed <- masks[[m]]
  t_q <- time_it(function() {
    tidygenclust:::radmix_update_q_cpp(dat, q, p, 1L, fixed)
  })
  t_em <- time_it(function() {
    tidygenclust:::radmix_em_cpp(dat, q, p, 1L, fixed)
  })
  res <- rbind(res, data.frame(
    mask = m, n_fixed = sum(fixed), update_q = t_q, em = t_em,
    load = readLines("/proc/loadavg")
  ))
}
t_p <- time_it(function() tidygenclust:::radmix_update_p_cpp(dat, q, p, 1L))
t_ll <- time_it(function() tidygenclust:::radmix_loglik_cpp(dat, q, p, 1L))
res$update_p <- t_p
res$loglik <- t_ll
print(res)
utils::write.csv(res, file.path(scratch, "time_kernels.csv"),
  row.names = FALSE)

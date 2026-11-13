# Starts of supervised fits: ADMIXTURE's random start + 5 EM steps against
# fastmixture's start (P from the allele frequencies of the reference
# populations), on EUR_subset (K = 3, references as in compare_eur.R) or
# HGDP pruned (K = 5, references as in hgdp_supervised.R, case "regions").
#
# Usage, from the root of the package:
#   Rscript data-raw/radmix/supervised/compare_starts.R <library> \
#     <scratch dir> <data: eur|hgdp> [n_cores] [seeds] [starts]
# (seeds: an R expression, e.g. 1:3; starts: a comma-separated subset of
# random,freq_em,freq)
# Starts compared, for each seed:
# - random: random Q and P, 5 EM steps (the default, as ADMIXTURE);
# - freq_em: P from the reference frequencies (init_p), Q of the others
#   from the 5 EM steps of a random start (what init_p alone gives);
# - freq: P from the reference frequencies, Q of the others random, no EM
#   (fastmixture's start).
# Prints and writes <scratch dir>/compare_starts_<data>.csv: log-likelihood,
# iterations and time of each fit.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
scratch <- args[2]
data_name <- args[3]
n_cores <- if (length(args) >= 4) as.integer(args[4]) else 1L
seeds <- if (length(args) >= 5) eval(parse(text = args[5])) else 1:3
starts <- if (length(args) >= 6) {
  strsplit(args[6], ",")[[1]]
} else {
  c("random", "freq_em", "freq")
}
library(tidygenclust, lib.loc = lib)

if (data_name == "eur") {
  x <- tidypopgen::gen_tibble(Sys.getenv("TGC_EUR_BED", file.path(tools::R_user_dir("tidygenclust", "cache"), "EUR_subset.bed")),
    backingfile = tempfile("eur_", tmpdir = scratch), quiet = TRUE
  )
  fixture <- readRDS("tests/testthat/fixtures/admixture_eur_reference.rds")
  q3 <- fixture$runs[[2]]$q_millionths / 1e6
  labels <- c("mid", "zeta", "alpha")[apply(q3, 1, which.max)]
  labels[apply(q3, 1, max) <= 0.9] <- NA
  k <- 3
} else {
  hgdp_dir <- Sys.getenv("TGC_HGDP_DIR", "vignettes/articles/data/hgdp")
  # (the gen_tibble made by hgdp_supervised.R, if it exists)
  gt_file <- file.path(scratch, "hgdp_bench.gt")
  if (file.exists(gt_file)) {
    x <- tidypopgen::gt_load(gt_file)
  } else {
    x <- tidypopgen::gen_tibble(file.path(hgdp_dir, "hgdp_bench.bed"),
      backingfile = file.path(scratch, "hgdp_bench"), quiet = TRUE
    )
  }
  lab <- utils::read.delim(file.path(hgdp_dir, "hgdp650_id_pop.txt"))
  lab <- lab[match(x$id, lab$Id), ]
  ref_regions <- c("Africa", "Europe", "Est_Asia", "America", "Oceania")
  labels <- ifelse(lab$Pop7Groups %in% ref_regions, lab$Pop7Groups, NA)
  k <- 5
}

# the frequency of the alternate allele in each reference population, in
# the order of the columns (first appearance), within [1e-5, 1 - 1e-5]
pops <- unique(labels[!is.na(labels)])
p_ref <- sapply(pops, function(pop) {
  tidypopgen::loci_alt_freq(x[which(labels == pop), ])
})
p_ref <- pmin(pmax(p_ref, 1e-5), 1 - 1e-5)

n_iter_of <- function(fit) {
  as.integer(sub(".* ([0-9]+) iterations.*", "\\1",
    grep("iterations\\.$", fit$log[[1]], value = TRUE)))
}
res <- NULL
for (seed in seeds) {
  for (start in starts) {
    init_p <- if (start == "random") NULL else p_ref
    n_em <- if (start == "freq") 0 else 5
    t_fit <- system.time(fit <- gt_radmix(x,
      k = k, supervised = labels, init_p = init_p, n_em_iter = n_em,
      seed = seed, n_cores = n_cores
    ))[["elapsed"]]
    res <- rbind(res, data.frame(
      data = data_name, start = start, seed = seed, loglik = fit$loglik,
      iterations = n_iter_of(fit), time_s = t_fit,
      load = readLines("/proc/loadavg")
    ))
    print(utils::tail(res, 1), digits = 12)
  }
}
out_csv <- file.path(scratch, paste0("compare_starts_", data_name, ".csv"))
utils::write.table(res, out_csv,
  sep = ",", row.names = FALSE,
  col.names = !file.exists(out_csv), append = file.exists(out_csv)
)

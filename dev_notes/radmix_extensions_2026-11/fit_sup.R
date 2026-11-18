# One supervised (or unsupervised) gt_radmix() fit of the pruned HGDP
# (hgdp_bench: 1043 individuals, 130,943 SNPs), K = 5, timed; one fit per
# call, so that the builds can be alternated in separate processes.
#
# References (labels from hgdp650_id_pop.txt), as in the kit of the tag
# archive/supervised (data-raw/radmix/supervised/hgdp_supervised.R):
# - "pops": the individuals of five well separated populations (Yoruba,
#   French, Han, Pima, Papuan: 139 references, 904 estimated);
# - "regions": the individuals of five Pop7Groups regions (Africa, Europe,
#   Est_Asia, America, Oceania: 660 references, 383 estimated);
# - "unsupervised": no references.
#
# Usage:
#   Rscript fit_sup.R <library> <case> <n_cores> <seed> <out csv> <rds dir>
# Environment: TGC_HGDP_GT (the .gt file of hgdp_bench), TGC_HGDP_DIR (the
# folder with hgdp650_id_pop.txt).
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
case <- args[2]
n_cores <- as.integer(args[3])
seed <- as.integer(args[4])
out_csv <- args[5]
rds_dir <- args[6]
library(tidygenclust, lib.loc = lib)

gt_file <- Sys.getenv(
  "TGC_HGDP_GT", "/home/andrea/tgc_newspeed_work/sup/hgdp_bench.gt"
)
hgdp_dir <- Sys.getenv(
  "TGC_HGDP_DIR", "/home/andrea/git/tidygenclust/vignettes/articles/data/hgdp"
)
hgdp <- tidypopgen::gt_load(gt_file)
lab <- utils::read.delim(file.path(hgdp_dir, "hgdp650_id_pop.txt"))
lab <- lab[match(hgdp$id, lab$Id), ]
labels <- switch(case,
  pops = {
    ref_pops <- c("Yoruba", "French", "Han", "Pima", "Papuan")
    ifelse(lab$population %in% ref_pops, lab$population, NA)
  },
  regions = {
    ref_regions <- c("Africa", "Europe", "Est_Asia", "America", "Oceania")
    ifelse(lab$Pop7Groups %in% ref_regions, lab$Pop7Groups, NA)
  },
  unsupervised = NULL
)

load_avg <- readLines("/proc/loadavg")
t_fit <- system.time(fit <- gt_radmix(hgdp,
  k = 5, supervised = labels, seed = seed, n_cores = n_cores
))[["elapsed"]]
log_lines <- fit$log[[1]]
n_iter <- as.integer(sub(
  ".* (in|after) ([0-9]+) iterations.*", "\\2",
  grep("( in [0-9]+ iterations| after [0-9]+ iterations)", log_lines,
    value = TRUE
  )
))
# rejected quasi-Newton steps (the basic point was computed)
n_rejected <- sum(grepl("^Iteration .*basic = -?[0-9]", log_lines))
res <- data.frame(
  build = basename(lib), case = case, seed = seed, n_cores = n_cores,
  loglik = fit$loglik, iterations = n_iter, rejected = n_rejected,
  time_s = t_fit, load = load_avg, date = format(Sys.time())
)
print(res, digits = 15)
utils::write.table(res, out_csv,
  sep = ",", row.names = FALSE,
  col.names = !file.exists(out_csv), append = file.exists(out_csv)
)
saveRDS(
  list(Q = unclass(fit$Q[[1]]), P = fit$P[[1]], loglik = fit$loglik,
       log = log_lines),
  file.path(rds_dir, paste0(case, "_", n_cores, "_", seed, "_",
                            basename(lib), ".rds"))
)

# Supervised gt_radmix() on HGDP pruned (1043 individuals, 130,943 SNPs),
# K = 5: timings and iterations against the same data unsupervised, and the
# estimates against ADMIXTURE --supervised (if its output is available).
#
# Two sets of references (labels from hgdp650_id_pop.txt):
# - "pops": all the individuals of five well separated populations (Yoruba,
#   French, Han, Pima, Papuan: 139 references), the other 904 estimated;
# - "regions": all the individuals of five of the Pop7Groups regions
#   (Africa, Europe, Est_Asia, America, Oceania: 660 references), the 383
#   individuals of Central_South_Asia and Middle_Est estimated (references
#   are most of the individuals).
#
# Usage (from the root of the package; one fit per call, so that timings of
# different builds can be alternated):
#   Rscript data-raw/radmix/supervised/hgdp_supervised.R <library> \
#     <scratch dir> <case> <start> [n_cores] [seed]
# case: "pops", "regions" or "unsupervised"; start: "random" (random start +
# 5 EM steps, as ADMIXTURE) or "freq" (P started from the allele frequencies
# of the reference populations, as fastmixture; Q of the others random).
# Writes the .pop files for ADMIXTURE in <scratch dir>/admix, and
# appends a line to <scratch dir>/hgdp_supervised.csv.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
scratch <- args[2]
case <- args[3]
start <- args[4]
n_cores <- if (length(args) >= 5) as.integer(args[5]) else 1L
seed <- if (length(args) >= 6) as.integer(args[6]) else 1L
library(tidygenclust, lib.loc = lib)

hgdp_dir <- Sys.getenv("TGC_HGDP_DIR", "vignettes/articles/data/hgdp")
bed <- file.path(hgdp_dir, "hgdp_bench.bed")
# one gen_tibble (backing file) per scratch dir, reused across calls
# (gen_tibble() saves it as hgdp_bench.gt)
gt_file <- file.path(scratch, "hgdp_bench.gt")
if (!file.exists(gt_file)) {
  hgdp <- tidypopgen::gen_tibble(bed,
    backingfile = file.path(scratch, "hgdp_bench"), quiet = TRUE
  )
} else {
  hgdp <- tidypopgen::gt_load(gt_file)
}
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

# the .pop files for ADMIXTURE (written once)
adm_dir <- file.path(scratch, "admix")
dir.create(adm_dir, showWarnings = FALSE)
if (!is.null(labels)) {
  pop_file <- file.path(adm_dir, paste0("hgdp_", case, ".pop"))
  if (!file.exists(pop_file)) {
    writeLines(ifelse(is.na(labels), "-", labels), pop_file)
  }
}

# the start from the allele frequencies of the reference populations (the
# frequency of the alternate allele among their observed genotypes, as in
# fastmixture's initP())
init_p <- NULL
if (start == "freq") {
  pops <- unique(labels[!is.na(labels)])
  init_p <- sapply(pops, function(pop) {
    sub <- hgdp[which(labels == pop), ]
    tidypopgen::loci_alt_freq(sub)
  })
  init_p <- pmin(pmax(init_p, 1e-5), 1 - 1e-5)
}

load_avg <- readLines("/proc/loadavg")
t_fit <- system.time(fit <- gt_radmix(hgdp,
  k = 5, supervised = labels, init_p = init_p, seed = seed,
  n_cores = n_cores
))[["elapsed"]]
n_iter <- as.integer(sub(".* in ([0-9]+) iterations.*", "\\1",
  grep("( in [0-9]+ iterations| after [0-9]+ iterations)", fit$log[[1]],
    value = TRUE)))
res <- data.frame(
  build = basename(lib), case = case, start = start, seed = seed,
  n_cores = n_cores, loglik = fit$loglik, iterations = n_iter,
  time_s = t_fit, load = load_avg
)
print(res, digits = 12)
out_csv <- file.path(scratch, "hgdp_supervised.csv")
utils::write.table(res, out_csv,
  sep = ",", row.names = FALSE,
  col.names = !file.exists(out_csv), append = file.exists(out_csv)
)
saveRDS(fit, file.path(scratch, paste0(
  "hgdp_", case, "_", start, "_", seed, "_", basename(lib), ".rds"
)))

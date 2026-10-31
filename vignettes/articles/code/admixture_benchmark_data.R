# Preparation of the data for the benchmark of gt_radmix() against
# ADMIXTURE and fastmixture on the HGDP. The results are shown in the article
# admixture_benchmark.Rmd (in the parent directory of this script), which
# also shows the code of the chunks of this script (the sections that start
# with a "## ---- label" line, read with knitr::read_chunk()).
#
# This script downloads the HGDP data, filters and LD prunes them as
# recommended in the manual of ADMIXTURE, and writes the PLINK bed file that
# all the tools are run on (data/hgdp/hgdp_bench.bed), as well as a small
# subset of it (data/hgdp/hgdp_bench_quick.bed), to check the benchmark
# scripts in less time. Run it once, from the directory of the article
# (vignettes/articles):
#
#   cd vignettes/articles
#   Rscript code/admixture_benchmark_data.R
#
# The files in data/hgdp are not part of the repository (see .gitignore).

## ---- data-setup ----
hgdp_dir <- "data/hgdp"
dir.create(hgdp_dir, recursive = TRUE, showWarnings = FALSE)
# the backing files of the gen_tibbles used to prepare the data (removed at
# the end); gen_tibble() would rename them if they already existed
gt_dir <- file.path(hgdp_dir, "gt_prepare")
unlink(gt_dir, recursive = TRUE)
dir.create(gt_dir)

## ---- download ----
# the HGDP data (1043 individuals typed at about 650k loci) and the
# population of each individual, as used in the benchmark of tidypopgen.
# Each file is downloaded under a temporary name, so that an interrupted
# download is not mistaken for a complete file when the script is run again
zenodo_url <- "https://zenodo.org/records/15582364/files/"
hgdp_files <- c(
  "hgdp650.qc.hg19.bed", "hgdp650.qc.hg19.bim", "hgdp650.qc.hg19.fam",
  "hgdp650_id_pop.txt"
)
options(timeout = max(3600, getOption("timeout")))
for (file in hgdp_files) {
  path <- file.path(hgdp_dir, file)
  if (!file.exists(path)) {
    download.file(paste0(zenodo_url, file), paste0(path, ".part"), mode = "wb")
    file.rename(paste0(path, ".part"), path)
  }
}

## ---- filter ----
library(tidypopgen)
hgdp <- gen_tibble(
  file.path(hgdp_dir, "hgdp650.qc.hg19.bed"),
  backingfile = file.path(gt_dir, "hgdp_raw"),
  quiet = TRUE
)
# loci with less than 5% of missing genotypes, as recommended in the manual
# of ADMIXTURE, and polymorphic (fastmixture stops on monomorphic loci)
hgdp <- hgdp |>
  select_loci_if(loci_missingness(genotypes) < 0.05) |>
  select_loci_if(loci_maf(genotypes) > 0)

## ---- ld-clumping ----
# LD clumping with the thresholds of the manual of ADMIXTURE (plink
# --indep-pairwise 50 10 0.1). Clumping needs complete genotypes, so the
# missing ones are imputed first (on a backing file with only the loci kept
# above). loci_ld_clump() uses the imputed genotypes, but the gen_tibble
# keeps using the raw ones, which are those written to the bed file below:
# missing genotypes remain missing for all the tools
hgdp <- gt_update_backingfile(
  hgdp,
  backingfile = file.path(gt_dir, "hgdp_filtered"),
  quiet = TRUE
)
set.seed(1)
hgdp <- gt_impute_simple(hgdp, method = "random")
hgdp <- hgdp |>
  select_loci_if(
    loci_ld_clump(genotypes, thr_r2 = 0.1, size = 50, use_positions = FALSE)
  )
gt_as_plink(hgdp, file = file.path(hgdp_dir, "hgdp_bench.bed"))

## ---- quick-data ----
# a small subset to check the benchmark scripts (on which ADMIXTURE takes
# about 30 s on one core): every 5th individual, and 25,000 of the loci that
# are still polymorphic among them
hgdp_quick <- hgdp[seq(1, nrow(hgdp), by = 5), ] |>
  select_loci_if(loci_maf(genotypes) > 0)
n_loci <- count_loci(hgdp_quick)
set.seed(2)
quick_loci <- seq_len(n_loci) %in% sample(n_loci, min(25000, n_loci))
hgdp_quick <- select_loci_if(hgdp_quick, quick_loci)
gt_as_plink(hgdp_quick, file = file.path(hgdp_dir, "hgdp_bench_quick.bed"))

## ---- data-cleanup ----
message(
  "HGDP for the benchmark: ", nrow(hgdp), " individuals and ",
  count_loci(hgdp), " loci; quick subset: ", nrow(hgdp_quick),
  " individuals and ", count_loci(hgdp_quick), " loci"
)
unlink(gt_dir, recursive = TRUE)

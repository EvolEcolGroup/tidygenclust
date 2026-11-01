# Gets the datasets of the verification (see tgc_data() in common.R), and
# checks them against the files used on the laptop (md5 sums):
# - eur: EUR_subset (379 x 54,051), downloaded from SnpArrays.jl into the
#   cache of the tests;
# - hgdp_quick and hgdp: made by the data script of the benchmark article
#   (vignettes/articles/code/admixture_benchmark_data.R: downloads the HGDP
#   from Zenodo, ~190 MB, filters and LD prunes it);
# - hgdp_full: the HGDP filtered as for the benchmark (missingness < 5%,
#   polymorphic) but not LD pruned (1043 x 642,794), for the scaling.
#
#   Rscript data-raw/radmix/perf_verification/scripts/prepare_data.R
source("data-raw/radmix/perf_verification/scripts/common.R")
suppressMessages(library(tidypopgen))
hgdp_dir <- file.path(tgc_repo, "vignettes", "articles", "data", "hgdp")

# EUR_subset (as in tests/testthat/helper-radmix-admixture.R)
eur <- Sys.getenv("TGC_EUR_BED", file.path(
  tools::R_user_dir("tidygenclust", which = "cache"), "EUR_subset.bed"
))
if (!file.exists(eur)) {
  dir.create(dirname(eur), recursive = TRUE, showWarnings = FALSE)
  base_url <- "https://github.com/OpenMendel/SnpArrays.jl/raw/master/data/"
  for (ext in c("bed", "bim", "fam")) {
    utils::download.file(paste0(base_url, "EUR_subset.", ext),
      sub("bed$", ext, eur),
      mode = "wb"
    )
  }
}

# HGDP of the benchmark article (and its quick subset)
if (!file.exists(file.path(hgdp_dir, "hgdp_bench.bed")) ||
  !file.exists(file.path(hgdp_dir, "hgdp_bench_quick.bed"))) {
  status <- withr::with_dir(
    file.path(tgc_repo, "vignettes", "articles"),
    system2("Rscript", "code/admixture_benchmark_data.R")
  )
  stopifnot(status == 0)
}

# HGDP without LD pruning
if (!file.exists(file.path(hgdp_dir, "hgdp_full.bed"))) {
  hgdp <- gen_tibble(file.path(hgdp_dir, "hgdp650.qc.hg19.bed"),
    backingfile = tempfile("hgdp_raw", tmpdir = tgc_tmp), quiet = TRUE
  )
  hgdp <- hgdp |>
    select_loci_if(loci_missingness(genotypes) < 0.05) |>
    select_loci_if(loci_maf(genotypes) > 0)
  gt_as_plink(hgdp, file = file.path(hgdp_dir, "hgdp_full.bed"))
}

# what each dataset holds, and whether it is the file used on the laptop
laptop_md5 <- c(
  eur = "f5551ba6bee452dd9e7840972f59fd9c",
  hgdp_quick = "cd6dbc76c87e613bba840d3c3944dc88",
  hgdp = "2fe120a08c91e6eacc3140e3ea6b7472",
  hgdp_full = "9ba6eb536d1a6906b0e1fb43b13e6131"
)
for (d in names(laptop_md5)) {
  bed <- tgc_data(d)
  n_indiv <- length(readLines(sub("bed$", "fam", bed)))
  n_loci <- length(readLines(sub("bed$", "bim", bed)))
  md5 <- unname(tools::md5sum(bed))
  cat(sprintf(
    "%-10s %5d individuals x %7d loci  md5 %s  %s\n", d, n_indiv, n_loci,
    md5, if (is.na(laptop_md5[[d]])) "" else if (md5 == laptop_md5[[d]]) {
      "(same file as on the laptop)"
    } else {
      "(DIFFERENT from the laptop)"
    }
  ))
}

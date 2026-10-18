# Download the EUR_subset PLINK files used in the OpenADMIXTURE.jl tests
# (379 individuals x 54,051 SNPs, from the SnpArrays.jl repository) and convert
# them into a gen_tibble. The other scripts in data-raw/openadmixture/ rely on the output.
#
# Run from the root of the package:
#   Rscript data-raw/openadmixture/download_eur_subset.R

data_dir <- file.path("data-raw", "openadmixture", "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

base_url <- "https://github.com/OpenMendel/SnpArrays.jl/raw/master/data/"
for (ext in c("bed", "bim", "fam")) {
  dest <- file.path(data_dir, paste0("EUR_subset.", ext))
  if (!file.exists(dest)) {
    utils::download.file(paste0(base_url, "EUR_subset.", ext), dest,
      mode = "wb"
    )
  }
}

# create the gen_tibble (and its backing file), overwriting any old copy
unlink(file.path(data_dir, c("eur.bk", "eur.rds", "eur.gt")))
eur <- tidypopgen::gen_tibble(
  file.path(data_dir, "EUR_subset.bed"),
  backingfile = file.path(data_dir, "eur"),
  quiet = TRUE
)
tidypopgen::gt_save(eur, quiet = TRUE)
message(
  "Saved a gen_tibble with ", nrow(eur), " individuals and ",
  nrow(tidypopgen::show_loci(eur)), " loci in ", data_dir
)

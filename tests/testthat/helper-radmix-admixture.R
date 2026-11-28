# Helpers for the test that validates gt_radmix() against stored
# ADMIXTURE results (test_radmix_vs_admixture.R): downloading the
# EUR_subset data, and matching clusters between two Q matrices. When the data
# cannot be downloaded, the test is skipped.

# Path to the EUR_subset bed file of SnpArrays.jl (379 individuals, 54,051
# SNPs), downloaded from the SnpArrays.jl repository the
# first time it is needed. Files are cached in RADMIX_TEST_DATA if set,
# or in the user cache directory of the package. Returns NULL if the files
# cannot be downloaded.
radmix_eur_subset_bed <- function() {
  data_dir <- Sys.getenv(
    "RADMIX_TEST_DATA",
    tools::R_user_dir("tidygenclust", which = "cache")
  )
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  base_url <- "https://github.com/OpenMendel/SnpArrays.jl/raw/master/data/"
  for (ext in c("bed", "bim", "fam")) {
    dest <- file.path(data_dir, paste0("EUR_subset.", ext))
    if (!file.exists(dest)) {
      ok <- tryCatch(
        utils::download.file(paste0(base_url, "EUR_subset.", ext), dest,
          mode = "wb", quiet = TRUE
        ) == 0,
        error = function(e) FALSE,
        warning = function(w) FALSE
      )
      if (!ok) {
        unlink(dest)
        return(NULL)
      }
    }
  }
  file.path(data_dir, "EUR_subset.bed")
}

# Helpers for the test that validates gt_openadmixture() against stored
# ADMIXTURE results (test_openadmixture_vs_admixture.R): downloading the
# EUR_subset data, and matching clusters between two Q matrices. When the data
# cannot be downloaded, the test is skipped.

# Path to the EUR_subset bed file used in the OpenADMIXTURE.jl tests (379
# individuals, 54,051 SNPs), downloaded from the SnpArrays.jl repository the
# first time it is needed. Files are cached in OPENADMIXTURE_TEST_DATA if set,
# or in the user cache directory of the package. Returns NULL if the files
# cannot be downloaded.
openadmixture_eur_subset_bed <- function() {
  data_dir <- Sys.getenv(
    "OPENADMIXTURE_TEST_DATA",
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

# Reorder the columns of q to best match q_ref (the labels of the clusters
# are arbitrary); returns the permutation
openadmixture_match_clusters <- function(q, q_ref) {
  k <- ncol(q)
  # all the permutations of 1:n, one per row, built recursively: each value
  # i in turn, followed by every permutation of the other values
  perms <- function(n) {
    if (n == 1) {
      return(matrix(1))
    }
    sub <- perms(n - 1)
    do.call(rbind, lapply(seq_len(n), function(i) {
      cbind(i, matrix(setdiff(seq_len(n), i)[sub], nrow = nrow(sub)))
    }))
  }
  all_perms <- perms(k)
  diffs <- apply(all_perms, 1, function(x) sum(abs(q[, x] - q_ref)))
  all_perms[which.min(diffs), ]
}

# Regenerate the ADMIXTURE results used by tests/testthat/test_radmix_vs_admixture.R
#
# ADMIXTURE (version 1.3.0) is no longer developed and the EUR_subset data do
# not change, so rather than running ADMIXTURE every time the tests run (it
# takes several minutes), we run it once and store its estimates in
# tests/testthat/fixtures/admixture_eur_reference.rds. Rerun this script only if the
# comparison needs to change (e.g. different values of K or seeds).
#
# Requirements: ADMIXTURE, on the PATH or installed with
# Rscript data-raw/radmix/setup_admixture.R, and the EUR_subset data, downloaded
# with Rscript data-raw/radmix/download_eur_subset.R.
#
# Run from the root of the package (this takes about 30-40 minutes, mostly
# because of the cross-validation):
#   Rscript data-raw/radmix/regenerate_admixture_fixtures.R

k_values <- 2:4
seed <- 1
n_cores <- min(4L, parallel::detectCores())

bed_file <- file.path("data-raw", "radmix", "data", "EUR_subset.bed")
if (!file.exists(bed_file)) {
  stop("Run data-raw/radmix/download_eur_subset.R first")
}
eur <- tidypopgen::gen_tibble(bed_file,
  backingfile = tempfile("eur_"),
  quiet = TRUE
)

# use the ADMIXTURE binary on the PATH if there is one, otherwise the conda
# environment created by tidygenclust
conda_env <- if (nzchar(Sys.which("admixture"))) "none" else "ctidygenclust"

runs <- lapply(k_values, function(k) {
  message("Running ADMIXTURE with K = ", k)
  adm <- tidypopgen::gt_admixture(eur,
    k = k, seed = seed, n_cores = n_cores,
    # ADMIXTURE's 5-fold cross-validation; this does not change the
    # estimates (Q, P and log-likelihood of the full data), but makes each
    # run take several times longer
    crossval = TRUE,
    conda_env = conda_env
  )
  # ADMIXTURE writes its estimates with 6 decimal places, so we can store
  # them exactly (and much more compactly) as integers in millionths
  list(
    k = k,
    loglik = adm$loglik,
    cv = adm$cv,
    q_millionths = matrix(
      as.integer(round(as.matrix(adm$Q[[1]]) * 1e6)),
      ncol = k
    ),
    p_millionths = matrix(
      as.integer(round(as.matrix(adm$P[[1]]) * 1e6)),
      ncol = k
    ),
    log = adm$log[[1]]
  )
})

version_line <- grep("ADMIXTURE Version", runs[[1]]$log, value = TRUE)
fixture <- list(
  description = paste(
    "ADMIXTURE estimates for the EUR_subset data of SnpArrays.jl, used by",
    "tests/testthat/test_radmix_vs_admixture.R, including the 5-fold cross-validation",
    "error (cv) of each run. Q and P are stored in millionths;",
    "P is the frequency of the second allele of the bim file (allele_ref",
    "for a gen_tibble), with loci in the order of the bim file."
  ),
  admixture_version = trimws(gsub("\\*", "", version_line)),
  seed = seed,
  # MD5 checksums of the input files, so that the tests can check that the
  # data have not changed
  md5 = tools::md5sum(
    paste0(sub("bed$", "", bed_file), c("bed", "bim", "fam"))
  ),
  created = Sys.Date(),
  runs = runs
)
names(fixture$md5) <- c("bed", "bim", "fam")

fixture_dir <- file.path("tests", "testthat", "fixtures")
dir.create(fixture_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(fixture, file.path(fixture_dir, "admixture_eur_reference.rds"),
  compress = "xz"
)
message(
  "Saved ", file.path(fixture_dir, "admixture_eur_reference.rds"), " (",
  round(file.size(file.path(fixture_dir, "admixture_eur_reference.rds")) / 1024),
  " KB)"
)

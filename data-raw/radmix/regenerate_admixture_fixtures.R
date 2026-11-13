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
# The script makes two fixtures: the unsupervised runs (K = 2 to 4,
# admixture_eur_reference.rds) and a supervised run (K = 3,
# admixture_eur_supervised.rds, which needs the first one to choose its
# reference individuals, see the end of the script). To make only one of
# them, give its name:
#   Rscript data-raw/radmix/regenerate_admixture_fixtures.R unsupervised
#   Rscript data-raw/radmix/regenerate_admixture_fixtures.R supervised

parts <- commandArgs(trailingOnly = TRUE)
if (length(parts) == 0) parts <- c("unsupervised", "supervised")

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
fixture_dir <- file.path("tests", "testthat", "fixtures")
dir.create(fixture_dir, showWarnings = FALSE, recursive = TRUE)
# MD5 checksums of the input files, so that the tests can check that the
# data have not changed
data_md5 <- tools::md5sum(
  paste0(sub("bed$", "", bed_file), c("bed", "bim", "fam"))
)
names(data_md5) <- c("bed", "bim", "fam")

if ("unsupervised" %in% parts) {
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
  md5 = data_md5,
  created = Sys.Date(),
  runs = runs
)

saveRDS(fixture, file.path(fixture_dir, "admixture_eur_reference.rds"),
  compress = "xz"
)
message(
  "Saved ", file.path(fixture_dir, "admixture_eur_reference.rds"), " (",
  round(file.size(file.path(fixture_dir, "admixture_eur_reference.rds")) / 1024),
  " KB)"
)
}

# Supervised analysis (ADMIXTURE --supervised), K = 3. EUR_subset has no
# population labels, so the reference individuals are chosen from the
# unsupervised ADMIXTURE fit at K = 3 stored above: every individual with
# an ancestry proportion above 0.9 in a cluster is a reference of that
# cluster (15, 59 and 28 individuals); the other 277 are estimated. The
# clusters 1, 2 and 3 are labelled "mid", "zeta" and "alpha", so that the
# order of first appearance in the .pop file (zeta, alpha, mid), which
# ADMIXTURE uses for its columns, differs from the alphabetical order.
# gt_admixture() has no supervised mode, so ADMIXTURE is called directly
# (with cross-validation, which ADMIXTURE allows in supervised mode).
if ("supervised" %in% parts) {
  unsup <- readRDS(file.path(fixture_dir, "admixture_eur_reference.rds"))
  q3 <- unsup$runs[[which(k_values == 3)]]$q_millionths / 1e6
  sup_k <- 3
  sup_threshold <- 0.9
  sup_names <- c("mid", "zeta", "alpha")
  is_ref <- apply(q3, 1, max) > sup_threshold
  labels <- sup_names[apply(q3, 1, which.max)]
  labels[!is_ref] <- NA

  # (the binary must be on the PATH, e.g. with the conda environment
  # ctidygenclust activated)
  admixture <- Sys.which("admixture")
  if (!nzchar(admixture)) {
    stop("The supervised fixture needs the admixture binary on the PATH")
  }
  run_dir <- tempfile("admixture_supervised_")
  dir.create(run_dir)
  for (ext in c("bed", "bim", "fam")) {
    file.copy(sub("bed$", ext, bed_file),
      file.path(run_dir, paste0("eur.", ext)))
  }
  # the .pop file: the population of each reference, "-" for the others
  writeLines(ifelse(is.na(labels), "-", labels),
    file.path(run_dir, "eur.pop"))
  old_wd <- setwd(run_dir)
  log <- system2(admixture,
    c("--supervised", "--cv", paste0("-s", seed), "eur.bed", sup_k,
      paste0("-j", n_cores)),
    stdout = TRUE, stderr = TRUE
  )
  setwd(old_wd)
  q_sup <- as.matrix(utils::read.table(
    file.path(run_dir, paste0("eur.", sup_k, ".Q"))
  ))
  sup_fixture <- list(
    description = paste(
      "ADMIXTURE --supervised estimates (K = 3) for the EUR_subset data of",
      "SnpArrays.jl, used by tests/testthat/test_radmix_vs_admixture.R.",
      "Reference individuals: q > 0.9 in a cluster of the unsupervised",
      "ADMIXTURE fit at K = 3 (admixture_eur_reference.rds), labelled",
      "'mid', 'zeta' and 'alpha' for clusters 1 to 3 (labels, NA for the",
      "individuals to estimate). Q is stored in millionths, with its columns",
      "in ADMIXTURE's order (first appearance in the .pop file); cv is the",
      "5-fold cross-validation error."
    ),
    admixture_version = trimws(gsub("\\*", "",
      grep("ADMIXTURE Version", log, value = TRUE))),
    seed = seed,
    md5 = data_md5,
    created = Sys.Date(),
    k = sup_k,
    labels = labels,
    loglik = as.numeric(sub("Loglikelihood: ", "",
      grep("^Loglikelihood:", log, value = TRUE))),
    n_iter = as.integer(sub("Converged in ([0-9]+) iterations.*", "\\1",
      grep("^Converged in", log, value = TRUE))),
    cv = as.numeric(sub(".*: ", "", grep("^CV error", log, value = TRUE))),
    q_millionths = matrix(as.integer(round(q_sup * 1e6)), ncol = sup_k),
    log = log
  )
  saveRDS(sup_fixture, file.path(fixture_dir, "admixture_eur_supervised.rds"),
    compress = "xz"
  )
  message("Saved ", file.path(fixture_dir, "admixture_eur_supervised.rds"))
}

# Regenerate the NGSadmix results used by tests/testthat/test_radmix_vs_ngsadmix.R
#
# With genotype likelihoods (given by another package through
# gt_radmix_data()), gt_radmix() fits the model of NGSadmix (Skotte et al.
# 2013). NGSadmix is
# distributed with ANGSD and is not under active development, so rather than
# running it every time the tests run, we run it once on simulated genotype
# likelihoods and store its estimates in
# tests/testthat/fixtures/ngsadmix_gl_reference.rds, together with the
# likelihoods themselves (as 16-bit codes, see lik_codes() below), so that
# the test depends neither on NGSadmix nor on R's random number generator.
#
# Requirements: NGSadmix (from ANGSD 0.940), on the PATH; set NGSADMIX to
# the path of the executable to use another one.
#
# Run from the root of the package (this takes a few seconds):
#   Rscript data-raw/radmix/regenerate_ngsadmix_fixtures.R

# simulation: 60 individuals and 600 loci from K = 3 ancestral populations,
# with reads at a mean depth of 3 (so that about 5% of the genotypes have no
# reads, and most of the others are uncertain)
n_indiv <- 60
n_loci <- 600
k_sim <- 3
depth <- 3
sim_seed <- 11
# NGSadmix is run for each K with several seeds, and the best run (by
# log-likelihood) is kept
k_values <- 2:3
ngs_seeds <- 1:10

# the simulation helpers of the tests (radmix_simulate_admixed() and
# radmix_simulate_gl())
source(file.path("tests", "testthat", "helper-radmix-data.R"))

# The likelihoods are stored compactly, as 16-bit codes: two codes per
# individual and locus, the log-ratios of the likelihoods of 1 and 2 copies
# of allele2 to that of 0 copies, in units of 1/1024 nat, offset by 32768,
# after clamping each log-likelihood to about 32 nats below the best. This
# encodes an array of natural-log likelihoods (individuals x loci x 3, -Inf
# allowed) into a matrix with two rows per individual (u1, then u2) and one
# column per locus; both codes 32768 mean no data (three equal
# likelihoods).
lik_codes <- function(ll) {
  floor_c <- -32767 / 1024
  mx <- pmax(ll[, , 1], ll[, , 2], ll[, , 3])
  # each log-likelihood relative to the best one, clamped (which also turns
  # -Inf into the floor)
  x <- lapply(1:3, function(gg) pmax(ll[, , gg] - mx, floor_c))
  codes <- matrix(0L, 2 * nrow(mx), ncol(mx))
  # (round() rounds half to even)
  code <- function(d) 32768L + as.integer(round(d * 1024))
  codes[c(TRUE, FALSE), ] <- code(x[[2]] - x[[1]])
  codes[c(FALSE, TRUE), ] <- code(x[[3]] - x[[1]])
  codes
}

sim <- radmix_simulate_admixed(
  n = n_indiv, m = n_loci, k = k_sim, missing = 0, seed = sim_seed
)
gl_ll <- radmix_simulate_gl(sim$genotypes, depth = depth)
# the codes (two rows per individual, one column per locus); NGSadmix reads
# the likelihoods that they decode to, as does the test, so that both
# programs fit the same data
codes <- lik_codes(gl_ll)

# the beagle file read by NGSadmix: for each individual, the likelihoods of
# the three genotypes (0, 1 and 2 copies of allele2) normalised to
# sum to 1, decoded from the codes, with 10 significant digits (far below
# the precision of the encoding, 5e-4); an individual without data has
# 1/3 for each
e1 <- exp((codes[c(TRUE, FALSE), ] - 32768) / 1024)
e2 <- exp((codes[c(FALSE, TRUE), ] - 32768) / 1024)
total <- 1 + e1 + e2
lik <- array(c(1 / total, e1 / total, e2 / total), c(dim(e1), 3))
values <- matrix(
  sprintf("%.10g", aperm(lik, c(3, 1, 2))), n_loci,
  byrow = TRUE
)
work_dir <- tempfile("ngsadmix_")
dir.create(work_dir)
beagle_file <- file.path(work_dir, "sim.beagle.gz")
con <- gzfile(beagle_file, "w")
writeLines(paste(
  c("marker", "allele1", "allele2", rep(paste0("Ind", seq_len(n_indiv) - 1),
    each = 3
  )),
  collapse = "\t"
), con)
writeLines(paste(paste0("1_", seq_len(n_loci)), 0, 3,
  apply(values, 1, paste, collapse = "\t"),
  sep = "\t"
), con)
close(con)

# NGSadmix: the one in NGSADMIX, or on the PATH
ngsadmix <- Sys.getenv("NGSADMIX", Sys.which("NGSadmix"))
if (!nzchar(ngsadmix) || !file.exists(ngsadmix)) {
  stop("NGSadmix not found: put it on the PATH, or set NGSADMIX")
}

# NGSadmix's options: no filter on the minor allele frequency (by default it
# drops loci below 0.05), stringent convergence criteria (by default it
# stops when the log-likelihood changes by less than 0.1 in 50 iterations),
# and a single thread
run_ngsadmix <- function(k, seed) {
  prefix <- file.path(work_dir, paste0("k", k, "_s", seed))
  out <- system2(ngsadmix,
    c(
      "-likes", beagle_file, "-K", k, "-minMaf", 0, "-seed", seed,
      "-tol", 1e-10, "-tolLike50", 1e-8, "-maxiter", 10000, "-P", 1,
      "-o", prefix
    ),
    stdout = TRUE, stderr = TRUE
  )
  # (without the temporary directory, in the paths of the files that the
  # log names)
  log <- gsub(paste0(work_dir, "/"), "", readLines(paste0(prefix, ".log")),
    fixed = TRUE
  )
  like_line <- grep("best like=", log, value = TRUE)
  list(
    seed = seed,
    loglik = as.numeric(sub(".*best like=([-0-9.e]+).*", "\\1", like_line)),
    q = unname(as.matrix(utils::read.table(paste0(prefix, ".qopt")))),
    f = unname(as.matrix(utils::read.table(
      gzfile(paste0(prefix, ".fopt.gz"))
    ))),
    log = log
  )
}

runs <- lapply(k_values, function(k) {
  message("Running NGSadmix with K = ", k)
  all_runs <- lapply(ngs_seeds, run_ngsadmix, k = k)
  logliks <- vapply(all_runs, function(r) r$loglik, numeric(1))
  best <- all_runs[[which.max(logliks)]]
  message(
    "  log-likelihoods: ",
    paste(format(sort(logliks, decreasing = TRUE), nsmall = 3), collapse = ", ")
  )
  c(list(k = k, all_logliks = logliks), best)
})

fixture <- list(
  description = paste(
    "NGSadmix estimates for simulated genotype likelihoods, used by",
    "tests/testthat/test_radmix_vs_ngsadmix.R. `codes` holds 16-bit codes",
    "of the likelihoods (two rows per individual, one column per locus; see",
    "lik_codes() in data-raw/radmix/regenerate_ngsadmix_fixtures.R);
    "NGSadmix read the normalised likelihoods that they decode to. For each",
    "K, the run with the best log-likelihood among the seeds: q",
    "(individuals x K), f (loci x K, the frequency of allele2, the minor",
    "allele) and NGSadmix's log-likelihood, which counts log(1/3) for each",
    "individual without data at a locus."
  ),
  ngsadmix = basename(ngsadmix),
  options = "-minMaf 0 -tol 1e-10 -tolLike50 1e-8 -maxiter 10000 -P 1",
  simulation = list(
    n_indiv = n_indiv, n_loci = n_loci, k = k_sim, depth = depth,
    seed = sim_seed
  ),
  created = Sys.Date(),
  codes = codes,
  runs = runs
)

fixture_dir <- file.path("tests", "testthat", "fixtures")
fixture_file <- file.path(fixture_dir, "ngsadmix_gl_reference.rds")
saveRDS(fixture, fixture_file, compress = "xz")
message(
  "Saved ", fixture_file, " (", round(file.size(fixture_file) / 1024), " KB)"
)
unlink(work_dir, recursive = TRUE)

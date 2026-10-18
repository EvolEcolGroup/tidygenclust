# Regenerate the OpenADMIXTURE.jl results used by tests/testthat/test_openadmixture_vs_julia.R
#
# OpenADMIXTURE.jl is no longer actively developed, so rather than running
# Julia every time the tests run, we run it once on a small simulated dataset
# and store its results, together with the data and the starting values, in
# tests/testthat/fixtures/openadmixture_julia_reference.rds. The tests then run the
# corresponding functions of gt_openadmixture() and compare. Rerun this script
# only if the comparison needs to change.
#
# Starting from the same values of Q and P, data-raw/openadmixture/julia/openadmixture_reference.jl
# computes with OpenADMIXTURE.jl: the log-likelihood, one EM step, one round
# of Newton updates of Q and P, and the full quasi-Newton algorithm (with the
# log-likelihood at each iteration).
#
# Requirements: Julia with OpenADMIXTURE.jl. The easiest way to get it is the
# conda environment created by data-raw/openadmixture/setup_julia_env.sh (the name of the
# environment can be changed with the environment variable OPENADMIXTURE_JULIA_ENV);
# alternatively, set OPENADMIXTURE_JULIA to a julia executable whose default
# environment includes OpenADMIXTURE.jl.
#
# Run from the root of the package:
#   Rscript data-raw/openadmixture/regenerate_julia_fixtures.R

# Julia command ---------------------------------------------------------------

# Julia is started with --threads=1,0 (see data-raw/openadmixture/julia/openadmixture_reference.jl);
# a single thread also makes the results of OpenADMIXTURE.jl deterministic
julia_args <- "--threads=1,0"
if (nzchar(Sys.getenv("OPENADMIXTURE_JULIA"))) {
  julia_cmd <- Sys.getenv("OPENADMIXTURE_JULIA")
} else {
  # run julia through `conda run`, so that the Julia depot of the conda
  # environment is used
  julia_cmd <- Sys.getenv("CONDA_EXE", Sys.which("conda"))
  if (!nzchar(julia_cmd)) {
    stop("conda not found; set OPENADMIXTURE_JULIA or CONDA_EXE")
  }
  julia_args <- c(
    "run", "--no-capture-output",
    "-n", Sys.getenv("OPENADMIXTURE_JULIA_ENV", "openadmixturer_julia"),
    "julia", julia_args
  )
}

# Data and starting values ----------------------------------------------------

# Simulated admixed genotypes: 150 individuals and 300 loci (more than the
# 64 x 64 tiles used by both implementations), K = 4, with 3% missing
# genotypes. The genotypes themselves are stored in the fixture, so the tests
# do not depend on the random number generator.
source(file.path("tests", "testthat", "helper-openadmixture-data.R"))
sim <- openadmixture_simulate_admixed(n = 150, m = 300, k = 4, missing = 0.03, seed = 21)
genotypes <- sim$genotypes
k <- 4

# starting values in the internal layout of both packages (K x I and K x J)
set.seed(5)
q0 <- matrix(stats::runif(k * nrow(genotypes)), k)
q0 <- sweep(q0, 2, colSums(q0), "/")
p0 <- matrix(stats::runif(k * ncol(genotypes), 0.05, 0.95), k)

# Run OpenADMIXTURE.jl --------------------------------------------------------

work_dir <- tempfile("openadmixture_julia_")
dir.create(work_dir)
gt <- openadmixture_test_gt(genotypes, backingfile = file.path(work_dir, "sim_gt"))
bed_file <- tidypopgen::gt_as_plink(gt,
  file = file.path(work_dir, "sim"),
  chromosomes_as_int = TRUE
)

# SnpArrays (and hence OpenADMIXTURE.jl) counts copies of the second allele
# of the bim file, which gt_as_plink() sets to allele_ref, whereas
# gt_openadmixture() counts allele_alt: Julia's genotypes are 2 - g, and its
# allele frequencies 1 - p. The likelihood is symmetric under this change, so
# we pass 1 - p0 to Julia and flip its P matrices back below; Q is unaffected.
write_csv_matrix <- function(x, name) {
  utils::write.table(x, file.path(work_dir, name),
    sep = ",",
    row.names = FALSE, col.names = FALSE
  )
}
read_csv_matrix <- function(name) {
  x <- utils::read.csv(file.path(work_dir, name), header = FALSE)
  unname(as.matrix(x))
}
write_csv_matrix(q0, "q0.csv")
write_csv_matrix(1 - p0, "p0.csv")

status <- system2(
  julia_cmd,
  c(
    julia_args,
    shQuote(file.path("data-raw", "openadmixture", "julia", "openadmixture_reference.jl")),
    shQuote(bed_file), shQuote(file.path(work_dir, "q0.csv")),
    shQuote(file.path(work_dir, "p0.csv")), shQuote(work_dir), "1000"
  )
)
if (status != 0) stop("OpenADMIXTURE.jl failed")

# versions, for the record
julia_version <- system2(julia_cmd, c(julia_args, "--version"),
  stdout = TRUE
)
openadmixture_jl_version <- system2(
  julia_cmd,
  c(julia_args, "-e", shQuote(paste(
    "using Pkg;",
    "info = Pkg.dependencies()[Base.UUID(",
    "\"425c031d-9e5b-441b-a07e-4acde9a966a3\")];",
    "print(info.version, \" \", info.git_revision)"
  ))),
  stdout = TRUE
)

# Save the fixture ------------------------------------------------------------

fixture <- list(
  description = paste(
    "Results of OpenADMIXTURE.jl on simulated data, used by",
    "tests/testthat/test_openadmixture_vs_julia.R. Matrices are in the internal layout",
    "(K x individuals for Q, K x loci for P); P is the frequency of",
    "allele_alt (Julia's values have been flipped)."
  ),
  julia_version = julia_version,
  openadmixture_jl_version = openadmixture_jl_version,
  created = Sys.Date(),
  genotypes = genotypes,
  q0 = q0,
  p0 = p0,
  loglik = drop(read_csv_matrix("loglik.csv")),
  em_q = read_csv_matrix("em_q.csv"),
  em_p = 1 - read_csv_matrix("em_p.csv"),
  update_q = read_csv_matrix("upd_q.csv"),
  update_p = 1 - read_csv_matrix("upd_p.csv"),
  run_q = read_csv_matrix("run_q.csv"),
  run_p = 1 - read_csv_matrix("run_p.csv"),
  run_loglik = drop(read_csv_matrix("run_ll.csv"))
)
unlink(work_dir, recursive = TRUE)

fixture_file <- file.path(
  "tests", "testthat", "fixtures", "openadmixture_julia_reference.rds"
)
dir.create(dirname(fixture_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(fixture, fixture_file, compress = "xz")
message(
  "Saved ", fixture_file, " (", round(file.size(fixture_file) / 1024),
  " KB)"
)

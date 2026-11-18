# Simulated genotype likelihoods for timing gt_radmix() on a gl_tibble:
# N individuals x M loci from K = 3 ancestral populations (1% missing
# genotypes), reads at a mean depth of 4 with an error rate of 0.01, as in
# tidygenclust's data-raw/radmix/gl/benchmark_gl.R. Writes a beagle file,
# imports it with tidypopgl::gl_tibble() (gl16 backing file) and saves it
# (gl_save()).
#
# Usage: Rscript make_gl_data.R <tidygenclust lib> <tidypopgl lib> <out dir>
# Environment: GL_N (400), GL_M (100000), TGC_REPO (tidygenclust checkout,
# for the simulation helpers of its tests).
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[2], args[1], .libPaths()))
out_dir <- args[3]
n_indiv <- as.integer(Sys.getenv("GL_N", "400"))
n_loci <- as.integer(Sys.getenv("GL_M", "100000"))
tgc_repo <- Sys.getenv("TGC_REPO", "/home/andrea/git/tidygenclust")
library(tidygenclust)
ns <- asNamespace("tidygenclust")
radmix_gl_dense_source_cpp <- ns$radmix_gl_dense_source_cpp
source(file.path(tgc_repo, "tests", "testthat", "helper-radmix-data.R"))

sim <- radmix_simulate_admixed(
  n = n_indiv, m = n_loci, k = 3, missing = 0.01, seed = 1
)
set.seed(2)
ll <- radmix_simulate_gl(sim$genotypes, depth = 4)
rm(sim)
invisible(gc())

# the beagle file: normalised likelihoods (8 significant digits), one row
# per locus, alleles A (0) and T (3), markers 1_<position>
beagle_file <- file.path(out_dir, sprintf("sim_%d_%d.beagle.gz", n_indiv,
                                          n_loci))
# the .fai file that gl_tibble() reads the lengths of the chromosomes from
writeLines(paste("1", n_loci * 100 + 1000, 0, 60, 61, sep = "\t"),
           sub("\\.beagle\\.gz$", ".fai", beagle_file))
con <- gzfile(beagle_file, "w")
writeLines(paste(
  c("marker", "allele1", "allele2",
    rep(paste0("ind", seq_len(n_indiv)), each = 3)),
  collapse = "\t"
), con)
for (j0 in seq(1, n_loci, by = 5000)) {
  js <- j0:min(n_loci, j0 + 4999)
  lk <- exp(ll[, js, , drop = FALSE])
  total <- lk[, , 1] + lk[, , 2] + lk[, , 3]
  lik <- lk / as.vector(total)
  values <- matrix(
    sprintf("%.8g", aperm(lik, c(3, 1, 2))), length(js),
    byrow = TRUE
  )
  writeLines(paste(sprintf("1_%d", as.integer(js * 100)), 0, 3,
    apply(values, 1, paste, collapse = "\t"),
    sep = "\t"
  ), con)
}
close(con)
message("beagle file written")

# the gl_tibble (its backing file next to the beagle file)
gl <- tidypopgl::gl_tibble(beagle_file,
  backingfile = file.path(out_dir, sprintf("sim_%d_%d", n_indiv, n_loci)),
  quiet = TRUE
)
tidypopgl::gl_save(gl, quiet = TRUE)
message("gl_tibble saved")

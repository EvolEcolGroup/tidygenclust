# Simulated genotypes for the part of the benchmark of gt_radmix() on data
# larger than the memory, shown in the article admixture_benchmark.Rmd (in
# the parent directory of this script), which also shows the code of the
# chunks of this script (the sections that start with a "## ---- label"
# line, read with knitr::read_chunk()).
#
# 220,000 individuals typed at 200,000 loci, admixed from 5 ancestral
# populations, written as a PLINK bed file and read into a gen_tibble: its
# backing file (one byte per genotype) takes 44 GB, and the bed file (two
# bits per genotype) 11 GB. Run it once, from the directory of the article
# (vignettes/articles), on a machine with enough disk space (about 100 GB
# while it runs, 55 GB at the end) and memory (about 10 GB per core):
#
#   cd vignettes/articles
#   Rscript code/admixture_benchmark_large_data.R
#
# The data are written to data/large (or to the directory given by the
# environment variable LARGE_DIR), which is not part of the repository. The
# number of cores used to simulate them (LARGE_CORES, default 8) does not
# change the data.

## ---- large-setup ----
large_dir <- Sys.getenv("LARGE_DIR", "data/large")
dir.create(large_dir, recursive = TRUE, showWarnings = FALSE)
n_cores <- as.integer(Sys.getenv("LARGE_CORES", "8"))
n_indiv <- 220000
n_loci <- 200000
k <- 5

## ---- large-simulate ----
# Five ancestral populations, with allele frequencies drawn under the
# Balding-Nichols model (Fst = 0.1 from frequencies uniform between 0.05 and
# 0.95 in their common ancestor), and five sampled populations of 44,000
# individuals each, whose ancestry proportions are drawn from a Dirichlet
# distribution centred on one ancestral population (80% on average, the
# rest shared among the other four). The genotypes are binomial draws from
# the individual allele frequencies, with 1% of them missing.
set.seed(1)
pop <- rep(seq_len(k), each = n_indiv / k)
alpha <- matrix(0.25, n_indiv, k)
alpha[cbind(seq_len(n_indiv), pop)] <- 4
q_true <- matrix(stats::rgamma(n_indiv * k, alpha), n_indiv, k)
q_true <- q_true / rowSums(q_true)
saveRDS(q_true, file.path(large_dir, "q_true.rds"))
# the genotypes, written straight into the backing file of an FBM in blocks
# of loci (one byte per genotype, 3 for a missing genotype), on several
# cores; each block of loci has its own seed, so the data do not depend on
# the number of cores. Each core runs BLAS (for its matrix products) on a
# single thread, set in a.FUN, so bigstatsr's check against two levels of
# parallelism is not needed
options(bigstatsr.check.parallel.blas = FALSE)
genotypes <- bigstatsr::FBM.code256(
  n_indiv, n_loci,
  code = bigsnpr::CODE_012,
  backingfile = file.path(large_dir, "sim_fbm"),
  init = NULL
)
block_size <- 250
bigstatsr::big_apply(
  genotypes,
  a.FUN = function(x, ind, q_true, k, block_size) {
    bigparallelr::set_blas_ncores(1)
    for (j0 in seq(1, length(ind), by = block_size)) {
      cols <- ind[j0:min(length(ind), j0 + block_size - 1)]
      set.seed(1e6 + cols[1])
      f <- 0.1
      p_anc <- stats::runif(length(cols), 0.05, 0.95)
      p <- matrix(
        stats::rbeta(
          k * length(cols),
          rep(p_anc * (1 - f) / f, each = k),
          rep((1 - p_anc) * (1 - f) / f, each = k)
        ),
        k, length(cols)
      )
      g <- stats::rbinom(nrow(x) * length(cols), 2, q_true %*% p)
      g[stats::runif(length(g)) < 0.01] <- 3L
      x[, cols] <- matrix(g, nrow(x))
    }
    NULL
  },
  a.combine = "c",
  block.size = block_size * 8,
  ncores = n_cores,
  q_true = q_true,
  k = k,
  block_size = block_size
)
genotypes$save()

## ---- large-bed ----
# a PLINK bed file of the genotypes, with the sampled population of each
# individual as its family ID, and the loci on 20 chromosomes
sim <- structure(
  list(
    genotypes = genotypes,
    fam = data.frame(
      family.ID = paste0("pop", pop),
      sample.ID = sprintf("ind%06d", seq_len(n_indiv)),
      paternal.ID = 0L,
      maternal.ID = 0L,
      sex = 0L,
      affection = -9L
    ),
    map = data.frame(
      chromosome = rep(1:20, each = n_loci / 20),
      marker.ID = sprintf("snp%06d", seq_len(n_loci)),
      genetic.dist = 0,
      physical.pos = rep(seq_len(n_loci / 20) * 1000L, 20),
      allele1 = "A",
      allele2 = "G"
    )
  ),
  class = "bigSNP"
)
bed_file <- bigsnpr::snp_writeBed(sim, file.path(large_dir, "sim_large.bed"))
# the FBM is no longer needed: the gen_tibble gets its own backing file
rm(sim, genotypes)
invisible(gc())
unlink(file.path(large_dir, c("sim_fbm.bk", "sim_fbm.rds")))

## ---- large-gen-tibble ----
# the gen_tibble that gt_radmix() reads, saved so that each run can load it
# without reading the bed file again
sim_gt <- tidypopgen::gen_tibble(
  bed_file,
  backingfile = file.path(large_dir, "sim_large"),
  quiet = TRUE
)
invisible(tidypopgen::gt_save(sim_gt, quiet = TRUE))

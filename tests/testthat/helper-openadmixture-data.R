# Small gen_tibbles used across tests

# limit number of threads for tests
data.table::setDTthreads(2)

openadmixture_test_gt <- function(genotypes,
                                  backingfile = tempfile("openadmixture_")) {
  n_loci <- ncol(genotypes)
  loci <- data.frame(
    name = paste0("rs", seq_len(n_loci)),
    chromosome = rep(c("chr1", "chr2"), length.out = n_loci),
    position = as.integer(seq_len(n_loci) * 10),
    genetic_dist = as.double(rep(0, n_loci)),
    allele_ref = rep("A", n_loci),
    allele_alt = rep("T", n_loci)
  )
  indiv_meta <- data.frame(
    id = paste0("ind", seq_len(nrow(genotypes))),
    population = rep(c("pop1", "pop2"), length.out = nrow(genotypes))
  )
  tidypopgen::gen_tibble(
    x = genotypes,
    loci = loci,
    indiv_meta = indiv_meta,
    backingfile = backingfile,
    quiet = TRUE
  )
}

# simulate admixed genotypes with known Q (n x K) and P (K x m)
openadmixture_simulate_admixed <- function(n = 60, m = 200, k = 3,
                                           missing = 0.02, seed = 1) {
  set.seed(seed)
  p <- matrix(stats::rbeta(k * m, 0.5, 0.5), k, m)
  q <- matrix(stats::rgamma(n * k, 0.3), n, k)
  q <- q / rowSums(q)
  prob <- q %*% p
  g <- matrix(stats::rbinom(n * m, 2, prob), n, m)
  g[sample(length(g), round(missing * length(g)))] <- NA
  list(genotypes = g, q = q, p = p)
}

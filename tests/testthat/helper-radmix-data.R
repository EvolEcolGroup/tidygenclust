# Small gen_tibbles used across tests

# limit number of threads for tests
data.table::setDTthreads(2)

# a gen_tibble (with its own backing file) from a matrix of alternate allele
# counts (individuals x loci, NA for missing), with loci rs1, rs2, ... on two
# chromosomes and individuals ind1, ind2, ... in two populations
radmix_test_gt <- function(genotypes,
                           backingfile = tempfile("radmix_")) {
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
# (this sets the seed of the global random number generator)
radmix_simulate_admixed <- function(n = 60, m = 200, k = 3,
                                    missing = 0.02, seed = 1) {
  set.seed(seed)
  # U-shaped allele frequencies, many close to 0 or 1, as in real data
  p <- matrix(stats::rbeta(k * m, 0.5, 0.5), k, m)
  # normalised Gamma(0.3) draws, i.e. Dirichlet(0.3): most individuals have
  # most of their ancestry from one population
  q <- matrix(stats::rgamma(n * k, 0.3), n, k)
  q <- q / rowSums(q)
  # alternate allele counts of diploids, Binomial(2, (QP)_ij); a random
  # fraction `missing` of them is then set to NA
  prob <- q %*% p
  g <- matrix(stats::rbinom(n * m, 2, prob), n, m)
  g[sample(length(g), round(missing * length(g)))] <- NA
  list(genotypes = g, q = q, p = p)
}

# a gen_tibble in which the individuals `hap` (indices of rows of
# `genotypes`) are pseudohaploid: their heterozygous genotypes are made
# homozygous for the alternate allele, so that tidypopgen::gt_pseudohaploid()
# finds them pseudohaploid
radmix_pseudohaploid_gt <- function(genotypes, hap) {
  genotypes[hap, ] <- 2 * (genotypes[hap, ] >= 1)
  tidypopgen::gt_pseudohaploid(radmix_test_gt(genotypes))
}

# Reorder the columns of q to best match q_ref (the labels of the clusters
# are arbitrary); returns the permutation
radmix_match_clusters <- function(q, q_ref) {
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

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

# Genotype likelihoods of the genotypes g (copies of the alternate allele,
# individuals x loci, NA = no data): reads at each genotype, Poisson(depth)
# of them (none where g is NA), each showing the alternate allele with
# probability g / 2 (1 - err) + (1 - g / 2) err. Returns the natural-log
# likelihoods of 0, 1 and 2 copies of the alternate allele, an array
# individuals x loci x 3; a genotype without reads has three equal
# likelihoods (no data).
radmix_simulate_gl <- function(g, depth = 3, err = 0.01) {
  n <- nrow(g)
  m <- ncol(g)
  dp <- matrix(stats::rpois(n * m, depth), n, m)
  dp[is.na(g)] <- 0
  g0 <- ifelse(is.na(g), 0, g)
  alt <- matrix(
    stats::rbinom(n * m, dp, g0 / 2 * (1 - err) + (1 - g0 / 2) * err),
    n, m
  )
  ll <- array(0, c(n, m, 3))
  for (gg in 0:2) {
    pa <- gg / 2 * (1 - err) + (1 - gg / 2) * err
    ll[, , gg + 1] <- alt * log(pa) + (dp - alt) * log(1 - pa)
  }
  ll
}

# A dataset of genotype likelihoods, as the method of gt_radmix_data() of a
# package that stores them returns, from
# an array of natural-log likelihoods (individuals x loci x 3, -Inf
# allowed; see radmix_simulate_gl()). The likelihoods reach the kernels
# through the example source of tidygenclust (src/radmix_gl_source.cpp),
# which holds them in memory: those of 1 and 2 copies of the alternate
# allele relative to 0 copies (e1, e2), each first clamped to about 32 nats
# below the best (so that a likelihood of 0 becomes a small, positive one),
# and whether there are data (three equal likelihoods mean none). The list
# also keeps them as matrices, for the R references.
radmix_test_gl <- function(ll) {
  mx <- pmax(ll[, , 1], ll[, , 2], ll[, , 3])
  x <- lapply(1:3, function(gg) pmax(ll[, , gg] - mx, -32))
  e1 <- exp(x[[2]] - x[[1]])
  e2 <- exp(x[[3]] - x[[1]])
  obs <- !(x[[1]] == x[[2]] & x[[1]] == x[[3]])
  c(
    list(kind = "gl_source", version = 1L),
    radmix_gl_dense_source_cpp(e1, e2, obs),
    list(
      n_indiv = nrow(e1), n_loci = ncol(e1),
      e1 = e1, e2 = e2, obs = obs
    )
  )
}

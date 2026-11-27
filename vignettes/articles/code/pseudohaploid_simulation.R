# Simulation of diploid and pseudohaploid individuals from a few populations,
# used in the article pseudohaploids.Rmd (in the parent directory of this
# script), which shows the code of the chunks of this script (the sections
# that start with a "## ---- label" line, read with knitr::read_chunk()).
# The function can also be used on its own, after
# source("code/pseudohaploid_simulation.R").

## ---- simulate-function ----
# Simulate diploid and pseudohaploid individuals from a few populations.
# n_diploid and n_pseudohaploid are named vectors with the number of
# individuals of each population (a population missing from one of them has
# none of that type), e.g. n_diploid = c(A = 20, B = 20) and
# n_pseudohaploid = c(A = 40). The allele frequencies of each population
# drift from those of a common ancestral population (with the given F_ST).
# Returns a gen_tibble with the columns id, population and type ("diploid" or
# "pseudohaploid"), in which the pseudohaploids are coded 0/2 but not
# declared as pseudohaploids, as when they are read from a file. The true
# allele frequencies of the populations (a matrix with one row per locus and
# one column per population) are in its attribute "frequencies"
simulate_pseudohaploids <- function(n_diploid, n_pseudohaploid = c(),
                                    n_loci = 5000, fst = 0.1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  populations <- union(names(n_diploid), names(n_pseudohaploid))
  # the number of individuals of a population in n (0 if it is not in n)
  count <- function(n, pop) if (pop %in% names(n)) n[[pop]] else 0
  # the frequency of the alternate allele of each locus in the ancestral
  # population
  p_ancestral <- stats::runif(n_loci, 0.05, 0.95)
  frequencies <- list()
  genotypes <- list()
  indiv_meta <- list()
  for (pop in populations) {
    # the frequencies in this population, drifted from the ancestral ones
    # (Balding-Nichols model: a beta distribution with mean p_ancestral and
    # variance fst * p_ancestral * (1 - p_ancestral))
    p_pop <- frequencies[[pop]] <- stats::rbeta(
      n_loci,
      p_ancestral * (1 - fst) / fst,
      (1 - p_ancestral) * (1 - fst) / fst
    )
    # n individuals with `copies` copies of each locus, each copy carrying
    # the alternate allele with probability p_pop (one row per individual)
    draw <- function(n, copies) {
      matrix(
        stats::rbinom(n * n_loci, copies, p_pop),
        ncol = n_loci, byrow = TRUE
      )
    }
    # a diploid carries two copies; a pseudohaploid carries one, stored as a
    # homozygous genotype (0 or 2), as in tidypopgen
    n_dip <- count(n_diploid, pop)
    n_pseudo <- count(n_pseudohaploid, pop)
    genotypes[[pop]] <- rbind(draw(n_dip, 2), 2 * draw(n_pseudo, 1))
    indiv_meta[[pop]] <- data.frame(
      population = rep(pop, n_dip + n_pseudo),
      type = rep(c("diploid", "pseudohaploid"), c(n_dip, n_pseudo))
    )
  }
  indiv_meta <- do.call(rbind, indiv_meta)
  indiv_meta <- data.frame(
    id = paste0("ind", seq_len(nrow(indiv_meta))),
    indiv_meta
  )
  loci <- data.frame(
    name = paste0("rs", seq_len(n_loci)),
    chromosome = 1L,
    position = seq_len(n_loci),
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
  gt <- tidypopgen::gen_tibble(
    do.call(rbind, genotypes),
    indiv_meta = indiv_meta,
    loci = loci,
    backingfile = tempfile("pseudohaploids_"),
    quiet = TRUE
  )
  attr(gt, "frequencies") <- do.call(cbind, frequencies)
  gt
}

# Why pseudohaploid samples must be declared as such
#
# Low-coverage (e.g. ancient DNA) samples are often "pseudohaploid": a single
# read is drawn at each site, so only one allele is observed. tidypopgen
# stores such genotypes as 0/2 (homozygous), and gt_pseudohaploid() marks the
# individuals concerned (ploidy 1), so that gt_radmix() treats each of
# their genotypes as a single allele (a binomial with n = 1). This script
# simulates data where forgetting to do so, i.e. treating pseudohaploids as
# diploids, gives the wrong answer.
#
# What goes wrong, and what does not
# ----------------------------------
# Treated as diploid, a pseudohaploid genotype of 0 or 2 counts as two
# identical allele copies. For a single individual with given allele
# frequencies this only squares its likelihood, which doubles the
# log-likelihood without moving its maximum: the ancestry proportions
# estimated at a given K are essentially unaffected. But the size of the
# likelihood matters for model selection: cross-validation scores held-out
# 0/2 genotypes as two allele copies that are never heterozygous, and with
# many undeclared pseudohaploids it fails to detect real population
# structure, preferring K = 1 over the true K = 2.
#
# Simulation
# ----------
# Two diploid populations, A and B (20 individuals each), with allele
# frequencies drawn from the Balding-Nichols model (F_ST = 0.05 from a
# common ancestral population) at 5,000 SNPs, plus 40 pseudohaploid
# individuals from population A (one allele drawn at each SNP, stored as
# 0/2). The true number of populations is 2. Each replicate is analysed with
# gt_radmix() for K = 1 to 4 with cross-validation, once treating
# every individual as diploid and once after gt_pseudohaploid(). Results
# from the development laptop (5 replicates):
# - treated as diploid: cross-validation picks K = 1 in every replicate;
# - pseudohaploid-aware: cross-validation picks K = 2 in every replicate;
# - at K = 2 both treatments estimate the same ancestry for the
#   pseudohaploids (about 96% from population A).
#
# The results are written to data-raw/radmix/pseudohaploid_vs_diploid.csv.
#
# Run from the root of the package, with tidygenclust installed (this takes
# a few minutes):
#   Rscript data-raw/radmix/pseudohaploid_vs_diploid.R

library(tidygenclust)
library(tidypopgen)

n_loci <- 5000
fst <- 0.05
n_a <- 20
n_b <- 20
n_pseudo <- 40
n_replicates <- 5
k_values <- 1:4
n_cores <- min(4L, parallel::detectCores())

# Simulate one dataset, returning a gen_tibble with a population column
simulate_dataset <- function(seed) {
  set.seed(seed)
  # allele frequencies: a common ancestral population, and two populations
  # drifted from it (Balding-Nichols model)
  p_anc <- stats::runif(n_loci, 0.05, 0.95)
  balding_nichols <- function(p) {
    stats::rbeta(length(p), p * (1 - fst) / fst, (1 - p) * (1 - fst) / fst)
  }
  p_a <- balding_nichols(p_anc)
  p_b <- balding_nichols(p_anc)
  # diploid genotypes: two allele copies per locus
  diploid <- function(n, p) t(replicate(n, stats::rbinom(length(p), 2, p)))
  # pseudohaploid genotypes from population A: a single allele per locus,
  # stored as 0/2 as tidypopgen does
  pseudohaploid <- t(replicate(n_pseudo, 2 * stats::rbinom(n_loci, 1, p_a)))
  genotypes <- rbind(diploid(n_a, p_a), diploid(n_b, p_b), pseudohaploid)
  population <- c(
    rep("A", n_a), rep("B", n_b), rep("A (pseudohaploid)", n_pseudo)
  )
  loci <- data.frame(
    name = paste0("rs", seq_len(n_loci)), chromosome = 1L,
    position = seq_len(n_loci), genetic_dist = 0,
    allele_ref = "A", allele_alt = "T"
  )
  gen_tibble(genotypes,
    loci = loci,
    indiv_meta = data.frame(
      id = paste0("ind", seq_along(population)), population = population
    ),
    backingfile = tempfile("pseudo_"), quiet = TRUE
  )
}

results <- list()
for (replicate in seq_len(n_replicates)) {
  gt <- simulate_dataset(seed = replicate)
  # gt_pseudohaploid() finds individuals without heterozygous genotypes and
  # sets their ploidy to 1; here it should find exactly the 40 simulated
  # pseudohaploids
  gt_pseudo <- gt_pseudohaploid(gt)
  stopifnot(sum(indiv_ploidy(gt_pseudo) == 1) == n_pseudo)
  is_pseudo <- gt$population == "A (pseudohaploid)"

  for (treatment in c("treated as diploid", "pseudohaploid-aware")) {
    x <- if (treatment == "treated as diploid") gt else gt_pseudo
    adm <- gt_radmix(x,
      k = k_values, seed = replicate, crossval = TRUE,
      n_cores = n_cores
    )
    # ancestry of the pseudohaploids at K = 2, from the component that
    # dominates in the diploid individuals of population A
    q2 <- unclass(get_q_matrix(adm, k = 2, run = 1))
    comp_a <- which.max(colMeans(q2[gt$population == "A", , drop = FALSE]))
    row <- data.frame(
      replicate = replicate,
      treatment = treatment,
      best_k = adm$k[which.min(adm$cv)],
      pseudo_ancestry_a_k2 = round(mean(q2[is_pseudo, comp_a]), 3)
    )
    # cross-validation error of each K
    for (i in seq_along(adm$k)) {
      row[[paste0("cv_k", adm$k[i])]] <- round(adm$cv[i], 4)
    }
    print(row, row.names = FALSE)
    results[[length(results) + 1]] <- row
  }
}

results <- do.call(rbind, results)
utils::write.csv(results,
  file.path("data-raw", "radmix", "pseudohaploid_vs_diploid.csv"),
  row.names = FALSE
)

# summary: how often each treatment picks each K
print(table(treatment = results$treatment, best_k = results$best_k))

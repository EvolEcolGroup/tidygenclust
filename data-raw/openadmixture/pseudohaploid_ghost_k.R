# Undeclared pseudohaploids create "ghost" clusters when K is chosen by AIC
#
# A companion to pseudohaploid_vs_diploid.R, which shows that treating
# pseudohaploid samples as diploids can make cross-validation miss real
# population structure (K = 1 instead of K = 2). This script shows the
# opposite error: with a criterion based on the log-likelihood, such as the
# Akaike information criterion (AIC), undeclared pseudohaploids make extra,
# spurious ("ghost") clusters look supported.
#
# Why
# ---
# Treated as diploid, a pseudohaploid genotype of 0 or 2 counts as two
# identical allele copies, so the log-likelihood of a pseudohaploid
# individual is exactly twice the one it should have: the model believes it
# has twice as much information about these individuals as it really has.
# Every additional cluster improves the fit to the noise of the sample, and
# for undeclared pseudohaploids this improvement is counted twice, which is
# enough to outweigh the penalty of AIC for the extra parameters
# (K x (loci + individuals)).
#
# Cross-validation, the criterion used by ADMIXTURE and by
# gt_openadmixture(crossval = TRUE), is not fooled in this way: it scores
# held-out genotypes, so an improvement that comes from fitting noise does
# not carry over to them. Doubling the weight of the pseudohaploids only
# changes how much each held-out genotype counts, and in our simulations it
# made cross-validation more reluctant, not less, to add clusters (see
# pseudohaploid_vs_diploid.R). We could not find a design in which
# cross-validation picks a ghost cluster because of undeclared
# pseudohaploids alone; the risk is for users who choose K from the
# log-likelihood, e.g. with AIC.
#
# Simulation
# ----------
# Two populations, A and B, with allele frequencies drawn from the
# Balding-Nichols model (F_ST = 0.1 from a common ancestral population) at
# 5,000 SNPs; pseudohaploid individuals have one allele drawn at each SNP,
# stored as 0/2. The true number of populations is 2. Two scenarios:
# 1. "pseudohaploids within A": as in pseudohaploid_vs_diploid.R, but with
#    stronger structure so that cross-validation finds K = 2 with both
#    treatments: 20 diploid individuals from each population, plus 40
#    pseudohaploid individuals from population A;
# 2. "B all pseudohaploid": 20 diploid individuals from A and 40
#    pseudohaploid individuals from B, e.g. a population known only from
#    ancient samples.
# Each replicate is analysed with gt_openadmixture() for K = 1 to 6 with
# cross-validation, once treating every individual as diploid and once
# after gt_pseudohaploid(). AIC is computed as -2 log-likelihood + 2 x the
# number of free parameters, K x loci (P) + (K - 1) x individuals (Q). At the
# K picked by AIC, we also count the clusters over which each group of
# individuals is spread (the components with a mean ancestry above 0.1 in
# that group). Results from the development container (5 replicates of each
# scenario):
# - treated as diploid: AIC keeps decreasing up to the largest K tested,
#   K = 6, in every replicate of both scenarios (cross-validation picks
#   K = 2). Where the pseudohaploids share population A with diploids
#   (scenario 1), the ghost clusters split population A, diploid and
#   pseudohaploid individuals alike, into arbitrary pieces; where a
#   population is entirely pseudohaploid (scenario 2), the ghost clusters
#   all fall within it, splitting it into 4 or 5 clusters, whereas the diploid
#   population A remains a single cluster;
# - pseudohaploid-aware: both AIC and cross-validation pick K = 2 in every
#   replicate, with one cluster per population.
#
# The results are written to data-raw/openadmixture/pseudohaploid_ghost_k.csv.
#
# Run from the root of the package, with tidygenclust installed (this takes
# a few minutes):
#   Rscript data-raw/openadmixture/pseudohaploid_ghost_k.R

library(tidygenclust)
library(tidypopgen)

n_loci <- 5000
fst <- 0.1
# number of diploid and pseudohaploid individuals from each population
scenarios <- list(
  "pseudohaploids within A" = c(
    n_dip_a = 20, n_dip_b = 20, n_pseudo_a = 40, n_pseudo_b = 0
  ),
  "B all pseudohaploid" = c(
    n_dip_a = 20, n_dip_b = 0, n_pseudo_a = 0, n_pseudo_b = 40
  )
)
n_replicates <- 5
k_values <- 1:6
n_cores <- min(4L, parallel::detectCores())

# Simulate one dataset, returning a gen_tibble with a population column
simulate_dataset <- function(seed, n) {
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
  diploid <- function(n, p) {
    if (n == 0) return(NULL)
    t(replicate(n, stats::rbinom(length(p), 2, p)))
  }
  # pseudohaploid genotypes: a single allele per locus, stored as 0/2 as
  # tidypopgen does
  pseudohaploid <- function(n, p) {
    if (n == 0) return(NULL)
    t(replicate(n, 2 * stats::rbinom(length(p), 1, p)))
  }
  # (pseudohaploids are drawn first, as in pseudohaploid_vs_diploid.R)
  pseudo_a <- pseudohaploid(n[["n_pseudo_a"]], p_a)
  pseudo_b <- pseudohaploid(n[["n_pseudo_b"]], p_b)
  genotypes <- rbind(
    diploid(n[["n_dip_a"]], p_a), diploid(n[["n_dip_b"]], p_b),
    pseudo_a, pseudo_b
  )
  population <- c(
    rep("A", n[["n_dip_a"]]), rep("B", n[["n_dip_b"]]),
    rep("A (pseudohaploid)", n[["n_pseudo_a"]]),
    rep("B (pseudohaploid)", n[["n_pseudo_b"]])
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
    backingfile = tempfile("ghost_"), quiet = TRUE
  )
}

groups <- c("A", "B", "A (pseudohaploid)", "B (pseudohaploid)")
results <- list()
for (scenario in names(scenarios)) {
  n <- scenarios[[scenario]]
  for (replicate in seq_len(n_replicates)) {
    gt <- simulate_dataset(seed = replicate, n = n)
    gt_pseudo <- gt_pseudohaploid(gt)
    stopifnot(
      sum(indiv_ploidy(gt_pseudo) == 1) == n[["n_pseudo_a"]] + n[["n_pseudo_b"]]
    )

    for (treatment in c("treated as diploid", "pseudohaploid-aware")) {
      x <- if (treatment == "treated as diploid") gt else gt_pseudo
      adm <- gt_openadmixture(x,
        k = k_values, seed = replicate, crossval = TRUE,
        n_cores = n_cores
      )
      n_par <- adm$k * n_loci + (adm$k - 1) * nrow(gt)
      aic <- -2 * adm$loglik + 2 * n_par
      best_k_aic <- adm$k[which.min(aic)]
      row <- data.frame(
        scenario = scenario,
        replicate = replicate,
        treatment = treatment,
        best_k_aic = best_k_aic,
        best_k_cv = adm$k[which.min(adm$cv)]
      )
      # at the K picked by AIC, the number of clusters over which each group
      # of individuals is spread (NA for groups absent from the scenario)
      q <- unclass(get_q_matrix(adm, k = best_k_aic, run = 1))
      for (group in groups) {
        in_group <- gt$population == group
        row[[paste0("n_clusters_", make.names(group))]] <- if (any(in_group)) {
          sum(colMeans(q[in_group, , drop = FALSE]) > 0.1)
        } else {
          NA
        }
      }
      # AIC of each K, relative to the best K
      for (i in seq_along(adm$k)) {
        row[[paste0("delta_aic_k", adm$k[i])]] <- round(aic[i] - min(aic), 1)
      }
      # cross-validation error of each K
      for (i in seq_along(adm$k)) {
        row[[paste0("cv_k", adm$k[i])]] <- round(adm$cv[i], 4)
      }
      print(row, row.names = FALSE)
      results[[length(results) + 1]] <- row
    }
  }
}

results <- do.call(rbind, results)
utils::write.csv(results,
  file.path("data-raw", "openadmixture", "pseudohaploid_ghost_k.csv"),
  row.names = FALSE
)

# summary: how often each treatment picks each K, by each criterion
for (scenario in names(scenarios)) {
  cat("\nScenario:", scenario, "\n")
  res <- results[results$scenario == scenario, ]
  print(table(treatment = res$treatment, best_k_aic = res$best_k_aic))
  print(table(treatment = res$treatment, best_k_cv = res$best_k_cv))
}

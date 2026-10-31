# Analyses of simulated datasets with pseudohaploids, treated as diploids
# (undeclared) and declared with gt_pseudohaploid(). The results are shown in
# the article pseudohaploids.Rmd (in the parent directory of this script),
# which also shows the code of the chunks of this script (the sections that
# start with a "## ---- label" line, read with knitr::read_chunk()).
#
# Run from the directory of the article (vignettes/articles), with
# tidygenclust installed with R CMD INSTALL (devtools::load_all() compiles
# its C++ code without optimisation):
#
#   cd vignettes/articles
#   Rscript code/pseudohaploid_runs.R
#
# It takes about half an hour on 4 cores, and saves its results in
# data/pseudohaploids, which the article reads.

## ---- runs-setup ----
library(tidygenclust)
library(tidypopgen)
source("code/pseudohaploid_simulation.R")
start_time <- Sys.time()
out_dir <- "data/pseudohaploids"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# each replicate simulates a dataset with its own seed (also used for the
# fits), so every dataset can be simulated again on its own
replicates <- 1:5
n_cores <- min(4L, parallel::detectCores())

## ---- runs-analyse ----
# Analyse a dataset for each value of k, with cross-validation, twice: with
# the pseudohaploids undeclared (treated as diploids) and declared with
# gt_pseudohaploid(). Returns the log-likelihood and cross-validation error
# of each fit (fits), and, if with_q is TRUE, the ancestry proportions of
# each individual for each k > 1 (q, in long format)
analyse <- function(gt, k, seed, with_q = FALSE) {
  fits <- list()
  q <- list()
  for (treatment in c("undeclared", "declared")) {
    x <- if (treatment == "declared") gt_pseudohaploid(gt) else gt
    adm <- gt_radmix(
      x,
      k = k,
      seed = seed,
      crossval = TRUE,
      n_cores = n_cores
    )
    fits[[treatment]] <- data.frame(
      treatment = treatment,
      k = adm$k,
      loglik = adm$loglik,
      cv = adm$cv
    )
    if (with_q) {
      for (i in which(adm$k > 1)) {
        q_i <- unclass(adm$Q[[i]])
        q[[paste(treatment, i)]] <- data.frame(
          treatment = treatment,
          k = adm$k[i],
          id = gt$id,
          population = gt$population,
          type = gt$type,
          cluster = rep(seq_len(ncol(q_i)), each = nrow(q_i)),
          q = round(as.vector(q_i), 4)
        )
      }
    }
  }
  list(fits = do.call(rbind, fits), q = do.call(rbind, q))
}

## ---- runs-scenarios ----
# three scenarios, all with 5,000 loci and two populations (so the true K is
# 2): pseudohaploids in one of the two populations, with moderate (F_ST =
# 0.1) or weak (F_ST = 0.05) structure; and a population known only from
# pseudohaploids (e.g. ancient samples), next to a diploid population
scenarios <- list(
  pseudohaploids_in_a = list(
    n_diploid = c(A = 20, B = 20), n_pseudohaploid = c(A = 40), fst = 0.1
  ),
  weak_structure = list(
    n_diploid = c(A = 20, B = 20), n_pseudohaploid = c(A = 40), fst = 0.05
  ),
  b_pseudohaploid = list(
    n_diploid = c(A = 20), n_pseudohaploid = c(B = 40), fst = 0.1
  )
)
# each scenario is analysed for K = 1 to 6
fits <- list()
q <- list()
for (scenario in names(scenarios)) {
  for (replicate in replicates) {
    message("Scenario ", scenario, ", replicate ", replicate)
    s <- scenarios[[scenario]]
    gt <- simulate_pseudohaploids(
      n_diploid = s$n_diploid,
      n_pseudohaploid = s$n_pseudohaploid,
      fst = s$fst,
      seed = replicate
    )
    res <- analyse(gt, k = 1:6, seed = replicate, with_q = TRUE)
    run <- paste(scenario, replicate)
    fits[[run]] <- data.frame(scenario, replicate, res$fits)
    q[[run]] <- data.frame(scenario, replicate, res$q)
  }
}
fits <- do.call(rbind, fits)
q <- do.call(rbind, q)

## ---- runs-grid ----
# when does cross-validation prefer K = 2 over K = 1? 20 diploids from each
# population, and 0 to 80 pseudohaploids from population A, with weak to
# moderate structure (K = 1 and 2 only)
grid <- expand.grid(
  fst = c(0.03, 0.05, 0.1),
  n_pseudohaploid = c(0, 20, 40, 80),
  replicate = replicates
)
grid_fits <- list()
for (i in seq_len(nrow(grid))) {
  message("Grid ", i, " of ", nrow(grid))
  gt <- simulate_pseudohaploids(
    n_diploid = c(A = 20, B = 20),
    n_pseudohaploid = c(A = grid$n_pseudohaploid[i]),
    fst = grid$fst[i],
    seed = grid$replicate[i]
  )
  res <- analyse(gt, k = 1:2, seed = grid$replicate[i])
  grid_fits[[i]] <- data.frame(grid[i, ], res$fits, row.names = NULL)
}
grid_fits <- do.call(rbind, grid_fits)

## ---- runs-sample-size ----
# why do the individuals of the larger population appear slightly admixed
# at K = 2? Two populations of different or equal sizes, with or without
# pseudohaploids, fitted at K = 2 with the allele frequencies estimated (with
# the pseudohaploids undeclared and declared), and with the allele
# frequencies fixed at their true values (pseudohaploids declared)
sample_size_designs <- list(
  equal = list(n_diploid = c(A = 20, B = 20), n_pseudohaploid = c()),
  a_larger = list(n_diploid = c(A = 60, B = 20), n_pseudohaploid = c()),
  pseudohaploids_in_a = list(
    n_diploid = c(A = 20, B = 20), n_pseudohaploid = c(A = 40)
  ),
  b_larger = list(n_diploid = c(A = 20, B = 60), n_pseudohaploid = c(A = 40))
)
# more replicates than above, as the differences are small
sample_size_replicates <- 1:10
sample_size <- list()
for (design in names(sample_size_designs)) {
  for (replicate in sample_size_replicates) {
    message("Design ", design, ", replicate ", replicate)
    d <- sample_size_designs[[design]]
    gt <- simulate_pseudohaploids(
      n_diploid = d$n_diploid,
      n_pseudohaploid = d$n_pseudohaploid,
      fst = 0.1,
      seed = replicate
    )
    fit <- function(x, frequencies) {
      adm <- gt_radmix(
        x,
        k = 2,
        seed = replicate,
        n_cores = n_cores,
        init_p = if (frequencies == "true") attr(gt, "frequencies"),
        fix_p = frequencies == "true"
      )
      q_i <- unclass(adm$Q[[1]])
      # the cluster of each population: with the true frequencies, the
      # column of the population; otherwise, the cluster with the largest
      # mean ancestry among its individuals
      own <- if (frequencies == "true") {
        c(A = 1, B = 2)
      } else {
        sapply(c(A = "A", B = "B"), function(pop) {
          unname(which.max(colMeans(q_i[gt$population == pop, ])))
        })
      }
      own <- own[gt$population]
      # the mean ancestry of each group from the other cluster
      other <- 1 - q_i[cbind(seq_len(nrow(q_i)), own)]
      data.frame(
        design = design,
        replicate = replicate,
        frequencies = frequencies,
        stats::aggregate(
          list(other = other),
          list(population = gt$population, type = gt$type),
          mean
        )
      )
    }
    runs <- list(
      fit(gt_pseudohaploid(gt), "estimated"),
      fit(gt_pseudohaploid(gt), "true")
    )
    runs[[1]]$treatment <- runs[[2]]$treatment <- "declared"
    if (length(d$n_pseudohaploid) > 0) {
      runs[[3]] <- fit(gt, "estimated")
      runs[[3]]$treatment <- "undeclared"
    }
    sample_size[[paste(design, replicate)]] <- do.call(rbind, runs)
  }
}
sample_size <- do.call(rbind, sample_size)

## ---- runs-save ----
write.csv(fits, file.path(out_dir, "fits.csv"), row.names = FALSE)
write.csv(grid_fits, file.path(out_dir, "grid.csv"), row.names = FALSE)
write.csv(
  sample_size,
  file.path(out_dir, "sample_size.csv"),
  row.names = FALSE
)
q_file <- gzfile(file.path(out_dir, "q.csv.gz"), "w")
write.csv(q, q_file, row.names = FALSE)
close(q_file)
session <- data.frame(
  item = c("date", "R", "tidygenclust", "tidypopgen"),
  value = c(
    format(Sys.Date()),
    R.version.string,
    as.character(utils::packageVersion("tidygenclust")),
    as.character(utils::packageVersion("tidypopgen"))
  )
)
write.csv(session, file.path(out_dir, "session.csv"), row.names = FALSE)
message(
  "Results saved in ", out_dir, ", after ",
  format(round(difftime(Sys.time(), start_time, units = "mins"), 1))
)

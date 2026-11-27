# Compare gt_clumppling() with the python implementation of clumppling on
# random scenarios, beyond the fixed cases stored for the tests: random
# subsets of the capeverde replicates (with shuffled individuals and random
# options), and simulated replicates with planted modes. For each scenario,
# the modes, their representatives and statistics, the costs across K and the
# aligned modes must be the same. The only exception is when a mode includes
# replicates with identical Q matrices (some capeverde runs are identical):
# their costs to the other replicates differ only by rounding errors, so the
# two implementations can pick different (but identical) representatives,
# which is reported as "same (tied representatives)".
#
# Requirements: the conda environment with clumppling 2.0.10, created with
#   bash data-raw/clumppling/setup_clumppling_env.sh
# Run from the root of the package (this takes a few minutes):
#   Rscript data-raw/clumppling/validate_vs_python.R
# The results are written to data-raw/clumppling/validate_vs_python.csv.

options(tidygenclust.pkg_root = normalizePath("."))
source("tests/testthat/helper-clumppling.R")
source("data-raw/clumppling/py_clumppling.R")
pkgload::load_all(".", quiet = TRUE)

n_scenarios <- 30
cd_methods <- c(
  "louvain", "markov_clustering", "infomap", "label_propagation",
  "walktrap"
)

# compare two gt_clumppling objects obtained from the gt_admix object x,
# returning the names of the elements that differ ("" if they are the same)
compare_clumppling <- function(res, ref, x, tolerance = 1e-10) {
  # representatives that differ, but have identical Q matrices
  q_data <- clumppling_prepare_q_list(x$Q)
  rep_res <- res$mode_stats$Representative
  rep_ref <- ref$mode_stats$Representative
  tied <- length(rep_res) == length(rep_ref) &&
    any(rep_res != rep_ref) &&
    all(vapply(which(rep_res != rep_ref), function(i) {
      identical(
        q_data$q_mats[[match(rep_res[i], q_data$q_names)]],
        q_data$q_mats[[match(rep_ref[i], q_data$q_names)]]
      )
    }, logical(1)))
  if (tied) {
    res$mode_stats$Representative <- rep_ref
  }
  checks <- list(
    K_range = identical(res$K_range, ref$K_range),
    N = identical(res$N, ref$N),
    mode_replicates = identical(
      lapply(res$mode_replicates, as.integer),
      lapply(ref$mode_replicates, as.integer)
    ),
    mode_names = identical(
      res$mode_stats[c("Mode", "Representative", "Size")],
      ref$mode_stats[c("Mode", "Representative", "Size")]
    ),
    mode_stats = isTRUE(all.equal(res$mode_stats, ref$mode_stats,
      tolerance = tolerance
    )),
    cost_acrossK = isTRUE(all.equal(res$cost_acrossK, ref$cost_acrossK,
      tolerance = tolerance
    )),
    aligned_modes = isTRUE(all.equal(res$aligned_modes, ref$aligned_modes,
      tolerance = tolerance
    ))
  )
  failed <- names(checks)[!unlist(checks)]
  if (tied && length(failed) == 0) {
    return("same (tied representatives)")
  }
  paste(failed, collapse = ";")
}

# simulated replicates: for each K, two modes (different Q matrices), and
# replicates that are noisy copies of one of them with permuted clusters
simulate_replicates <- function(n_ind, k_values, n_reps) {
  q_list <- list()
  for (k in k_values) {
    modes <- list(
      clumppling_test_random_q(n_ind, k, alpha = 0.3),
      clumppling_test_random_q(n_ind, k, alpha = 0.3)
    )
    for (i in seq_len(n_reps)) {
      q <- modes[[sample(2, 1, prob = c(0.7, 0.3))]]
      q <- q[, sample(k), drop = FALSE] +
        matrix(stats::runif(n_ind * k, 0, 0.02), n_ind, k)
      q_list[[length(q_list) + 1]] <- q / rowSums(q)
    }
  }
  clumppling_test_admix(q_list[sample(length(q_list))])
}

capeverde <- clumppling_test_capeverde()
set.seed(20261023)
results <- list()
for (i in seq_len(2 * n_scenarios)) {
  simulated <- i > n_scenarios
  if (simulated) {
    k_values <- sort(sample(2:6, sample(2:4, 1)))
    x <- simulate_replicates(sample(50:200, 1), k_values, sample(3:15, 1))
    data <- "simulated"
  } else {
    # 2 to 4 values of K, with a random number of replicates for each
    k_values <- sort(sample(2:5, sample(2:4, 1)))
    idx <- unlist(lapply(k_values, function(k) {
      runs <- which(capeverde$k == k)
      sample(runs, sample(c(1, 2, 5:50), 1))
    }))
    q_list <- capeverde$Q[sample(idx)]
    if (stats::runif(1) < 0.5) {
      ind_order <- sample(399)
      q_list <- lapply(q_list, function(q) q[ind_order, ])
    }
    x <- clumppling_test_admix(q_list)
    data <- "capeverde"
  }
  args <- list(
    use_rep = stats::runif(1) < 0.5,
    merge = stats::runif(1) < 0.5,
    use_best_pair = stats::runif(1) < 0.5,
    cd_method = sample(cd_methods, 1, prob = c(0.6, 0.1, 0.1, 0.1, 0.1)),
    test_comm = stats::runif(1) < 0.2,
    seed = i
  )
  py_res <- do.call(py_clumppling, c(list(x), args))
  res <- suppressMessages(do.call(gt_clumppling, c(list(x), args)))
  results[[i]] <- data.frame(
    scenario = i,
    data = data,
    k_values = paste(k_values, collapse = ","),
    n_replicates = length(x$Q),
    n_individuals = nrow(x$Q[[1]]),
    use_rep = args$use_rep,
    merge = args$merge,
    use_best_pair = args$use_best_pair,
    cd_method = args$cd_method,
    test_comm = args$test_comm,
    n_modes = length(res$aligned_modes),
    differences = compare_clumppling(res, py_res, x)
  )
  message(
    "Scenario ", i, " (", data, "): ",
    if (results[[i]]$differences == "") "same" else results[[i]]$differences
  )
}
results <- do.call(rbind, results)
same <- results$differences %in% c("", "same (tied representatives)")
utils::write.csv(results, "data-raw/clumppling/validate_vs_python.csv",
  row.names = FALSE
)
message(
  sum(same), " of ", nrow(results), " scenarios give the same results (",
  sum(results$differences != ""), " with tied representatives)"
)

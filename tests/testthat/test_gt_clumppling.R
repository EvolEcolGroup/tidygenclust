# Tests of gt_clumppling(), comparing it with the results of the python
# implementation of clumppling stored in
# fixtures/clumppling_python_reference.rds (see
# data-raw/clumppling/regenerate_python_fixtures.R), and of the methods for
# gt_clumppling objects

ref <- clumppling_test_py_reference()
capeverde <- clumppling_test_capeverde()

# the modes implied by community labels (0-based, one per replicate of each K)
communities_to_modes <- function(communities) {
  modes <- list()
  for (k_label in names(communities)) {
    comm <- as.integer(communities[[k_label]])
    for (m in sort(unique(comm))) {
      modes[[paste0(k_label, "M", m + 1)]] <- which(comm == m)
    }
  }
  modes
}

test_that("gt_clumppling reproduces the python implementation", {
  res <- gt_clumppling(capeverde, seed = 42)
  expect_same_clumppling(res, ref$capeverde_list)
  # the aligned modes are identical to python's up to rounding
  expect_equal(res$aligned_modes, ref$capeverde_list$aligned_modes,
    tolerance = 1e-14
  )
})

test_that("options of gt_clumppling reproduce the python implementation", {
  expect_same_clumppling(
    gt_clumppling(capeverde, use_rep = FALSE, seed = 42),
    ref$capeverde_avg,
    fingerprint = TRUE
  )
  expect_same_clumppling(
    gt_clumppling(capeverde, use_best_pair = FALSE, seed = 42),
    ref$capeverde_major,
    fingerprint = TRUE
  )
  expect_same_clumppling(
    gt_clumppling(capeverde, merge = FALSE, seed = 42),
    ref$capeverde_nomerge,
    fingerprint = TRUE
  )
})

test_that("STRUCTURE results with non-consecutive K match python", {
  # python read the original STRUCTURE files, whose rows (rounded to 3
  # decimals) are rescaled to sum to 1 in the gt_admix object.
  # The replicates with K = 17 and K = 19 do not have a clear community
  # structure, and python's louvain itself gives different modes depending on
  # the seed; with this seed, R's random numbers give the same modes as
  # python with seed 42
  res <- gt_clumppling(clumppling_test_chicken(), seed = 7)
  expect_identical(res$K_range, c(17L, 19L, 21L))
  expect_same_clumppling(res, ref$chicken, fingerprint = TRUE)
})

test_that("K values with a single replicate match the python implementation", {
  anole_admix <- tidypopgen::read_q_files(
    testthat::test_path("..", "data", "anole_q")
  )
  res <- gt_clumppling(anole_admix)
  expect_same_clumppling(res, ref$anole)
  expect_equal(res$aligned_modes, ref$anole$aligned_modes, tolerance = 1e-14)
})

test_that("shuffled individuals and runs with K = 1 match python", {
  shuffled <- clumppling_test_shuffled()
  args <- list(
    default = list(),
    avg_major = list(use_rep = FALSE, use_best_pair = FALSE),
    nomerge = list(merge = FALSE)
  )
  for (i in names(args)) {
    msgs <- capture_messages(
      res <- do.call(gt_clumppling, c(list(shuffled, seed = 1), args[[i]]))
    )
    expect_identical(msgs, c(
      "Ignored Q matrix 1 with K=1.\n", "Ignored Q matrix 42 with K=1.\n"
    ))
    expect_same_clumppling(res, ref$capeverde_shuffled[[i]],
      fingerprint = TRUE
    )
  }
})

test_that("community detection methods give the same modes as python", {
  methods <- c(
    "markov_clustering", "infomap", "label_propagation", "walktrap"
  )
  for (method in methods) {
    if (method != "markov_clustering") {
      skip_if_not_installed("igraph")
    }
    res <- gt_clumppling(capeverde, cd_method = method, seed = 1)
    py_modes <- communities_to_modes(lapply(
      ref$capeverde_within_k,
      function(x) x$communities[[paste0(method, "_False")]]
    ))
    expect_identical(lapply(res$mode_replicates, as.integer), py_modes,
      info = method
    )
  }
  # with the test for community structure (not available in tidygenclust)
  res <- gt_clumppling(capeverde, test_comm = TRUE, seed = 42)
  py_modes <- communities_to_modes(lapply(
    ref$capeverde_within_k,
    function(x) x$communities[["louvain_True"]]
  ))
  expect_identical(res$mode_replicates, py_modes)
  # a custom function, here putting all the replicates in one mode
  res <- gt_clumppling(capeverde, cd_method = function(adj) {
    rep(1, nrow(adj))
  })
  expect_named(res$aligned_modes, c("K2M1", "K3M1", "K4M1", "K5M1"))
  expect_identical(res$mode_stats$Size, c(50, 50, 50, 50))
})

test_that("the seed makes results reproducible, and is restored", {
  q_sub <- clumppling_test_admix(capeverde$Q[c(1:10, 151:160)])
  set.seed(123)
  seed_before <- .Random.seed
  res <- gt_clumppling(q_sub, seed = 1)
  expect_identical(.Random.seed, seed_before)
  expect_identical(gt_clumppling(q_sub, seed = 1), res)
  withr::with_preserve_seed({
    rm(".Random.seed", envir = globalenv())
    gt_clumppling(q_sub, seed = 1)
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  })
})

test_that("gt_clumppling checks its inputs", {
  expect_error(
    gt_clumppling(clumppling_test_admix(capeverde$Q[1])),
    "'x' must contain more than one Q matrix"
  )
  expect_error(
    # Q matrices are ordered by run and then by K, so these are two runs
    # with the same K
    gt_clumppling(clumppling_test_admix(capeverde$Q[c(1, 5)])),
    "At least two K values are required for alignment across K: K=2"
  )
  # Q files are read by tidypopgen, and lists of Q matrices are not accepted
  expect_error(gt_clumppling(42), "'x' must be a gt_admix object")
  expect_error(gt_clumppling(capeverde$Q), "'x' must be a gt_admix object")
  expect_error(
    gt_clumppling(testthat::test_path("..", "data", "anole_q")),
    "tidypopgen::read_q_files"
  )
  expect_error(
    gt_clumppling(capeverde, cd_method = "kmeans"),
    "'cd_method' must be one of"
  )
  expect_warning(
    gt_clumppling(
      clumppling_test_admix(capeverde$Q[c(1, 51)]),
      output_path = tempfile()
    ),
    "'output_path' is deprecated"
  )
})

# the methods for gt_clumppling objects
clump_res <- gt_clumppling(capeverde, seed = 42)

test_that("gt_clumppling works as expected", {
  expect_true(inherits(clump_res, "gt_clumppling"))
  # known number of indivs in capeverde data
  expect_true(clump_res$N == 399)
  # known k values in capeverde data
  expect_true(all(clump_res$K_range == c(2, 3, 4, 5)))
  expect_true(all(c(
    "N", "K_range", "mode_replicates", "cost_acrossK",
    "aligned_modes", "mode_stats"
  ) %in%
    names(clump_res)))

  get_modes <- tidy(clump_res, matrix = "modes")
  expect_identical(get_modes$label, names(clump_res$aligned_modes))
  get_major_modes <- tidy(clump_res, matrix = "major_modes")
  expect_identical(nrow(get_major_modes), 4L)
  get_q_modes <- tidy(clump_res, matrix = "q_modes")
  expect_identical(nrow(get_q_modes$K5M1), 399L * 5L)
  get_q_major_modes <- tidy(clump_res, matrix = "q_major_modes")
  expect_length(get_q_major_modes, 4)

  # Plotting the results
  expect_s3_class(autoplot(clump_res, type = "modes_within_k", k = 3), "ggplot")
  expect_s3_class(autoplot(clump_res, type = "major_modes"), "ggplot")
  expect_s3_class(autoplot(clump_res, type = "all_modes"), "ggplot")
  expect_s3_class(autoplot(clump_res, type = "modes"), "ggplot")
  # with the populations stored in the gt_admix object
  expect_s3_class(
    autoplot(clump_res, type = "all_modes", group = capeverde$group),
    "ggplot"
  )
})

test_that("subsetting gt_clumppling", {
  sub_clump_res <- subset_gt_clumppling(clump_res, indivs = 1:100)

  expect_warning(
    autoplot(sub_clump_res, type = "modes_within_k", k = 3),
    "modes will not be representative of this subset"
  )
  expect_warning(
    autoplot(sub_clump_res, type = "major_modes"),
    "modes will not be representative of this subset"
  )
  expect_warning(
    autoplot(sub_clump_res, type = "all_modes"),
    "modes will not be representative of this subset"
  )
  expect_warning(
    autoplot(sub_clump_res, type = "modes"),
    "modes will not be representative of this subset"
  )

  sub_clump_res <- subset_gt_clumppling(clump_res, k = 2:4)

  autoplot(sub_clump_res, type = "modes_within_k", k = 3)
  autoplot(sub_clump_res, type = "major_modes")
  autoplot(sub_clump_res, type = "all_modes")
  autoplot(sub_clump_res, type = "modes")

  sub_clump_res <- subset_gt_clumppling(clump_res, indivs = 1:100, k = 2:4)

  expect_warning(
    autoplot(sub_clump_res, type = "modes_within_k", k = 3),
    "modes will not be representative of this subset"
  )
  expect_warning(
    autoplot(sub_clump_res, type = "major_modes"),
    "modes will not be representative of this subset"
  )
  expect_warning(
    autoplot(sub_clump_res, type = "all_modes"),
    "modes will not be representative of this subset"
  )
  expect_warning(
    autoplot(sub_clump_res, type = "modes"),
    "modes will not be representative of this subset"
  )

  expect_error(
    subset_gt_clumppling(clump_res, k = c(2, 4)),
    "k values must be consecutive from min to max"
  )

  expect_error(
    subset_gt_clumppling(clump_res, k = c(4, 3, 2)),
    "k values must be consecutive from min to max"
  )

  expect_error(
    subset_gt_clumppling(clump_res, indivs = c(1, 3, 5)),
    "indivs values must be strictly increasing and consecutive"
  )

  expect_error(
    subset_gt_clumppling(clump_res, indivs = c(5, 4, 3)),
    "indivs values must be strictly increasing and consecutive"
  )
})

test_that("subsetting gt_clumppling by k also subsets mode_stats", {
  sub_clump_res <- subset_gt_clumppling(clump_res, k = 2:4)
  expect_identical(sub_clump_res$mode_stats$Mode,
                   names(sub_clump_res$aligned_modes))
})

test_that("clusters keep their colour and position across K, for K >= 10", {
  # a single individual whose membership to cluster c is proportional to c,
  # so that each segment of the bar identifies its cluster by its height
  cluster_fills <- function(k) {
    q_tidied <- tibble::tibble(
      id = 1L, q = as.character(seq_len(k)),
      percentage = seq_len(k) / sum(seq_len(k))
    )
    bar <- ggplot2::layer_data(plot_q(q_tidied, "K", NULL))
    bar <- bar[order(bar$ymin), ]
    cluster <- round((bar$ymax - bar$ymin) * sum(seq_len(k)))
    # clusters are stacked in order, the first one at the bottom
    expect_identical(as.integer(cluster), seq_len(k))
    bar$fill
  }
  fill_11 <- cluster_fills(11)
  fill_12 <- cluster_fills(12)
  expect_identical(fill_12[1:11], fill_11)
  expect_identical(cluster_fills(3), fill_11[1:3])
})

test_that("the options of gt_clumppling are validated", {
  # a non-finite resolution used to make louvain loop forever
  for (bad in list(NA, NaN, Inf, -1, 0, c(1, 2), "1")) {
    expect_error(gt_clumppling(capeverde, cd_res = bad), "'cd_res' must be")
  }
  expect_error(gt_clumppling(capeverde, merge = "no"), "'merge' must be")
  expect_error(gt_clumppling(capeverde, use_rep = NA), "'use_rep' must be")
  expect_error(
    gt_clumppling(capeverde, use_best_pair = c(TRUE, FALSE)),
    "'use_best_pair' must be"
  )
  expect_error(gt_clumppling(capeverde, test_comm = 1), "'test_comm' must be")
  expect_error(
    gt_clumppling(capeverde, comm_min = NA),
    "'comm_min' and 'comm_max'"
  )
  expect_error(
    gt_clumppling(capeverde, comm_max = -1),
    "'comm_min' and 'comm_max'"
  )
})

# Tests of the arguments of gt_adamixture(), which are checked before
# python is started, and of what is passed to python and returned, with a
# fake python wrapper, so that these tests do not need ADAMIXTURE

# a fake python wrapper, which records the arguments of each fit and returns
# Q, P and the cross-validation error of the right dimensions
local_fake_adamixture <- function(env = parent.frame()) {
  fits <- new.env()
  fits$args <- list()
  fake_module <- list(adamixture_fit = function(...) {
    args <- list(...)
    fits$args[[length(fits$args) + 1]] <- args
    n <- if (args$plink) 3 else length(args$rows)
    m <- if (args$plink) 2 else length(args$cols)
    list(
      Q = matrix(1 / args$K, n, args$K),
      P = matrix(0.5, m, args$K),
      cv = if (is.null(args$cv)) NULL else 0.1
    )
  })
  testthat::local_mocked_bindings(
    .ensure_adamixture_python = function() fake_module,
    .env = env
  )
  fits
}

test_that("invalid arguments are rejected before starting python", {
  local_mocked_bindings(
    .ensure_adamixture_python = function() stop("python started")
  )
  expect_error(gt_adamixture("plink_prefix", k = 0), "'k' should be")
  expect_error(gt_adamixture("plink_prefix", k = 2.5), "'k' should be")
  expect_error(gt_adamixture("plink_prefix", k = 2, n_runs = 0), "'n_runs'")
  expect_error(gt_adamixture("plink_prefix", k = 2, threads = NA), "'threads'")
  expect_error(
    gt_adamixture("plink_prefix", k = 2, max_iter = 0), "'max_iter'"
  )
  for (bad_tol in list(-1, 0, NA, c(0.1, 0.2), "a")) {
    expect_error(
      gt_adamixture("plink_prefix", k = 2, tol = bad_tol),
      "'tol' should be a single positive number"
    )
  }
  expect_error(
    gt_adamixture("plink_prefix", k = 2, n_runs = 3, seed = 1:2),
    "'seed' should be a single number or a vector of length 'n_runs'"
  )
  expect_error(
    gt_adamixture("plink_prefix", k = 2, seed = NA),
    "'seed' should only contain numbers"
  )
  for (bad_cv in list(1, 2.5, c(5, 5), NA)) {
    expect_error(
      gt_adamixture("plink_prefix", k = 2, cv = bad_cv),
      "cv must be a single integer greater than 1"
    )
  }
  expect_error(gt_adamixture(42, k = 2), "data must be a gen_tibble")
})

test_that("each fit gets its own value of k and seed", {
  fits <- local_fake_adamixture()
  # the default seed (42) is incremented for each run, and the same seeds
  # are used for each value of k
  adm <- gt_adamixture("plink_prefix", k = 2:3, n_runs = 2)
  expect_identical(
    vapply(fits$args, function(a) a$seed, integer(1)),
    c(42L, 43L, 42L, 43L)
  )
  expect_identical(
    vapply(fits$args, function(a) a$K, integer(1)), c(2L, 2L, 3L, 3L)
  )
  expect_identical(adm$k, c(2L, 2L, 3L, 3L))
  expect_length(adm$Q, 4)
  expect_s3_class(adm$Q[[1]], "q_matrix")
  expect_identical(adm$algorithm, "adamixture")
  # neither P nor cv by default
  expect_false("P" %in% names(adm))
  expect_false("cv" %in% names(adm))
})

test_that("the settings are passed to python as the right types", {
  fits <- local_fake_adamixture()
  adm <- gt_adamixture("plink_prefix",
    k = 2, threads = 4, seed = 7, tol = 0.01, max_iter = 50, cv = 5,
    no_freqs = FALSE, verbose = FALSE
  )
  args <- fits$args[[1]]
  expect_identical(args$bfile, "plink_prefix")
  expect_true(args$plink)
  expect_identical(args$threads, 4L)
  expect_identical(args$seed, 7L)
  expect_identical(args$tol, 0.01)
  expect_identical(args$max_iter, 50L)
  expect_identical(args$cv, 5L)
  expect_false(args$verbose)
  # P and cv are kept when asked for
  expect_length(adm$P, 1)
  expect_identical(adm$cv, 0.1)
  # without cross-validation, cv is passed as NULL (None)
  gt_adamixture("plink_prefix", k = 2)
  expect_null(fits$args[[2]]$cv)
})

test_that("a gen_tibble is read through the indices of its backing file", {
  test_gt <- radmix_test_gt(matrix(c(0, 1, 2, 1, 0, 2, 2, 2, 0), 3, 3))
  fits <- local_fake_adamixture()
  # a subset, with the individuals and loci reordered
  sub_gt <- tidypopgen::select_loci(test_gt[c(3, 1), ], c(3, 1))
  adm <- gt_adamixture(sub_gt, k = 2, no_freqs = FALSE)
  args <- fits$args[[1]]
  expect_false(args$plink)
  expect_identical(args$bfile, tidypopgen::gt_get_file_names(test_gt)[2])
  expect_identical(args$fbm_dim, c(3L, 3L))
  expect_identical(args$rows, c(2L, 0L))
  expect_identical(args$cols, c(2L, 0L))
  expect_identical(adm$id, sub_gt$id)
  expect_identical(dim(adm$Q[[1]]), c(2L, 2L))
  expect_identical(dim(adm$P[[1]]), c(2L, 2L))
})

test_that("the grouping of a gen_tibble is checked before starting python", {
  test_gt <- radmix_test_gt(matrix(c(0, 1, 2, 1, 0, 2), 3, 2))
  grouped_2 <- test_gt %>%
    dplyr::mutate(pop2 = .data$population) %>%
    dplyr::group_by(.data$population, .data$pop2)
  # the fits are not run
  local_mocked_bindings(
    .ensure_adamixture_python = function() stop("python started")
  )
  expect_error(
    gt_adamixture(grouped_2, k = 2),
    "grouped by a single variable \\(it is grouped by population, pop2\\)"
  )
  # a single grouping variable is stored with the ids of the individuals
  local_fake_adamixture()
  adm <- gt_adamixture(dplyr::group_by(test_gt, .data$population), k = 2)
  expect_identical(adm$id, test_gt$id)
  expect_identical(adm$group, test_gt$population)
  expect_null(gt_adamixture(test_gt, k = 2)$group)
})

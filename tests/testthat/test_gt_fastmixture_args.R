# Tests of the arguments of gt_fastmixture(), which are checked before
# python is started, so that these tests do not need fastmixture

test_that("invalid arguments are rejected before starting python", {
  local_mocked_bindings(
    .ensure_fastmixture_python = function() stop("python started")
  )
  expect_error(gt_fastmixture("plink_prefix", k = 1), "'k' should be")
  expect_error(gt_fastmixture("plink_prefix", k = 2.5), "'k' should be")
  expect_error(gt_fastmixture("plink_prefix", k = 2, n_runs = 0), "'n_runs'")
  expect_error(gt_fastmixture("plink_prefix", k = 2, threads = NA), "'threads'")
  expect_error(
    gt_fastmixture("plink_prefix", k = 2, n_runs = 3, seed = 1:2),
    "'seed' should be a single number or a vector of length 'n_runs'"
  )
  expect_error(
    gt_fastmixture("plink_prefix", k = 2, seed = NA),
    "'seed' should only contain numbers"
  )
  for (bad_cv in list(NA, 1, 2.5)) {
    expect_error(
      gt_fastmixture("plink_prefix", k = 2, cv = bad_cv),
      "cv must be a single integer greater than 1"
    )
  }
  expect_error(gt_fastmixture(42, k = 2), "data must be a gen_tibble")
})

test_that("a single seed is incremented for each run", {
  # a fake python wrapper, which records the seed of each run
  seeds <- c()
  fake_module <- list(fastmixture_run = function(args) {
    seeds <<- c(seeds, args$seed)
    matrix(1 / args$K, 3, args$K)
  })
  local_mocked_bindings(.ensure_fastmixture_python = function() fake_module)
  local_mocked_bindings(
    import = function(module, ...) list(Namespace = function(...) list(...)),
    .package = "reticulate"
  )
  # the default seed (42) used to fail with n_runs > 1
  adm <- gt_fastmixture("plink_prefix", k = 2:3, n_runs = 2)
  expect_identical(seeds, c(42L, 43L, 42L, 43L))
  expect_identical(adm$k, c(2L, 2L, 3L, 3L))
})

test_that("'safety' is deprecated and not passed to python", {
  # a fake python wrapper, which records the arguments of each run
  run_args <- list()
  fake_module <- list(fastmixture_run = function(args) {
    run_args[[length(run_args) + 1]] <<- args
    matrix(1 / args$K, 3, args$K)
  })
  local_mocked_bindings(.ensure_fastmixture_python = function() fake_module)
  local_mocked_bindings(
    import = function(module, ...) list(Namespace = function(...) list(...)),
    .package = "reticulate"
  )
  expect_warning(
    gt_fastmixture("plink_prefix", k = 2, safety = TRUE),
    "'safety' is deprecated and ignored"
  )
  expect_no_warning(gt_fastmixture("plink_prefix", k = 2))
  expect_false(any(vapply(run_args, function(a) "safety" %in% names(a), NA)))
})

test_that("the SVD-based initialisation is the default", {
  # a fake python wrapper, which records the arguments of each run
  run_args <- list()
  fake_module <- list(fastmixture_run = function(args) {
    run_args[[length(run_args) + 1]] <<- args
    matrix(1 / args$K, 3, args$K)
  })
  local_mocked_bindings(.ensure_fastmixture_python = function() fake_module)
  local_mocked_bindings(
    import = function(module, ...) list(Namespace = function(...) list(...)),
    .package = "reticulate"
  )
  gt_fastmixture("plink_prefix", k = 2)
  gt_fastmixture("plink_prefix", k = 2, random_init = TRUE)
  expect_identical(
    vapply(run_args, function(a) a$random_init, NA), c(FALSE, TRUE)
  )
})

test_that("the grouping of a gen_tibble is checked before starting python", {
  test_gt <- radmix_test_gt(matrix(c(0, 1, 2, 1, 0, 2), 3, 2))
  grouped_2 <- test_gt %>%
    dplyr::mutate(pop2 = .data$population) %>%
    dplyr::group_by(.data$population, .data$pop2)
  # the fits are not run (they used to run before the error)
  local_mocked_bindings(
    .ensure_fastmixture_python = function() stop("python started")
  )
  expect_error(
    gt_fastmixture(grouped_2, k = 2),
    "grouped by a single variable \\(it is grouped by population, pop2\\)"
  )
  # a single grouping variable is stored with the ids of the individuals
  fake_module <- list(fastmixture_run = function(args) {
    matrix(1 / args$K, length(args$rows), args$K)
  })
  local_mocked_bindings(.ensure_fastmixture_python = function() fake_module)
  local_mocked_bindings(
    import = function(module, ...) list(Namespace = function(...) list(...)),
    .package = "reticulate"
  )
  adm <- gt_fastmixture(dplyr::group_by(test_gt, .data$population), k = 2)
  expect_identical(adm$id, test_gt$id)
  expect_identical(adm$group, test_gt$population)
  expect_null(gt_fastmixture(test_gt, k = 2)$group)
})

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

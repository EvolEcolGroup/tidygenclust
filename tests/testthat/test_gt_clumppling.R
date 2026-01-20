skip_on_cran()

# skip if the conda environment does not exist
skip_if(!reticulate::condaenv_exists("cclumppling"))

input_path <- system.file("extdata/capeverde.zip", package = "tidygenclust")
clump_res <- input_path %>% gt_clumppling()
# TODO we need to check that what we return here is in the form we expect

test_that("gt_clumppling works as expected", {
  expect_true(inherits(clump_res, "gt_clumppling"))
  # known number of indivs in capeverde data
  expect_true(clump_res$N == 399)
  # known k values in capeverde data
  expect_true(all(clump_res$K_range == c(2,3,4,5)))
  expect_true(all(c("N", "K_range", "mode_replicates", "cost_acrossK",
                    "aligned_modes") %in%
                    names(clump_res)))

  get_modes <- tidy(clump_res, matrix = "modes")
  get_major_modes <- tidy(clump_res, matrix = "major_modes")
  get_q_modes <- tidy(clump_res, matrix = "q_modes")
  get_q_major_modes <- tidy(clump_res, matrix = "q_major_modes")

  # Plotting the results
  autoplot(clump_res, type = "modes_within_k", k = 3)
  autoplot(clump_res, type = "major_modes")
  autoplot(clump_res, type = "all_modes")
  autoplot(clump_res, type = "modes")
})


test_that("subsetting gt_clumppling",{
  sub_clump_res <- subset_gt_clumppling(clump_res, indivs = 1:100)

  expect_warning(autoplot(sub_clump_res, type = "modes_within_k", k = 3),
                 "modes will not be representative of this subset")
  expect_warning(autoplot(sub_clump_res, type = "major_modes"),
                 "modes will not be representative of this subset")
  expect_warning(autoplot(sub_clump_res, type = "all_modes"),
                 "modes will not be representative of this subset")
  expect_warning(autoplot(sub_clump_res, type = "modes"),
                 "modes will not be representative of this subset")

  sub_clump_res <- subset_gt_clumppling(clump_res, k = 2:4)

  autoplot(sub_clump_res, type = "modes_within_k", k = 3)
  autoplot(sub_clump_res, type = "major_modes")
  autoplot(sub_clump_res, type = "all_modes")
  autoplot(sub_clump_res, type = "modes")

  sub_clump_res <- subset_gt_clumppling(clump_res, indivs = 1:100, k = 2:4)

  expect_warning(autoplot(sub_clump_res, type = "modes_within_k", k = 3),
                 "modes will not be representative of this subset")
  expect_warning(autoplot(sub_clump_res, type = "major_modes"),
                 "modes will not be representative of this subset")
  expect_warning(autoplot(sub_clump_res, type = "all_modes"),
                 "modes will not be representative of this subset")
  expect_warning(autoplot(sub_clump_res, type = "modes"),
                 "modes will not be representative of this subset")

  expect_error(subset_gt_clumppling(clump_res, k = c(2,4)),
               "k values must be consecutive from min to max")

  expect_error(subset_gt_clumppling(clump_res, k = c(4,3,2)),
               "k values must be consecutive from min to max")

  expect_error(subset_gt_clumppling(clump_res, indivs = c(1,3,5)),
               "indivs values must be strictly increasing and consecutive")

  expect_error(subset_gt_clumppling(clump_res, indivs = c(5,4,3)),
               "indivs values must be strictly increasing and consecutive")
})


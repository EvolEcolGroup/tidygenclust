skip_on_cran()
# ADAMIXTURE is an optional install
skip_if(!tidygenclust:::.tgc_tool_available("adamixture"))

library(tidypopgen)

# set up the gentibble
vcf_path <- system.file(
  "/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
  package = "tidypopgen"
)
anole_gt <- gen_tibble(
  vcf_path,
  quiet = TRUE,
  backingfile = tempfile("anolis_"),
  parser = "cpp"
)
pops_path <- system.file(
  "/extdata/anolis/plot_order_punctatus_n46.csv",
  package = "tidypopgen"
)
pops <- readr::read_csv(pops_path, show_col_types = FALSE)
anole_gt <- anole_gt %>% mutate(id = gsub("punc_", "", .data$id, ))
anole_gt <- anole_gt %>%
  mutate(population = pops$pop[match(pops$ID, .data$id)])


test_that("gt_adamixture", {
  # Multiple k with no P matrices
  k <- c(2:3)
  anole_qmat <- gt_adamixture(anole_gt, k, no_freqs = TRUE)
  expect_true(inherits(anole_qmat, c("list", "gt_admix")))
  expect_true(inherits(anole_qmat$Q[[1]], "q_matrix"))
  expect_identical(anole_qmat$algorithm, "adamixture")
  expect_identical(anole_qmat$id, anole_gt$id)
  expect_equal(anole_qmat$k, k)
  expect_true(!"P" %in% names(anole_qmat))
  expect_true(length(anole_qmat$Q) == length(k))
  expect_true(ncol(anole_qmat$Q[[1]]) == 2)
  expect_true(ncol(anole_qmat$Q[[2]]) == 3)
  expect_equal(nrow(anole_qmat$Q[[1]]), nrow(anole_gt))
  expect_equal(unname(rowSums(anole_qmat$Q[[2]])), rep(1, nrow(anole_gt)),
    tolerance = 1e-6
  )

  # Single k with P matrices
  k <- 2
  anole_qmat <- gt_adamixture(anole_gt, k, no_freqs = FALSE)
  expect_true(length(anole_qmat$P) == length(k))
  expect_true(ncol(anole_qmat$P[[1]]) == 2)
  expect_equal(nrow(anole_qmat$P[[1]]), count_loci(anole_gt))

  # grouped gen_tibble adds group info
  anole_qmat <- gt_adamixture(anole_gt %>% group_by(population), 2)
  expect_identical(anole_qmat$group, anole_gt$population)

  # same seed gives the same result
  q1 <- gt_adamixture(anole_gt, 3, seed = 1)
  q2 <- gt_adamixture(anole_gt, 3, seed = 1)
  expect_equal(q1$Q[[1]], q2$Q[[1]])

  # n_runs creates multiple repeats per k
  anole_runs <- gt_adamixture(anole_gt, 2, n_runs = 2, seed = c(1, 2))
  expect_equal(length(anole_runs$Q), 2)
  expect_equal(anole_runs$k, c(2, 2))
})

test_that("gt_adamixture with cv", {
  # Multiple k with no P matrices
  k <- c(2:4)
  anole_qmat <- gt_adamixture(anole_gt, k, cv = 5)
  expect_true("cv" %in% names(anole_qmat))
  expect_true(is.numeric(anole_qmat$cv))
  expect_equal(length(anole_qmat$cv), length(k))

  # Single k with P matrices
  anole_qmat <- gt_adamixture(anole_gt, 3, cv = 5, no_freqs = FALSE)
  expect_true("cv" %in% names(anole_qmat))
  expect_equal(length(anole_qmat$cv), 1)
  expect_true(ncol(anole_qmat$P[[1]]) == 3)
})

test_that("gt_adamixture output works with tidypopgen methods", {
  anole_qmat <- gt_adamixture(anole_gt %>% group_by(population), 2:3, cv = 5)
  expect_s3_class(
    autoplot(anole_qmat, type = "cv"), "ggplot"
  )
  expect_s3_class(
    autoplot(anole_qmat, type = "barplot", k = 3, run = 1),
    "ggplot"
  )
  tidy_q <- tidy(anole_qmat$Q[[2]], anole_gt)
  expect_equal(nrow(tidy_q), nrow(anole_gt) * 3)
})

test_that("gt_adamixture input checks", {
  expect_error(gt_adamixture(1:3, 2), "data must be a gen_tibble")
  expect_error(gt_adamixture(anole_gt, 2, cv = 1), "cv must be")
})

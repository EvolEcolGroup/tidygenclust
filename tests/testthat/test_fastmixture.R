# skip if the conda environment does not exist
skip_if(!reticulate::condaenv_exists("ctidygenclust"))

test_that("gt_fastmixture", {
  library(tidypopgen)

  #set up the gentibble
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
  pops <- readr::read_csv(pops_path)
  anole_gt <- anole_gt %>% mutate(id = gsub('punc_', "", .data$id, ))
  anole_gt <- anole_gt %>%
    mutate(population = pops$pop[match(pops$ID, .data$id)])

  # Multiple k and one repeat with no P matrices

  k <- c(2:3)

  anole_qmat <- gt_fastmixture(
    anole_gt,
    k,
    n_runs = 1,
    threads = 1,
    seed = 42,
    iter = 1000,
    tole = 0.5,
    batches = 32,
    supervised = NULL,
    check = 5,
    power = 11,
    chunk = 8192,
    als_iter = 1000,
    als_tole = 1e-4,
    no_freqs = TRUE,
    random_init = TRUE
  )
  #correct structure
  expect_true(inherits(anole_qmat, c("list", "gt_admix")))
  expect_true(all(names(anole_qmat) == names(anole_qmat)))
  expect_true(inherits(anole_qmat$Q[[1]], "q_matrix"))
  #correct number of repeats
  expect_true(sum(anole_qmat$k == 2) == 1)
  expect_true(sum(anole_qmat$k == 3) == 1)
  #check no p matrices
  expect_true(!"P" %in% names(anole_qmat))
  #correct number of q matrices
  expect_true(length(anole_qmat$Q) == length(k))
  #check q-matrix correspond to correct k
  index_k2 <- which(anole_qmat$k == 2)
  index_k3 <- which(anole_qmat$k == 3)
  expect_true(ncol(anole_qmat$Q[[index_k2]]) == 2)
  expect_true(ncol(anole_qmat$Q[[index_k3]]) == 3)

  # Multiple repeats and one k with P matrices

  k <- 2

  anole_qmat <- gt_fastmixture(
    anole_gt,
    k,
    n_runs = 2,
    threads = 1,
    seed = c(42, 31),
    iter = 1000,
    tole = 0.5,
    batches = 32,
    supervised = NULL,
    check = 5,
    power = 11,
    chunk = 8192,
    als_iter = 1000,
    als_tole = 1e-4,
    no_freqs = FALSE,
    random_init = TRUE
  )
  #correct structure
  expect_true(inherits(anole_qmat, c("list", "gt_admix")))
  expect_true(all(names(anole_qmat) == names(anole_qmat)))
  expect_true(inherits(anole_qmat$Q[[1]], "q_matrix"))
  #correct number of k
  expect_true(length(anole_qmat$Q) == length(k) * 2)
  #correct number of repeats
  expect_true(sum(anole_qmat$k == 2) == 2)
  #correct number of q matrices
  expect_true(length(anole_qmat$Q) == length(k) * 2)
  #correct number of p matrices
  expect_true(length(anole_qmat$P) == length(k) * 2)
  #check q-matrix are indexed as correct k
  index_k2 <- which(anole_qmat$k == 2)
  expect_true(ncol(anole_qmat$Q[[index_k2[1]]]) == 2)
  expect_true(ncol(anole_qmat$Q[[index_k2[2]]]) == 2)
  #check p-matrix are indexed as correct k
  expect_true(ncol(anole_qmat$P[[index_k2[1]]]) == 2)
  expect_true(ncol(anole_qmat$P[[index_k2[2]]]) == 2)

  # Multiple k and multiple repeats with P matrices

  k <- c(2:3)

  anole_qmat <- gt_fastmixture(
    anole_gt,
    k,
    n_runs = 2,
    threads = 1,
    seed = c(42, 31),
    iter = 1000,
    tole = 0.5,
    batches = 32,
    supervised = NULL,
    check = 5,
    power = 11,
    chunk = 8192,
    als_iter = 1000,
    als_tole = 1e-4,
    no_freqs = FALSE,
    random_init = TRUE
  )
  #correct structure
  expect_true(inherits(anole_qmat, c("list", "gt_admix")))
  expect_true(all(names(anole_qmat) == names(anole_qmat)))
  expect_true(inherits(anole_qmat$Q[[1]], "q_matrix"))
  #correct number of k
  expect_true(length(anole_qmat$Q) == length(k) * 2)
  #correct number of repeats
  expect_true(sum(anole_qmat$k == 2) == 2)
  expect_true(sum(anole_qmat$k == 3) == 2)
  #correct number of q matrices
  expect_true(length(anole_qmat$Q) == length(k) * 2)
  #correct number of p matrices
  expect_true(length(anole_qmat$P) == length(k) * 2)
  #check q-matrix are indexed as correct k
  index_k2 <- which(anole_qmat$k == 2)
  expect_true(ncol(anole_qmat$Q[[index_k2[1]]]) == 2)
  expect_true(ncol(anole_qmat$Q[[index_k2[2]]]) == 2)
  index_k3 <- which(anole_qmat$k == 3)
  expect_true(ncol(anole_qmat$Q[[index_k3[1]]]) == 3)
  expect_true(ncol(anole_qmat$Q[[index_k3[2]]]) == 3)
  #check p-matrix are indexed as correct k
  expect_true(ncol(anole_qmat$P[[index_k2[1]]]) == 2)
  expect_true(ncol(anole_qmat$P[[index_k2[2]]]) == 2)
  expect_true(ncol(anole_qmat$P[[index_k3[1]]]) == 3)
  expect_true(ncol(anole_qmat$P[[index_k3[2]]]) == 3)

  # Single K and one repeat with P matrices

  k <- 2

  anole_qmat <- gt_fastmixture(
    anole_gt,
    k,
    n_runs = 1,
    threads = 1,
    seed = 42,
    iter = 1000,
    tole = 0.5,
    batches = 32,
    supervised = NULL,
    check = 5,
    power = 11,
    chunk = 8192,
    als_iter = 1000,
    als_tole = 1e-4,
    no_freqs = FALSE,
    random_init = TRUE
  )
  #correct structure
  expect_true(inherits(anole_qmat, c("list", "gt_admix")))
  expect_true(all(names(anole_qmat) == names(anole_qmat)))
  expect_true(inherits(anole_qmat$Q[[1]], "q_matrix"))
  #correct number of k
  expect_true(length(anole_qmat$Q) == length(k))
  #correct number of repeats
  expect_true(sum(anole_qmat$k == 2) == 1)
  #correct number of q matrices
  expect_true(length(anole_qmat$Q) == length(k))
  #correct number of p matrices
  expect_true(length(anole_qmat$P) == length(k))
  #check q-matrix are indexed as correct k
  index_k2 <- which(anole_qmat$k == 2)
  expect_true(ncol(anole_qmat$Q[[index_k2[1]]]) == 2)
  #check p-matrix are indexed as correct k
  expect_true(ncol(anole_qmat$P[[index_k2[1]]]) == 2)
})

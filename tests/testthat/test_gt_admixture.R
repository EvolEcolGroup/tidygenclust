skip_on_cran()

# set the input file
vcf_path <-
  system.file(
    "/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
    package = "tidypopgen"
  )
anole_gt <- tidypopgen::gen_tibble(
  vcf_path,
  quiet = TRUE,
  backingfile = tempfile("anolis_")
)
pops_path <- system.file(
  "/extdata/anolis/plot_order_punctatus_n46.csv",
  package = "tidypopgen"
)
pops <- readr::read_csv(pops_path, show_col_types = FALSE)
anole_gt <- anole_gt %>% dplyr::mutate(id = gsub("punc_", "", .data$id, ))
anole_gt <- anole_gt %>%
  dplyr::mutate(population = pops$pop[match(pops$ID, .data$id)])

test_that("run admixture as single run", {
  # run admixture
  anole_adm <- tidypopgen::gt_admixture(
    anole_gt,
    k = 3,
    crossval = FALSE,
    n_cores = 1,
    seed = 123,
    conda_env = "auto"
  )
  # check the output
  expect_true(nrow(anole_adm$Q[[1]]) == nrow(anole_gt))
  expect_true(ncol(anole_adm$Q[[1]]) == 3)
  expect_true(is.null(anole_adm$cv))
})

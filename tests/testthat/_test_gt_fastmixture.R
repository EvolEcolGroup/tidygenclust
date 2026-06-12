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
pops <- readr::read_csv(pops_path)
anole_gt <- anole_gt %>% mutate(id = gsub("punc_", "", .data$id, ))
anole_gt <- anole_gt %>%
  mutate(population = pops$pop[match(pops$ID, .data$id)])

test_that("cannot reshape array error", {
  anole_gt <- gt_impute_simple(anole_gt)
  # subset anole_gt using loci_ld_clump
  id_rm <- loci_ld_clump(
    anole_gt,
    thr_r2 = 0.2,
    return_id = TRUE
  )
  anole_gt_sub <- anole_gt %>%
    select_loci_if(loci_ld_clump(genotypes, thr_r2 = 0.2, size = 10))
  k <- c(3)
  anole_qmat <- gt_fastmixture(
    anole_gt_sub,
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
    random_init = TRUE,
    cv = 7,
    cv_tole = 1e-7
  )

  # Error in py_call_impl(callable, call_args$unnamed, call_args$named) :
  #  ValueError: cannot reshape array of size 149454 into shape (2402,46)
})

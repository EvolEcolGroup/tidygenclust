# Regenerate the reference results of the python implementation of
# clumppling used to test the native port
# (tests/testthat/fixtures/clumppling_python_reference.rds).
#
# Requirements: the conda environment with clumppling 2.0.10, created with
#   bash data-raw/clumppling/setup_clumppling_env.sh
# Run from the root of the package:
#   Rscript data-raw/clumppling/regenerate_python_fixtures.R
# It takes a few minutes, as the python implementation is slow.

options(tidygenclust.pkg_root = normalizePath("."))
source("tests/testthat/helper-clumppling.R")
source("data-raw/clumppling/py_clumppling.R")

drv <- reticulate::import_from_path(
  "clumppling_reference",
  path = "data-raw/clumppling"
)

# replace the aligned modes of a gt_clumppling object with their fingerprints
compact <- function(x) {
  x$aligned_modes <- lapply(x$aligned_modes, clumppling_test_fingerprint)
  x
}

ref <- list()
ref$versions <- c(
  as.list(py_clumppling_versions()),
  date = as.character(Sys.Date())
)

# unit-level comparisons
inputs <- clumppling_test_unit_inputs()
ref$ilp <- drv$ilp_cases(inputs$ilp)
ref$graphs <- drv$graph_cases(inputs$graphs, seeds = list(1L, 2L, 3L))
ref$detect_cases <- drv$detect_communities_cases(list(
  matrix(0),
  matrix(c(0, 5e-7, 5e-7, 0), 2),
  matrix(c(0, 0.5, 0.5, 0), 2),
  (1 - diag(4)) * 5e-7,
  (1 - diag(4)) * 0.5,
  (1 - diag(4)) * 0.005
))

# alignment within K and mode detection on the capeverde data
capeverde <- clumppling_test_capeverde()
ref$capeverde_within_k <- drv$within_k_cases(
  lapply(capeverde$Q, unclass), py_rclumppling$`_prepare_q_list`,
  methods = list(
    "louvain", "infomap", "markov_clustering", "label_propagation",
    "walktrap"
  ),
  seed = 42L
)

# complete runs
ref$capeverde_list <- py_clumppling(capeverde, seed = 42)
ref$capeverde_avg <- compact(py_clumppling(capeverde,
  use_rep = FALSE, seed = 42
))
ref$capeverde_major <- compact(py_clumppling(capeverde,
  use_best_pair = FALSE, seed = 42
))
ref$capeverde_nomerge <- compact(py_clumppling(capeverde,
  merge = FALSE, seed = 42
))
# clumppling reads the original STRUCTURE files (with rows that do not sum to
# 1 because of rounding, which the example gt_admix object rescales)
ref$chicken <- compact(py_clumppling(
  "data-raw/clumppling/chicken_gapk.zip",
  input_format = "structure", extension = "_f", seed = 42
))
ref$anole <- py_clumppling("tests/data/anole_q",
  input_format = "generalQ", extension = ".Q"
)
# shuffled individuals, and runs with K = 1 (which are ignored)
shuffled <- clumppling_test_shuffled()
ref$capeverde_shuffled <- list(
  default = compact(py_clumppling(shuffled, seed = 1)),
  avg_major = compact(py_clumppling(shuffled,
    use_rep = FALSE, use_best_pair = FALSE, seed = 1
  )),
  nomerge = compact(py_clumppling(shuffled, merge = FALSE, seed = 1))
)

saveRDS(ref, "tests/testthat/fixtures/clumppling_python_reference.rds",
  compress = "xz"
)

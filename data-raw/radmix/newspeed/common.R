# Settings and helpers shared by the R scripts of the newspeed kit
# (data-raw/radmix/newspeed), sourced by each of them. Run the scripts from
# the root of the repository. Environment variables:
# - TGC_WORK: working directory of the kit, with the sources and libraries
#   of the builds (build_variants.sh) and the cached gen_tibbles (data/);
#   default ../tgc_newspeed_work
# - TGC_EUR_BED: the EUR_subset bed file (default: the cache of the tests)
# - TGC_HGDP_DIR: the directory of the HGDP bed files (default:
#   vignettes/articles/data/hgdp of the repository)

tgc_repo <- normalizePath(".")
stopifnot(file.exists(file.path(tgc_repo, "DESCRIPTION")))
tgc_work <- Sys.getenv(
  "TGC_WORK", file.path(dirname(tgc_repo), "tgc_newspeed_work")
)

# the bed file of a dataset:
# - eur: EUR_subset, 379 x 54,051, no missing genotypes
# - hgdp_quick: 209 x 25,000, with missing genotypes (toy data: debugging
#   only, never for decisions)
# - hgdp: the HGDP of the benchmark article, 1043 x 130,943 (LD pruned)
# - hgdp_full: the HGDP without LD pruning, 1043 x 642,794
tgc_bed <- function(name) {
  hgdp_dir <- Sys.getenv(
    "TGC_HGDP_DIR",
    file.path(tgc_repo, "vignettes", "articles", "data", "hgdp")
  )
  path <- switch(name,
    eur = Sys.getenv("TGC_EUR_BED", file.path(
      tools::R_user_dir("tidygenclust", which = "cache"), "EUR_subset.bed"
    )),
    hgdp_quick = file.path(hgdp_dir, "hgdp_bench_quick.bed"),
    hgdp = file.path(hgdp_dir, "hgdp_bench.bed"),
    hgdp_full = file.path(hgdp_dir, "hgdp_full.bed"),
    stop("unknown dataset ", name)
  )
  if (!file.exists(path)) {
    stop(path, " not found (see the README of the kit for how to make it)")
  }
  path
}

# load the build of tidygenclust installed in `lib`, check that it is the
# one loaded, and return its namespace (for the internal functions)
tgc_load <- function(lib) {
  suppressMessages({
    library(tidypopgen)
    library(tidygenclust, lib.loc = lib)
  })
  ns <- asNamespace("tidygenclust")
  loaded <- dirname(getNamespaceInfo(ns, "path"))
  stopifnot(normalizePath(loaded) == normalizePath(lib))
  ns
}

# a gen_tibble of a dataset. The first call converts the bed file into a
# backing file under $TGC_WORK/data and saves the gen_tibble there; later
# calls (from any process) reload it, which takes a second instead of up to
# a minute for hgdp_full. The kernels only read the backing file.
tgc_read <- function(name) {
  dir <- file.path(tgc_work, "data")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  gt_file <- file.path(dir, paste0(name, ".gt"))
  if (file.exists(gt_file)) {
    return(tidypopgen::gt_load(gt_file))
  }
  x <- tidypopgen::gen_tibble(tgc_bed(name),
    backingfile = file.path(dir, name), quiet = TRUE
  )
  tidypopgen::gt_save(x, quiet = TRUE)
  x
}

# "eur:3,hgdp:5" -> list(list(data = "eur", k = 3L), ...)
tgc_cases <- function(s) {
  lapply(strsplit(strsplit(s, ",")[[1]], ":"), function(z) {
    list(data = z[1], k = as.integer(z[2]))
  })
}

tgc_ints <- function(s) as.integer(strsplit(s, ",")[[1]])

# the 1-minute load average (Linux; NA elsewhere)
tgc_load_avg <- function() {
  if (!file.exists("/proc/loadavg")) {
    return(NA_real_)
  }
  as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1])
}

# Counts from the log of a fit (radmix_fit_cpp()): the iterations, the
# rejected quasi-Newton points (iterations whose log-likelihood of x2,
# "basic", was computed), and the number of passes over the genotypes of
# the main loop, with and without the fused pass of commit "fewer passes"
# (an iteration makes 4 passes of block updates and 1 log-likelihood of the
# quasi-Newton point, plus 1 for a rejected point; the fused pass saves the
# first update of every iteration but the first; the start costs 1 pass
# either way, and EM priming is not counted).
tgc_log_counts <- function(log) {
  iter_lines <- grep("^Iteration [0-9]+:", log, value = TRUE)
  n_iter <- length(iter_lines)
  n_rej <- sum(!grepl("basic = not computed", iter_lines))
  c(
    iterations = n_iter, rejected = n_rej,
    passes_before = 1 + 5 * n_iter + n_rej,
    passes_fused = 1 + 4 * n_iter + n_rej
  )
}

# append a data frame to a csv file (with a header if the file is new)
tgc_append <- function(df, file) {
  new <- !file.exists(file)
  utils::write.table(df, file,
    sep = ",", append = !new, col.names = new,
    row.names = FALSE
  )
}

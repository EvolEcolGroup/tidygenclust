# Times the two BLAS-bound steps of gt_radmix(init = "svd") with different
# numbers of BLAS threads (RhpcBLASctl::blas_set_num_threads()) and OpenMP
# threads (n_cores): the SVD (radmix_svd_factors(), whose
# gt_pca_partialSVD() gets its threads from the BLAS only) and the start
# (radmix_svd_start(): 3 ALS solutions, whose products call the BLAS and
# whose bounded least squares use OpenMP). See section 2.5 of HANDOVER.md.
#
#   Rscript svd_blas_threads.R <lib> <out dir> <dataset> <blas threads> \
#     <n_cores> [ks] [label]
#
# e.g. Rscript svd_blas_threads.R $TGC_WORK/libs/comb out hgdp 1,2,4,8 1,8
# blas threads: comma separated; 0 = leave the BLAS as it is (its default).
# Appends to <out>/svd_blas_threads.csv (label: a free text to tell runs
# apart, e.g. how the process was pinned).
args <- commandArgs(TRUE)
lib <- args[1]
out <- args[2]
dataset <- args[3]
source("data-raw/radmix/perf_verification/scripts/common.R")
blas <- tgc_ints(args[4])
cores <- tgc_ints(args[5])
ks <- if (length(args) >= 6) tgc_ints(args[6]) else c(5L, 10L)
label <- if (length(args) >= 7) args[7] else ""
ns <- tgc_load(lib)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
default_blas <- RhpcBLASctl::blas_get_num_procs()

x <- tgc_read(dataset)
if (any(tidypopgen::loci_missingness(x) > 0)) {
  x <- tidypopgen::gt_impute_simple(x, method = "mean0")
}
data <- ns$radmix_prepare(x)
rank <- max(ks) + 5
for (b in blas) {
  RhpcBLASctl::blas_set_num_threads(if (b == 0) default_blas else b)
  nb <- RhpcBLASctl::blas_get_num_procs()
  for (n in cores) {
    Sys.sleep(2)
    t_svd <- system.time(
      factors <- ns$radmix_svd_factors(x, rank = rank, n_cores = n)
    )
    row <- data.frame(
      label = label, dataset = dataset, blas_threads = nb, n_cores = n,
      step = "svd", k = NA, seconds = t_svd[["elapsed"]],
      user = t_svd[["user.self"]], system = t_svd[["sys.self"]],
      als_iterations = NA
    )
    tgc_append(row, file.path(out, "svd_blas_threads.csv"))
    print(row, row.names = FALSE)
    for (k in ks) {
      Sys.sleep(2)
      set.seed(1)
      t_start <- system.time(
        st <- ns$radmix_svd_start(data, factors = factors, k = k, n_cores = n)
      )
      iters <- sub(".*: ([0-9]+) iterations.*", "\\1",
        grep("^ALS solution", st$log, value = TRUE))
      row <- data.frame(
        label = label, dataset = dataset, blas_threads = nb, n_cores = n,
        step = "start", k = k, seconds = t_start[["elapsed"]],
        user = t_start[["user.self"]], system = t_start[["sys.self"]],
        als_iterations = paste(iters, collapse = "+")
      )
      tgc_append(row, file.path(out, "svd_blas_threads.csv"))
      print(row, row.names = FALSE)
    }
  }
}

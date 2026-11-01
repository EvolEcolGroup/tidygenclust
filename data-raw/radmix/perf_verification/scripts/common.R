# Settings and helpers shared by the R scripts of the verification kit
# (data-raw/radmix/perf_verification/scripts), sourced by each of them. Run
# the scripts from the root of the repository. Environment variables:
# - TGC_WORK: working directory with the sources and libraries of the
#   builds (made by build_variants.sh); default ../tgc_perf_work
# - TGC_TMP: where the backing files of the gen_tibbles are written
#   (default: R's tempdir(); set it if /tmp is a small tmpfs, since the
#   full HGDP takes ~0.7 GB)
# - TGC_EUR_BED: the EUR_subset bed file (default: the cache of the tests)

tgc_repo <- normalizePath(".")
stopifnot(file.exists(file.path(tgc_repo, "DESCRIPTION")))
tgc_work <- Sys.getenv("TGC_WORK", file.path(dirname(tgc_repo), "tgc_perf_work"))

tgc_tmp <- Sys.getenv("TGC_TMP", "")
if (nzchar(tgc_tmp)) {
  # a directory of this process, removed when R exits
  tgc_tmp <- file.path(tgc_tmp, paste0("tgc_", Sys.getpid()))
  dir.create(tgc_tmp, recursive = TRUE, showWarnings = FALSE)
  reg.finalizer(globalenv(), function(e) unlink(tgc_tmp, recursive = TRUE),
    onexit = TRUE
  )
} else {
  tgc_tmp <- tempdir()
}

# the bed file of a dataset (see prepare_data.R):
# - eur: EUR_subset, 379 x 54,051, no missing genotypes
# - hgdp_quick: 209 x 25,000, with missing genotypes (called "hgdp" in the
#   laptop results)
# - hgdp: the HGDP of the benchmark article, 1043 x 130,943 (LD pruned)
# - hgdp_full: the HGDP without LD pruning, 1043 x 642,794
tgc_data <- function(name) {
  hgdp_dir <- file.path(tgc_repo, "vignettes", "articles", "data", "hgdp")
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
    stop(path, " not found: run scripts/prepare_data.R first")
  }
  path
}

# load the build of tidygenclust installed in `lib`, check that it is the
# one loaded, and return its namespace (for the internal kernels)
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

# a gen_tibble of a dataset, with its own backing file
tgc_read <- function(name) {
  tidypopgen::gen_tibble(tgc_data(name),
    backingfile = tempfile(name, tmpdir = tgc_tmp), quiet = TRUE
  )
}

# "eur:3,hgdp_quick:5" -> list of list(data = "eur", k = 3L), ...
tgc_cases <- function(s) {
  lapply(strsplit(strsplit(s, ",")[[1]], ":"), function(z) {
    list(data = z[1], k = as.integer(z[2]))
  })
}

tgc_ints <- function(s) as.integer(strsplit(s, ",")[[1]])

# the number of iterations of the first run of a gt_radmix() fit, from its
# log
tgc_iterations <- function(fit) {
  line <- grep("Converged in|Failed to converge", fit$log[[1]], value = TRUE)
  as.integer(sub(".* ([0-9]+) iterations.*", "\\1", line[length(line)]))
}

# append a data frame to a csv file (with a header if the file is new)
tgc_append <- function(df, file) {
  new <- !file.exists(file)
  utils::write.table(df, file,
    sep = ",", append = !new, col.names = new,
    row.names = FALSE
  )
}

# median elapsed time of n calls of f
tgc_median_time <- function(f, n) {
  stats::median(replicate(n, system.time(f())[["elapsed"]]))
}

# reorder the columns of q to match those of q_ref (cluster labels are
# arbitrary): exhaustive for K <= 7, greedy on the cost matrix above
tgc_match_clusters <- function(q, q_ref) {
  k <- ncol(q)
  cost <- outer(seq_len(k), seq_len(k), Vectorize(function(a, b) {
    sum(abs(q[, a] - q_ref[, b]))
  }))
  if (k <= 7) {
    perms <- function(n) {
      if (n == 1) {
        return(matrix(1))
      }
      sub <- perms(n - 1)
      do.call(rbind, lapply(seq_len(n), function(i) {
        cbind(i, matrix(setdiff(seq_len(n), i)[sub], nrow = nrow(sub)))
      }))
    }
    all <- perms(k)
    # all[r, b] = the column of q matched to column b of q_ref
    tot <- apply(all, 1, function(p) sum(cost[cbind(p, seq_len(k))]))
    return(q[, all[which.min(tot), ], drop = FALSE])
  }
  perm <- integer(k)
  free_a <- seq_len(k)
  free_b <- seq_len(k)
  while (length(free_a)) {
    sub <- cost[free_a, free_b, drop = FALSE]
    w <- which(sub == min(sub), arr.ind = TRUE)[1, ]
    perm[free_b[w[2]]] <- free_a[w[1]]
    free_a <- free_a[-w[1]]
    free_b <- free_b[-w[2]]
  }
  q[, perm, drop = FALSE]
}

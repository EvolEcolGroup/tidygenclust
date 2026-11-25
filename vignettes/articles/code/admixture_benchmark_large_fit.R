# One run of gt_radmix() on the simulated data of
# admixture_benchmark_large_data.R, from the command line, for the part of
# the benchmark on data larger than the memory (admixture_benchmark_large.sh
# runs it in memory and with the memory capped). It loads the gen_tibble,
# runs gt_radmix(), and appends a line to the results: its time, number of
# iterations and log-likelihood, the bytes it read from disk, its memory
# (the peak of its resident memory, and the cap and peak of its cgroup),
# whether the data were larger than the memory available to it (and so
# streamed), and how far its ancestry proportions are from those
# simulated. Usage, from the
# directory of the article (vignettes/articles):
#
#   Rscript code/admixture_benchmark_large_fit.R <run name> <K> <seed> \
#     <cores> <results file> <memory cap> <read throughput>
#
# (the memory cap, e.g. "8G" or "none", and the read throughput of the
# disk in MB/s, 0 for the disk as it is, are only recorded with the
# results: the caller sets them)

## ---- large-fit ----
args <- commandArgs(trailingOnly = TRUE)
run <- args[1]
k <- as.integer(args[2])
seed <- as.integer(args[3])
n_cores <- as.integer(args[4])
results_file <- args[5]
memory_cap <- args[6]
read_mbps <- as.numeric(args[7])
large_dir <- Sys.getenv("LARGE_DIR", "data/large")
library(tidygenclust)
sim_gt <- tidypopgen::gt_load(file.path(large_dir, "sim_large.gt"))
# the bytes read from disk by this process so far (Linux)
read_bytes <- function() {
  io <- readLines("/proc/self/io")
  as.numeric(sub(".*: ", "", grep("^read_bytes", io, value = TRUE)))
}
# whether gt_radmix() will stream the data: it does when the backing file is
# larger than 0.8 times the memory available to it (the smaller of the free
# memory and the room under its memory cap), and then reads it in chunks
stream <- tidygenclust:::radmix_stream_mode(gt_radmix_data(sim_gt))
mem_available <- tidygenclust:::radmix_mem_available()
# the peak of the resident memory of this process (Linux): its own memory
# plus the pages of the backing file mapped in it
peak_rss <- function() {
  x <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", x)) * 1024
}
# a memory file of the cgroup of this process (cgroups v2): memory.max is
# its cap, and memory.peak the peak of its memory, page cache included
# (NA without a cap: the cgroup is then the user's session)
cgroup_bytes <- function(name) {
  cg <- sub("^0::", "", grep("^0::", readLines("/proc/self/cgroup"),
    value = TRUE
  ))
  f <- file.path("/sys/fs/cgroup", cg, name)
  cap <- file.path("/sys/fs/cgroup", cg, "memory.max")
  if (!file.exists(f) || readLines(cap) == "max") {
    return(NA_real_)
  }
  as.numeric(readLines(f))
}
read_before <- read_bytes()
start <- proc.time()[["elapsed"]]
fit <- gt_radmix(sim_gt, k = k, seed = seed, n_cores = n_cores)
seconds <- proc.time()[["elapsed"]] - start
read_after <- read_bytes()

## ---- large-check ----
# the ancestry proportions estimated against those simulated: the clusters
# are numbered arbitrarily, so the columns are matched (by the permutation
# with the smallest sum of squared differences), then the root mean square
# difference is computed
q_true <- readRDS(file.path(large_dir, "q_true.rds"))
q <- fit$Q[[1]]
perms <- as.matrix(expand.grid(rep(list(seq_len(k)), k)))
perms <- perms[apply(perms, 1, function(x) length(unique(x)) == k), ]
sse <- apply(perms, 1, function(x) sum((q[, x] - q_true)^2))
rmse <- sqrt(min(sse) / length(q))
converged_line <- grep("^Converged in", fit$log[[1]], value = TRUE)[1]
result <- data.frame(
  run = run,
  memory_cap = memory_cap,
  read_mbps = read_mbps,
  k = k,
  n_cores = n_cores,
  seconds = seconds,
  iterations = as.integer(
    sub("^Converged in ([0-9]+) iterations.*", "\\1", converged_line)
  ),
  loglik = fit$loglik,
  read_gb = (read_after - read_before) / 1e9,
  peak_rss_gb = peak_rss() / 1e9,
  cgroup_max_gb = cgroup_bytes("memory.max") / 1e9,
  cgroup_peak_gb = cgroup_bytes("memory.peak") / 1e9,
  data_gb = stream$data_bytes / 1e9,
  mem_available_gb = mem_available / 1e9,
  streamed = stream$stream,
  chunk_mb = stream$chunk_bytes / 2^20,
  rmse_q = rmse
)
write.table(
  result, results_file,
  sep = ",", row.names = FALSE,
  col.names = !file.exists(results_file),
  append = file.exists(results_file)
)
print(result)

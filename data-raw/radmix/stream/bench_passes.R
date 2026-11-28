# Step 1 of dev_notes/radmix_streaming.md: what each pass of gt_radmix()
# costs when the backing file does not fit in the page cache. Runs the
# passes of a fit in the order of the fitting loop (EM, then the fused Q
# update, then P, Q, P, fused Q for each iteration), each timed on its own,
# with the bytes read from disk, the read throughput, the CPU used (all
# threads) and the major page faults. Run under bench_passes.sh (memory cap,
# read throttle, cold cache).
# Usage: Rscript bench_passes.R <prefix> <mode: base/FALSE/TRUE> <n_cores>
#          <k> <n_iter> <chunk_mb> <label>
args <- commandArgs(TRUE)
prefix <- args[1]
mode <- args[2]
n_cores <- as.integer(args[3])
k <- as.integer(args[4])
n_iter <- as.integer(args[5])
chunk_mb <- as.numeric(args[6])
label <- args[7]
fbm <- bigstatsr::big_attach(paste0(prefix, ".rds"))
lut <- rep(NA_real_, 256)
lut[1:3] <- 0:2
data <- list(
  fbm = fbm, rows = seq_len(nrow(fbm)), cols = seq_len(ncol(fbm)),
  kind = "code256", lut = lut, n_trials = rep(2, nrow(fbm)),
  scale = rep(1, nrow(fbm)), n_indiv = nrow(fbm), n_loci = ncol(fbm)
)
if (mode != "base") {
  data$stream <- as.logical(mode)
  data$chunk_bytes <- chunk_mb * 2^20
}
ns <- asNamespace("tidygenclust")
set.seed(1)
q <- matrix(stats::runif(k * nrow(fbm)), k)
q <- sweep(q, 2, colSums(q), "/")
p <- matrix(stats::runif(k * ncol(fbm), 0.05, 0.95), k)
io <- function() {
  x <- readLines("/proc/self/io")
  as.numeric(sub(".*: ", "", grep("^read_bytes", x, value = TRUE)))
}
majflt <- function() {
  as.numeric(strsplit(readLines("/proc/self/stat"), " ")[[1]][12])
}
file_mb <- file.size(paste0(prefix, ".bk")) / 2^20
# the passes of a fit; each returns the new q or p (or nothing)
passes <- list(
  em = function() {
    r <- ns$radmix_em_cpp(data, q, p, n_cores)
    q <<- r$q
    p <<- r$p
  },
  update_q_loglik = function() ns$radmix_update_q_loglik_cpp(data, q, p, n_cores),
  update_p = function() p <<- ns$radmix_update_p_cpp(data, q, p, n_cores),
  update_q = function() q <<- ns$radmix_update_q_cpp(data, q, p, n_cores),
  loglik = function() ns$radmix_loglik_cpp(data, q, p, n_cores)
)
order <- c("em", "update_q_loglik", rep(c(
  "update_p", "update_q", "update_p", "update_q_loglik"
), n_iter))
for (i in seq_along(order)) {
  r0 <- io()
  f0 <- majflt()
  t <- system.time(invisible(passes[[order[i]]]()))
  r1 <- io()
  f1 <- majflt()
  read_mb <- (r1 - r0) / 2^20
  cpu <- t[["user.self"]] + t[["sys.self"]]
  cat(sprintf("%s,%s,%d,%s,%.2f,%.0f,%.2f,%.0f,%.0f,%.2f\n",
    label, mode, i, order[i], t[["elapsed"]], read_mb,
    read_mb / file_mb, read_mb / t[["elapsed"]],
    100 * cpu / t[["elapsed"]] / n_cores, (f1 - f0) / 1000))
}

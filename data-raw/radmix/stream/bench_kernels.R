# Times of the passes of gt_radmix() (warm page cache), with and without
# streaming, for bench_stream.sh's datasets.
# Usage: Rscript bench_kernels.R <prefix> <modes, e.g. FALSE,TRUE or base>
#          <n_cores> <k> <reps> <chunk_mb>
args <- commandArgs(TRUE)
prefix <- args[1]
modes <- strsplit(args[2], ",")[[1]]
n_cores <- as.integer(args[3])
k <- as.integer(args[4])
reps <- as.integer(args[5])
chunk_mb <- as.numeric(args[6])
fbm <- bigstatsr::big_attach(paste0(prefix, ".rds"))
lut <- rep(NA_real_, 256)
lut[1:3] <- 0:2
data0 <- list(
  fbm = fbm, rows = seq_len(nrow(fbm)), cols = seq_len(ncol(fbm)),
  kind = "code256", lut = lut, n_trials = rep(2, nrow(fbm)),
  scale = rep(1, nrow(fbm)), n_indiv = nrow(fbm), n_loci = ncol(fbm)
)
set.seed(1)
q <- matrix(stats::runif(k * nrow(fbm)), k)
q <- sweep(q, 2, colSums(q), "/")
p <- matrix(stats::runif(k * ncol(fbm), 0.05, 0.95), k)
ns <- asNamespace("tidygenclust")
# warm the page cache
invisible(ns$radmix_loglik_cpp(data0, q, p, n_cores))
res <- NULL
for (r in seq_len(reps)) {
  for (m in modes) {
    data <- data0
    # modes: base, FALSE, TRUE, or TRUE_nopf (streaming without prefetch)
    if (m != "base") {
      data$stream <- m != "FALSE"
      data$stream_prefetch <- m != "TRUE_nopf"
      data$chunk_bytes <- chunk_mb * 2^20
    }
    tq <- system.time(ns$radmix_update_q_cpp(data, q, p, n_cores))[[3]]
    tql <- system.time(ns$radmix_update_q_loglik_cpp(data, q, p, n_cores))[[3]]
    tp <- system.time(ns$radmix_update_p_cpp(data, q, p, n_cores))[[3]]
    tl <- system.time(ns$radmix_loglik_cpp(data, q, p, n_cores))[[3]]
    te <- system.time(ns$radmix_em_cpp(data, q, p, n_cores))[[3]]
    res <- rbind(res, data.frame(mode = m, rep = r, update_q = tq,
      update_q_loglik = tql, update_p = tp, loglik = tl, em = te))
  }
}
print(aggregate(. ~ mode, res[, -2], median))

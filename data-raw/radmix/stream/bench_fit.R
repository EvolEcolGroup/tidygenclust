# One fit of gt_radmix()'s C++ loop on an FBM made by make_data.R, timed,
# with the bytes read from disk and a hash of the estimates, for
# bench_stream.sh (the elapsed time, then at the end the user and system
# CPU time of all the threads).
# Usage: Rscript bench_fit.R <prefix> <stream: TRUE/FALSE/base> <n_cores>
#          <k> <max_iter> <n_em_iter> <chunk_mb> <label>
args <- commandArgs(TRUE)
prefix <- args[1]
stream <- args[2]
n_cores <- as.integer(args[3])
k <- as.integer(args[4])
max_iter <- as.integer(args[5])
n_em_iter <- as.integer(args[6])
chunk_mb <- as.numeric(args[7])
label <- args[8]
fbm <- bigstatsr::big_attach(paste0(prefix, ".rds"))
lut <- rep(NA_real_, 256)
lut[1:3] <- 0:2
data <- list(
  fbm = fbm, rows = seq_len(nrow(fbm)), cols = seq_len(ncol(fbm)),
  kind = "code256", lut = lut, n_trials = rep(2, nrow(fbm)),
  scale = rep(1, nrow(fbm)), n_indiv = nrow(fbm), n_loci = ncol(fbm)
)
# modes: "base" (the build without streaming, which ignores these
# elements), FALSE (the usual layout), TRUE, TRUE_nopf (streaming without
# prefetching), TRUE_nozz (streaming without alternating directions) and
# CHUNKQ (not streaming, the Q side in chunks: chunk_q)
if (stream != "base") {
  data$stream <- !stream %in% c("FALSE", "CHUNKQ")
  data$chunk_q <- stream == "CHUNKQ"
  data$stream_prefetch <- stream != "TRUE_nopf"
  data$stream_zigzag <- stream != "TRUE_nozz"
  data$chunk_bytes <- chunk_mb * 2^20
}
io <- function() {
  x <- readLines("/proc/self/io")
  as.numeric(sub(".*: ", "", grep("^read_bytes", x, value = TRUE)))
}
majflt <- function() as.numeric(strsplit(readLines("/proc/self/stat"), " ")[[1]][12])
set.seed(42)
r0 <- io()
f0 <- majflt()
tm <- system.time(
  # (rtol 0 and no step_tol: exactly max_iter iterations)
  fit <- tidygenclust:::radmix_fit_cpp(data, k, NULL, NULL, FALSE, FALSE,
    max_iter, 0, Inf, n_em_iter, 3L, n_cores, FALSE, NULL)
)
t <- tm[["elapsed"]]
r1 <- io()
f1 <- majflt()
h <- rlang::hash(list(fit$q, fit$p, fit$loglik_history))
cat(sprintf("%s,%s,%d,%.2f,%.0f,%.0f,%s,%.10f,%.1f,%.1f\n", label, stream, n_cores, t,
  (r1 - r0) / 2^20, f1 - f0, h, fit$loglik, tm[["user.self"]], tm[["sys.self"]]))

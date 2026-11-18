# Single passes of gt_radmix() over the genotype likelihoods of a gl_tibble
# (through tidypopgl's source), and optionally over the same likelihoods held
# in memory as doubles (tidygenclust's example source, which decodes
# nothing), plus a whole fit. One build of tidypopgl per call, so that builds
# can be alternated in separate processes.
#
# Usage:
#   Rscript time_gl_passes.R <tidygenclust lib> <tidypopgl lib> <label> \
#     <n_cores> <reps> <out csv> [dense] [fit]
# dense: also time the in-memory example source ("dense" to do it);
# fit: also time a whole fit ("fit").
# Environment: GL_DATA (the .gl file), GL_K (3).
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[2], args[1], .libPaths()))
label <- args[3]
n_cores <- as.integer(args[4])
reps <- as.integer(args[5])
out_csv <- args[6]
do_dense <- "dense" %in% args[-(1:6)]
do_fit <- "fit" %in% args[-(1:6)]
gl_data <- Sys.getenv(
  "GL_DATA", "/home/andrea/tgc_supfused_work/data/sim_400_100000.gl"
)
k <- as.integer(Sys.getenv("GL_K", "3"))
suppressPackageStartupMessages({
  library(tidygenclust)
  library(tidypopgl)
})
ns <- asNamespace("tidygenclust")
gl <- gl_load(gl_data, quiet = TRUE)
sources <- list(gl16 = gt_radmix_data(gl))
if (do_dense) {
  # the decoded likelihoods, as doubles in memory: e1 = E(u1), e2 = E(u2)
  # with E(u) = exp((u - 32768) / 1024) (the gl16 decoding), obs = not both
  # codes missing
  fbm <- tidypopgl:::.gl_kernel_args(gl)$fbm
  codes <- fbm[, ]
  u1 <- codes[seq(1, nrow(codes), by = 2), , drop = FALSE]
  u2 <- codes[seq(2, nrow(codes), by = 2), , drop = FALSE]
  rm(codes)
  obs <- !(u1 == 32768 & u2 == 32768)
  e1 <- exp((u1 - 32768) / 1024)
  e2 <- exp((u2 - 32768) / 1024)
  rm(u1, u2)
  sources$dense <- c(
    list(kind = "gl_source", version = 1L),
    ns$radmix_gl_dense_source_cpp(e1, e2, obs),
    list(n_indiv = nrow(e1), n_loci = ncol(e1), e1 = e1, e2 = e2, obs = obs)
  )
  invisible(gc())
}
n_indiv <- nrow(gl)
n_loci <- sources$gl16$n_loci

# the same estimates for every build: random q on the simplex, p in
# (0.05, 0.95)
set.seed(42)
q <- matrix(stats::runif(k * n_indiv), k, n_indiv)
q <- sweep(q, 2, colSums(q), "/")
p <- matrix(stats::runif(k * n_loci, 0.05, 0.95), k, n_loci)
passes <- list(
  loglik = function(d) ns$radmix_loglik_cpp(d, q, p, n_cores),
  em = function(d) ns$radmix_em_cpp(d, q, p, n_cores),
  update_q = function(d) ns$radmix_update_q_cpp(d, q, p, n_cores),
  update_q_loglik = function(d) ns$radmix_update_q_loglik_cpp(d, q, p, n_cores),
  update_p = function(d) ns$radmix_update_p_cpp(d, q, p, n_cores)
)
load_avg <- function() strsplit(readLines("/proc/loadavg"), " ")[[1]][1]
rows <- list()
check <- list()
for (src in names(sources)) {
  d <- sources[[src]]
  for (pass in names(passes)) {
    f <- passes[[pass]]
    res <- f(d)
    if (src == "gl16") check[[pass]] <- res
    times <- vapply(seq_len(reps), function(r) {
      system.time(f(d))[["elapsed"]]
    }, numeric(1))
    rows[[length(rows) + 1]] <- data.frame(
      label = label, source = src, what = pass, n_cores = n_cores,
      seconds = stats::median(times), min_seconds = min(times),
      max_seconds = max(times), iterations = NA_integer_, loglik = NA_real_,
      load = load_avg(), date = format(Sys.time())
    )
  }
  if (do_fit && src == "gl16") {
    t_fit <- system.time(
      fit <- gt_radmix(gl, k = k, seed = 1, n_cores = n_cores)
    )[["elapsed"]]
    n_iter <- as.integer(sub(
      ".* (in|after) ([0-9]+) iterations.*", "\\2",
      grep("( in [0-9]+ iterations| after [0-9]+ iterations)",
        fit$log[[1]], value = TRUE
      )
    ))
    check$fit_q <- unclass(fit$Q[[1]])
    rows[[length(rows) + 1]] <- data.frame(
      label = label, source = src, what = "fit", n_cores = n_cores,
      seconds = t_fit, min_seconds = t_fit, max_seconds = t_fit,
      iterations = n_iter, loglik = fit$loglik, load = load_avg(),
      date = format(Sys.time())
    )
  }
}
res <- do.call(rbind, rows)
print(res[, c("label", "source", "what", "n_cores", "seconds", "iterations")],
      digits = 4)
utils::write.table(res, out_csv,
  sep = ",", row.names = FALSE,
  col.names = !file.exists(out_csv), append = file.exists(out_csv)
)
# the results of the passes, to check that builds agree
saveRDS(check, file.path(dirname(out_csv),
                         paste0("check_", label, "_", n_cores, ".rds")))

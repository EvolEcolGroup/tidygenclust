# Random starts against SVD starts (init = "svd", the posterior dosages) for
# gt_radmix() on the simulated gl_tibble (400 x 100,000, K = 3 simulated):
# time, iterations and log-likelihood of each fit, K = 3 and 5, seeds 1-3.
# Usage: Rscript svd_eval.R <tidygenclust lib> <tidypopgl lib> <n_cores> <out csv>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[2], args[1], .libPaths()))
n_cores <- as.integer(args[3])
suppressPackageStartupMessages({
  library(tidygenclust)
  library(tidypopgl)
})
gl <- gl_load(Sys.getenv("GL_DATA",
  "/home/andrea/tgc_supfused_work/data/sim_400_100000.gl"), quiet = TRUE)
rows <- list()
for (k in c(3, 5)) {
  for (seed in 1:3) {
    for (init in c("random", "svd")) {
      t <- system.time(fit <- gt_radmix(gl, k = k, seed = seed, init = init,
                                        n_cores = n_cores))[["elapsed"]]
      log <- fit$log[[1]]
      n_iter <- as.integer(sub(".* (in|after) ([0-9]+) iterations.*", "\\2",
        grep("( in [0-9]+ iterations| after [0-9]+ iterations)", log,
             value = TRUE)))
      rows[[length(rows) + 1]] <- data.frame(
        k = k, seed = seed, init = init, n_cores = n_cores, time_s = t,
        iterations = n_iter, loglik = fit$loglik,
        load = strsplit(readLines("/proc/loadavg"), " ")[[1]][1]
      )
      message(k, " ", seed, " ", init, ": ", round(t, 2), " s, ", n_iter,
              " iterations, loglik ", format(fit$loglik, digits = 12))
    }
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, args[4], row.names = FALSE)

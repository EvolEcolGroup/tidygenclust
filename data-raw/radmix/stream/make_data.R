# Simulated genotypes for the benchmark of streaming (bench_stream.sh):
# an FBM.code256 of n_indiv individuals and n_loci loci, admixed from 3
# populations, with 1% missing genotypes, saved as <prefix>.rds/.bk.
# Usage: Rscript make_data.R <n_indiv> <n_loci> <prefix>
args <- commandArgs(TRUE)
n_indiv <- as.integer(args[1])
n_loci <- as.integer(args[2])
prefix <- args[3]
set.seed(1)
k <- 3
q <- matrix(stats::rgamma(n_indiv * k, 0.3), n_indiv, k)
q <- q / rowSums(q)
fbm <- bigstatsr::FBM.code256(n_indiv, n_loci,
  code = c(0, 1, 2, NA, rep(NA, 252)),
  backingfile = prefix, init = NULL
)
block <- max(1L, floor(5e7 / n_indiv))
for (j0 in seq(1, n_loci, by = block)) {
  cols <- j0:min(n_loci, j0 + block - 1)
  p <- matrix(stats::rbeta(k * length(cols), 0.5, 0.5), k, length(cols))
  p <- pmin(pmax(p, 0.01), 0.99)
  qp <- q %*% p
  g <- stats::rbinom(length(qp), 2, qp)
  g[stats::runif(length(g)) < 0.01] <- 3L
  fbm[, cols] <- matrix(g, n_indiv)
}
fbm$save()
cat("written", prefix, n_indiv, "x", n_loci, "\n")

# Smoke test of the article's scripts for data larger than the memory, at a
# small size (20,000 x 75,000: 1.5 GB of genotypes, a 1 GiB cap), all three
# runs (the throttled one with sudo): three lines of results, the capped
# runs streamed and stayed under their cap, and the same log-likelihood.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
set -e
export R_LIBS=$TGC_WORK/libs/art:$W/lib_tpg LARGE_DIR=$W/large_smoke
export LARGE_CORES=16 LARGE_N_INDIV=20000 LARGE_N_LOCI=75000 LARGE_MEMORY_CAP=1G
rm -rf "$LARGE_DIR"
cd "$TGC_WORK/src/art/vignettes/articles"
Rscript code/admixture_benchmark_large_data.R
taskset -c 0-7 bash code/admixture_benchmark_large.sh
mkdir -p "$RES/large_smoke"
cp data/admixture_benchmark/large_*.csv "$RES/large_smoke/"
Rscript -e 'r <- read.csv("data/admixture_benchmark/large_runs.csv")
print(r); print(read.csv("data/admixture_benchmark/large_memory.csv"))
stopifnot(nrow(r) == 3, !r$streamed[1], all(r$streamed[-1]),
  all(r$cgroup_peak_gb[-1] <= r$cgroup_max_gb[-1]),
  length(unique(r$loglik)) == 1)'
rm -rf "$LARGE_DIR"

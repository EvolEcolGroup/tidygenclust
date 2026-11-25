# The article's fits: in memory, capped at 8 GiB on the NVMe, and capped
# with the reads throttled to 500 MB/s (sudo), on the P-cores; the results
# are written into the source copy, then copied to results/large.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
set -e
export R_LIBS=$TGC_WORK/libs/art:$W/lib_tpg LARGE_DIR=$W/large
cd "$TGC_WORK/src/art/vignettes/articles"
taskset -c 0-7 bash code/admixture_benchmark_large.sh
mkdir -p "$RES/large"
cp data/admixture_benchmark/large_*.csv "$RES/large/"
cat "$RES/large/large_runs.csv" "$RES/large/large_memory.csv"

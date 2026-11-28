# The 5 GB data in memory, then capped at 4 GiB (a little less than the
# data) and 1 GiB (a fifth), on the NVMe and throttled to 500 and 150 MB/s
# (sudo); the slowest run (1 GiB at 150 MB/s) last.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
set -e
export R_LIBS=$TGC_WORK/libs/art:$W/lib_tpg LARGE_DIR=$W/mid LARGE_NAME=mid
export LARGE_MEMORY_CAPS="4G 1G" LARGE_READ_MBPS="0 500 150"
cd "$TGC_WORK/src/art/vignettes/articles"
taskset -c 0-7 bash code/admixture_benchmark_large.sh
mkdir -p "$RES/mid"
cp data/admixture_benchmark/mid_*.csv "$RES/mid/"
cat "$RES/mid/mid_runs.csv" "$RES/mid/mid_memory.csv"

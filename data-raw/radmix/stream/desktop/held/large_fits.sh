# The article's fits: in memory and capped at 8 GiB, on the P-cores; the
# results are written into the source copy, then copied to results/large.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
export R_LIBS=$TGC_WORK/libs/new2:$W/lib_tpg LARGE_DIR=$W/large
cd "$TGC_WORK/src/new2/vignettes/articles"
taskset -c 0-7 bash code/admixture_benchmark_large.sh
mkdir -p "$RES/large"
cp data/admixture_benchmark/large_*.csv "$RES/large/"

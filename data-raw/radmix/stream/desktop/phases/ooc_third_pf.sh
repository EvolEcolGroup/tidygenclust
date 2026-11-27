# Rerun of ooc_third with the prefetch fixed (new2, b07f087: madvise in
# pieces of 128 KB): streaming only (all on, no prefetch, no alternation);
# base is unchanged, its rows are in desktop_third.csv.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
[ -e "$W/done/throttle_check" ] || { echo "throttle_check not done"; exit 1; }
cd "$BENCH_DATA"
LIB_NEW=$TGC_WORK/libs/new2:$W/lib_tpg NO_BASE=1 \
  taskset -c 0-7 sh "$STREAM/run_ooc_third.sh" "$RES/desktop_third_pf128.csv" 32

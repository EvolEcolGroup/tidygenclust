# Phase 2-3: full fits with 67% and 92% of the file cacheable, the NVMe as it
# is and throttled to 500 and 150 MB/s.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
[ -e "$W/done/throttle_check" ] || { echo "throttle_check not done"; exit 1; }
cd "$BENCH_DATA"
taskset -c 0-7 sh "$STREAM/run_ooc.sh" "$RES/desktop_ooc.csv" 2 1 64 ""

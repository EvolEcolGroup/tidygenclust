# Phase 4: the cost of each pass, base and streaming, in memory and with 33%,
# 67% and 92% of the file cacheable, on the NVMe as it is.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
cd "$BENCH_DATA"
taskset -c 0-7 sh "$STREAM/run_passes.sh" "$RES/desktop_passes.csv" "mem:0 33:0 67:0 92:0"

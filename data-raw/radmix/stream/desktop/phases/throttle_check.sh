# The caps work: a read of the 1.5 GB file throttled to 150 MB/s must take
# at least 9 s (10 s at the set rate), and a read under a 600 MB cap must
# be in a run-*.scope with memory.max set. The phases with throttled rows
# refuse to start unless this phase is done.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
cd "$BENCH_DATA"
set -e
t0=$(date +%s.%N)
sh "$STREAM/run_capped.sh" 0 150 big_a.bk "" sh -c \
  'cg=$(cut -d: -f3 /proc/self/cgroup); echo $cg; echo "io.max: $(cat /sys/fs/cgroup$cg/io.max)"; dd if=big_a.bk of=/dev/null bs=4M 2>&1 | tail -1'
t1=$(date +%s.%N)
echo "throttled read: $(echo "$t1 - $t0" | bc) s"
awk -v a="$t0" -v b="$t1" 'BEGIN { exit !(b - a >= 9) }'
sh "$STREAM/run_capped.sh" 600 0 big_a.bk "" sh -c \
  'cg=$(cut -d: -f3 /proc/self/cgroup); echo $cg; cat /sys/fs/cgroup$cg/memory.max; for i in 1 2; do dd if=big_a.bk of=/dev/null bs=4M 2>&1 | tail -1; done' |
  tee -a /dev/stderr | grep -q '^629145600$'
echo "caps OK"

#!/bin/sh
# Step 1 of dev_notes/radmix_streaming.md: the passes of dev (and, for
# comparison, of streaming) in memory and with 60% or 85% of the file
# cacheable, on the disk as it is and throttled to 500 MB/s.
# The memory caps are the peak memory of bench_passes.R itself
# (anonp_<prefix>.txt, measured by sampling the cgroup's rss) plus the
# fraction of the file (1430 MB).
# Usage: run_passes.sh <out.csv> [configs, default "mem:0 60:0 60:500 85:0"]
OUT=$1
CFGS=${2:-"mem:0 60:0 60:500 85:0"}
here=$(cd "$(dirname "$0")" && pwd)
cd /home/user/benchdata
[ -f "$OUT" ] || echo "label,mode,i,pass,elapsed_s,read_mb,file_reads,read_mbps,cpu_pct,major_faults_k,mem_cap_mb,throttle_mbps" > "$OUT"
for p in big_a big_b; do
  anon=$(cat anonp_$p.txt)
  for cfg in $CFGS; do
    frac=${cfg%%:*}; bps=${cfg##*:}
    cap=0; [ $frac = mem ] || cap=$((anon + 1430 * frac / 100))
    for m in base TRUE; do
      lib=/home/user/rlib_new; [ $m = base ] && lib=/home/user/rlib_base
      sh $here/bench_passes.sh $lib "$p:$frac:$bps" $p $m 4 4 2 64 $cap $bps >> "$OUT"
      echo "done $p $cfg $m"
    done
  done
done
echo PASSES_DONE

#!/bin/sh
# Step 1 of dev_notes/radmix_streaming.md: the passes of dev (and, for
# comparison, of streaming) in memory and with 60% or 85% of the file
# cacheable, on the disk as it is and throttled to 500 MB/s.
# The memory caps are the peak memory of bench_passes.R itself
# (anonp_<prefix>.txt, from measure_peak.sh) plus the fraction of the file.
# Environment: BENCH_DATA (the folder of the data of make_data.R, default
# /home/user/benchdata), LIB_BASE and LIB_NEW (the libraries of the builds of
# dev and of the branch, default /home/user/rlib_base and rlib_new),
# N_CORES (default 4), FILE_MB (the size of the backing files, default 1430).
# Usage: run_passes.sh <out.csv> [configs, default "mem:0 60:0 60:500 85:0"]
OUT=$1
CFGS=${2:-"mem:0 60:0 60:500 85:0"}
here=$(cd "$(dirname "$0")" && pwd)
cd "${BENCH_DATA:-/home/user/benchdata}"
LIB_BASE=${LIB_BASE:-/home/user/rlib_base}
LIB_NEW=${LIB_NEW:-/home/user/rlib_new}
N=${N_CORES:-4}
FILE_MB=${FILE_MB:-1430}
[ -f "$OUT" ] || echo "label,mode,i,pass,elapsed_s,read_mb,file_reads,read_mbps,cpu_pct,major_faults_k,mem_cap_mb,throttle_mbps" > "$OUT"
for p in big_a big_b; do
  anon=$(cat anonp_$p.txt)
  for cfg in $CFGS; do
    frac=${cfg%%:*}; bps=${cfg##*:}
    cap=0; [ $frac = mem ] || cap=$((anon + FILE_MB * frac / 100))
    for m in base TRUE; do
      lib=$LIB_NEW; [ $m = base ] && lib=$LIB_BASE
      sh $here/bench_passes.sh $lib "$p:$frac:$bps" $p $m $N 4 2 64 $cap $bps >> "$OUT"
      echo "done $p $cfg $m"
    done
  done
done
echo PASSES_DONE

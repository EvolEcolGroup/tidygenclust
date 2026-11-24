#!/bin/sh
# Data about 3 times the memory: the full fit of bench_fit.R with a third
# of the file cacheable (memory cap: the peak memory of the fit, in
# anon_<prefix>.txt from measure_peak.sh, plus a third of the file), on
# the disk as it is and throttled to 500 and
# 150 MB/s, with streaming (all on, without prefetching, without
# alternating directions) and dev ("base"; not on the slowest disk for
# 10,000 individuals, where it would take about 35 minutes).
# Environment: BENCH_DATA (the folder of the data of make_data.R, default
# /home/user/benchdata), LIB_BASE and LIB_NEW (the libraries of the builds of
# dev and of the branch, default /home/user/rlib_base and rlib_new),
# N_CORES (default 4), FILE_MB (the size of the backing files, default 1430).
# Usage: run_ooc_third.sh <out.csv> <chunk_mb>
OUT=$1; chunk=$2
here=$(cd "$(dirname "$0")" && pwd)
cd "${BENCH_DATA:-/home/user/benchdata}"
LIB_BASE=${LIB_BASE:-/home/user/rlib_base}
LIB_NEW=${LIB_NEW:-/home/user/rlib_new}
N=${N_CORES:-4}
FILE_MB=${FILE_MB:-1430}
for p in big_a big_b; do
  peak=$(cat anon_$p.txt)
  cap=$((peak + FILE_MB / 3))
  for bps in 0 500 150; do
    modes="TRUE TRUE_nopf TRUE_nozz base"
    [ $p = big_a ] && [ $bps = 150 ] && modes="TRUE TRUE_nopf TRUE_nozz"
    for m in $modes; do
      lib=$LIB_NEW; [ $m = base ] && lib=$LIB_BASE
      line=$(sh $here/bench_stream.sh $lib "$p:33" $p $m $N 4 2 1 $chunk $cap $bps)
      echo "$line" >> "$OUT"
      echo "$line"
    done
  done
done
echo THIRD_DONE

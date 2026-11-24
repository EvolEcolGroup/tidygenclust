#!/bin/sh
# Data about 3 times the memory: the full fit of bench_fit.R with a third
# of the file cacheable (memory cap: the peak memory of the fit, 166 MB and
# 549 MB, plus 477 MB), on the disk as it is and throttled to 500 and
# 150 MB/s, with streaming (all on, without prefetching, without
# alternating directions) and dev ("base"; not on the slowest disk for
# 10,000 individuals, where it would take about 35 minutes).
# Usage: run_ooc_third.sh <out.csv> <chunk_mb>
OUT=$1; chunk=$2
here=$(cd "$(dirname "$0")" && pwd)
cd /home/user/benchdata
for p in big_a big_b; do
  peak=166; [ $p = big_b ] && peak=549
  cap=$((peak + 477))
  for bps in 0 500 150; do
    modes="TRUE TRUE_nopf TRUE_nozz base"
    [ $p = big_a ] && [ $bps = 150 ] && modes="TRUE TRUE_nopf TRUE_nozz"
    for m in $modes; do
      lib=/home/user/rlib_new; [ $m = base ] && lib=/home/user/rlib_base
      line=$(sh $here/bench_stream.sh $lib "$p:33" $p $m 4 4 2 1 $chunk $cap $bps)
      echo "$line" >> "$OUT"
      echo "$line"
    done
  done
done
echo THIRD_DONE

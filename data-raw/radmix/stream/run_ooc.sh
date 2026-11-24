#!/bin/sh
# The out-of-core benchmark of dev_notes/radmix_streaming.md: each dataset
# under memory caps that leave the page cache a fraction of the data (60%
# and 85%), on the disk as it is and throttled (500 and 150 MB/s), with the
# build of dev ("base") and streaming (all on, without prefetching, and
# without alternating directions). Appends CSV lines to <out.csv>.
# Usage: run_ooc.sh <out.csv> <iter> <em> <chunk_mb> <skip_base>
# (skip_base: dataset:frac:bps triples whose base run is already measured)
OUT=$1; iter=$2; em=$3; chunk=$4; skip=$5
here=$(cd "$(dirname "$0")" && pwd)
cd /home/user/benchdata
B=$here/bench_stream.sh
for p in big_a big_b; do
  # the memory of the R process (measured, plus 100 MB) and a fraction of
  # the data (1430 MB)
  anon=$(cat anon_$p.txt)
  for frac in 60 85; do
    cap=$((anon + 1430 * frac / 100))
    for bps in 0 500 150; do
      modes="TRUE TRUE_nopf TRUE_nozz base"
      # (on the slowest disk, streaming with and without prefetching only)
      [ $bps = 150 ] && modes="TRUE TRUE_nopf"
      for m in $modes; do
        lib=/home/user/rlib_new; [ $m = base ] && lib=/home/user/rlib_base
        if [ $m = base ] && echo " $skip " | grep -q " $p:$frac:$bps "; then
          continue
        fi
        line=$(sh $B $lib "$p:$frac" $p $m 4 4 $iter $em $chunk $cap $bps)
        echo "$line" >> "$OUT"
        echo "$line"
      done
    done
  done
done
echo OOC_DONE

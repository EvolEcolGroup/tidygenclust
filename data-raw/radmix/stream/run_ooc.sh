#!/bin/sh
# The out-of-core benchmark of dev_notes/radmix_streaming.md: each dataset
# under memory caps that leave the page cache a fraction of the data (60%
# and 85%), on the disk as it is and throttled (500 and 150 MB/s), with the
# build of dev ("base") and streaming (all on, without prefetching, and
# without alternating directions). Appends CSV lines to <out.csv>.
# Environment: BENCH_DATA (the folder of the data of make_data.R, default
# /home/user/benchdata), LIB_BASE and LIB_NEW (the libraries of the builds of
# dev and of the branch, default /home/user/rlib_base and rlib_new),
# N_CORES (default 4), FILE_MB (the size of the backing files, default 1430).
# Usage: run_ooc.sh <out.csv> <iter> <em> <chunk_mb> <skip_base>
# (skip_base: dataset:frac:bps triples whose base run is already measured)
OUT=$1; iter=$2; em=$3; chunk=$4; skip=$5
here=$(cd "$(dirname "$0")" && pwd)
cd "${BENCH_DATA:-/home/user/benchdata}"
LIB_BASE=${LIB_BASE:-/home/user/rlib_base}
LIB_NEW=${LIB_NEW:-/home/user/rlib_new}
N=${N_CORES:-4}
FILE_MB=${FILE_MB:-1430}
B=$here/bench_stream.sh
for p in big_a big_b; do
  # the peak memory of the fit (anon_<prefix>.txt, from measure_peak.sh)
  # and a fraction of the file
  anon=$(cat anon_$p.txt)
  for frac in 60 85; do
    cap=$((anon + FILE_MB * frac / 100))
    for bps in 0 500 150; do
      modes="TRUE TRUE_nopf TRUE_nozz base"
      # (on the slowest disk, streaming with and without prefetching only)
      [ $bps = 150 ] && modes="TRUE TRUE_nopf"
      for m in $modes; do
        lib=$LIB_NEW; [ $m = base ] && lib=$LIB_BASE
        if [ $m = base ] && echo " $skip " | grep -q " $p:$frac:$bps "; then
          continue
        fi
        line=$(sh $B $lib "$p:$frac" $p $m $N 4 $iter $em $chunk $cap $bps)
        echo "$line" >> "$OUT"
        echo "$line"
      done
    done
  done
done
echo OOC_DONE

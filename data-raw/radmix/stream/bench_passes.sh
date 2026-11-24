#!/bin/sh
# bench_passes.R from a cold page cache, under a memory cap and a read
# throttle (run_capped.sh; 0 for none). Prints its CSV lines, with the caps.
# Usage: bench_passes.sh <lib> <label> <prefix> <mode> <n_cores> <k>
#          <n_iter> <chunk_mb> <mem_cap_mb> <read_mbps>
lib=$1; label=$2; prefix=$3; mode=$4; cores=$5; k=$6; iter=$7; chunk=$8
cap=$9; shift 9; bps=$1
here=$(cd "$(dirname "$0")" && pwd)
R_LIBS=$lib sh "$here/run_capped.sh" "$cap" "$bps" "$prefix.bk" "" \
  Rscript "$here/bench_passes.R" "$prefix" "$mode" "$cores" "$k" "$iter" "$chunk" "$label" |
  sed "s/\$/,$cap,$bps/"

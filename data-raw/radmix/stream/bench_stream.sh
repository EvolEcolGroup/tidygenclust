#!/bin/sh
# Benchmark of streaming in gt_radmix() (dev_notes/radmix_streaming.md):
# one fit of bench_fit.R, from a cold page cache (unless keep_cache), under
# a memory cap that also bounds the page cache and a cap on the read
# throughput of the disk (run_capped.sh), to emulate data larger than the
# memory on a slower disk. Prints the CSV line of bench_fit.R, with the
# caps.
# Usage: bench_stream.sh <lib> <label> <prefix> <stream> <n_cores> <k>
#          <max_iter> <n_em_iter> <chunk_mb> <mem_cap_mb or 0>
#          <read_mbps or 0> [keep_cache]
lib=$1; label=$2; prefix=$3; stream=$4; cores=$5; k=$6; iter=$7; em=$8
chunk=$9; shift 9; cap=$1; bps=$2; keep=$3
here=$(cd "$(dirname "$0")" && pwd)
out=$(R_LIBS=$lib sh "$here/run_capped.sh" "$cap" "$bps" "$prefix.bk" "$keep" \
  Rscript "$here/bench_fit.R" "$prefix" "$stream" "$cores" "$k" "$iter" "$em" "$chunk" "$label")
echo "$out,$cap,$bps"

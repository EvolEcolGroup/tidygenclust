#!/bin/sh
# Benchmark of streaming in gt_radmix() (dev_notes/radmix_streaming.md):
# one fit of bench_fit.R, from a cold page cache, optionally under a memory
# cap (cgroup v1 memory.limit_in_bytes, which also bounds the page cache)
# and a cap on the read throughput of the disk (blkio throttle), to emulate
# data larger than the memory on a slower disk. Needs root (drop_caches and
# the cgroups). Prints the CSV line of bench_fit.R, with the caps.
# Usage: bench_stream.sh <lib> <label> <prefix> <stream> <n_cores> <k>
#          <max_iter> <n_em_iter> <chunk_mb> <mem_cap_mb or 0>
#          <read_mbps or 0> [keep_cache]
lib=$1; label=$2; prefix=$3; stream=$4; cores=$5; k=$6; iter=$7; em=$8
chunk=$9; shift 9; cap=$1; bps=$2; keep=$3
here=$(cd "$(dirname "$0")" && pwd)
memcg=$(awk -F: '$2 ~ /(^|,)memory(,|$)/ {print $3}' /proc/self/cgroup)
memcg=/sys/fs/cgroup/memory$memcg/radbench
blkcg=/sys/fs/cgroup/blkio/radbench
dev=$(lsblk -no MAJ:MIN "$(df --output=source "$prefix.bk" | tail -1)" 2>/dev/null | head -1 | tr -d ' ')
[ -n "$dev" ] || dev=254:0
mkdir -p "$memcg" "$blkcg"
if [ "$cap" -gt 0 ]; then
  echo $((cap * 1024 * 1024)) > "$memcg/memory.limit_in_bytes"
else
  echo -1 > "$memcg/memory.limit_in_bytes"
fi
if [ "$bps" -gt 0 ]; then
  echo "$dev $((bps * 1000 * 1000))" > "$blkcg/blkio.throttle.read_bps_device"
else
  echo "$dev 0" > "$blkcg/blkio.throttle.read_bps_device" 2>/dev/null
fi
[ -n "$keep" ] || { sync; echo 1 > /proc/sys/vm/drop_caches; }
out=$(sh -c "echo \$\$ > $memcg/cgroup.procs; echo \$\$ > $blkcg/cgroup.procs; \
  R_LIBS=$lib exec Rscript $here/bench_fit.R $prefix $stream $cores $k $iter $em $chunk $label")
echo "$out,$cap,$bps"

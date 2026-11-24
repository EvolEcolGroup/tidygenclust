#!/bin/sh
# bench_passes.R from a cold page cache, under a memory cap (cgroup v1,
# which bounds the page cache too; 0 for none) and a read throttle (blkio;
# 0 for none). Needs root.
# Usage: bench_passes.sh <lib> <label> <prefix> <mode> <n_cores> <k>
#          <n_iter> <chunk_mb> <mem_cap_mb> <read_mbps>
lib=$1; label=$2; prefix=$3; mode=$4; cores=$5; k=$6; iter=$7; chunk=$8
cap=$9; shift 9; bps=$1
here=$(cd "$(dirname "$0")" && pwd)
memcg=/sys/fs/cgroup/memory$(awk -F: '$2 ~ /(^|,)memory(,|$)/ {print $3}' /proc/self/cgroup)/radbench
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
sync; echo 1 > /proc/sys/vm/drop_caches
sh -c "echo \$\$ > $memcg/cgroup.procs; echo \$\$ > $blkcg/cgroup.procs; \
  R_LIBS=$lib exec Rscript $here/bench_passes.R $prefix $mode $cores $k $iter $chunk $label" |
  sed "s/\$/,$cap,$bps/"

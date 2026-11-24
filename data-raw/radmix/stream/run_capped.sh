#!/bin/sh
# Run a command from a cold page cache (unless keep_cache), under a memory
# cap that also bounds the page cache, and a cap on the read throughput of
# the disk that holds <data_file>, for the benchmarks of
# dev_notes/radmix_streaming.md. Two ways, chosen automatically:
# - cgroup v1, as root (the cloud container): a child cgroup "radbench" of
#   this shell's memory and blkio cgroups (memory.limit_in_bytes,
#   blkio.throttle.read_bps_device);
# - systemd with cgroup v2 (the desktop): a transient scope. Without a
#   throttle, a scope of the user's systemd (`systemd-run --user --scope -p
#   MemoryMax=...`), which needs no root, as the memory controller is
#   delegated to users; with a throttle, a scope of the system's
#   (`sudo -n systemd-run --scope -p IOReadBandwidthMax=...`), in which the
#   command runs as the calling user (sudo -u), since the io controller is
#   not delegated. sudo must not ask for a password in the middle of a long
#   matrix: run `sudo -v` first in the same terminal and keep it alive, or
#   the throttled runs fail at once (sudo -n).
# The cold page cache: as root, the whole page cache is dropped; otherwise
# the pages of <data_file> alone (dd iflag=nocache, i.e.
# posix_fadvise(POSIX_FADV_DONTNEED), which needs no root and is all that
# matters here: the other files read are small).
# Usage: run_capped.sh <mem_cap_mb or 0> <read_mbps or 0> <data_file>
#          <keep_cache: "" or keep> <command> [args...]
cap=$1; bps=$2; data_file=$3; keep=$4; shift 4
if [ -z "$keep" ]; then
  sync
  if [ "$(id -u)" = 0 ]; then
    echo 1 > /proc/sys/vm/drop_caches
  else
    dd if="$data_file" iflag=nocache count=0 status=none
  fi
fi
v1mem=$(awk -F: '$2 ~ /(^|,)memory(,|$)/ {print $3}' /proc/self/cgroup)
if [ -n "$v1mem" ] && [ -d /sys/fs/cgroup/memory ] && [ "$(id -u)" = 0 ]; then
  memcg=/sys/fs/cgroup/memory$v1mem/radbench
  blkcg=/sys/fs/cgroup/blkio/radbench
  dev=$(lsblk -no MAJ:MIN "$(df --output=source "$data_file" | tail -1)" 2>/dev/null | head -1 | tr -d ' ')
  [ -n "$dev" ] || dev=254:0
  mkdir -p "$memcg" "$blkcg"
  if [ "$cap" -gt 0 ]; then
    echo $((cap * 1024 * 1024)) > "$memcg/memory.limit_in_bytes"
  else
    echo -1 > "$memcg/memory.limit_in_bytes"
  fi
  echo "$dev $((bps * 1000 * 1000))" > "$blkcg/blkio.throttle.read_bps_device" 2>/dev/null
  # (the shell moves itself into the cgroups, then becomes the command)
  exec sh -c 'echo $$ > "$1/cgroup.procs"; echo $$ > "$2/cgroup.procs"; shift 2; exec "$@"' \
    sh "$memcg" "$blkcg" "$@"
elif command -v systemd-run > /dev/null && [ -d /run/systemd/system ]; then
  if [ "$bps" -gt 0 ]; then
    # a system scope, with the command run as the calling user, and the
    # environment it needs (sudo resets it)
    set -- sudo -n -u "${SUDO_USER:-$(id -un)}" env "R_LIBS=$R_LIBS" \
      "PATH=$PATH" "TMPDIR=${TMPDIR:-/tmp}" \
      "OMP_PROC_BIND=${OMP_PROC_BIND:-false}" \
      "OMP_PLACES=${OMP_PLACES:-threads}" "$@"
    # (IOReadBandwidthMax takes a file, and limits the disk that holds it)
    set -- -p "IOReadBandwidthMax=$data_file ${bps}M" "$@"
    [ "$cap" -gt 0 ] && set -- -p "MemoryMax=${cap}M" -p "MemorySwapMax=0" "$@"
    exec sudo -n systemd-run --scope --quiet "$@"
  else
    # a scope of the user's systemd, which keeps the environment
    [ "$cap" -gt 0 ] && set -- -p "MemoryMax=${cap}M" -p "MemorySwapMax=0" "$@"
    exec systemd-run --user --scope --quiet "$@"
  fi
else
  echo "run_capped.sh: neither cgroup v1 as root nor systemd: cannot cap memory" >&2
  exit 1
fi

#!/bin/bash
# The part of the benchmark of gt_radmix() on data larger than the memory,
# shown in the article admixture_benchmark.Rmd (in the parent directory of
# this script), which also shows the code of the chunks of this script (the
# sections that start with a "## ---- label" line, read with
# knitr::read_chunk()).
#
# gt_radmix() is run three times on the simulated data of
# admixture_benchmark_large_data.R (a gen_tibble with a 50 GB backing file),
# with admixture_benchmark_large_fit.R: with the backing file in memory (in
# the page cache); with the memory of the process, page cache included,
# capped at 8 GiB, under a sixth of the data, reading from the disk as it
# is (an NVMe SSD); and with the same cap and the reads from the disk
# throttled to 500 MB/s. The cap is a transient scope of systemd (Linux,
# with cgroups v2), which needs no root: `systemd-run --user --scope -p
# MemoryMax=8G`; the throttle needs a scope of the system's systemd, so
# the third run uses sudo. Run it from the directory of the article
# (vignettes/articles), on an otherwise idle machine with enough memory to
# hold the data (about 60 GB), with tidygenclust installed with R CMD
# INSTALL, after `sudo -v`:
#
#   cd vignettes/articles
#   bash code/admixture_benchmark_large.sh
#
# It saves its results in data/admixture_benchmark (large_*.csv).
# LARGE_MEMORY_CAP sets another cap, e.g. to try the script on smaller data.

## ---- large-runs ----
set -euo pipefail
large_dir=${LARGE_DIR:-data/large}
out_dir=data/admixture_benchmark
mkdir -p "$out_dir"
k=5
seed=1
cores=8
memory_cap=${LARGE_MEMORY_CAP:-8G}
read_mbps=500
bk=$(readlink -f "$large_dir/sim_large.bk")
fit=(Rscript code/admixture_benchmark_large_fit.R)
runs=$out_dir/large_runs.csv
peaks=$out_dir/large_memory.csv
rm -f "$runs"
echo "run,peak_anon_gb" > "$peaks"

# the processes started by a process, those they started, and so on
descendants() {
  local c
  for c in $(pgrep -P "$1"); do
    echo "$c"
    descendants "$c"
  done
}
# run a command, and record the peak of the memory of its processes that is
# not the page cache (RssAnon, sampled every 0.2 s), in GB
with_peak() {
  local run=$1
  shift
  "$@" &
  local pid=$! peak=0 kb p
  while kill -0 "$pid" 2> /dev/null; do
    for p in $pid $(descendants "$pid"); do
      kb=$(awk '/^RssAnon:/ {print $2}' "/proc/$p/status" 2> /dev/null || true)
      if [ -n "$kb" ] && [ "$kb" -gt "$peak" ]; then peak=$kb; fi
    done
    sleep 0.2
  done
  wait "$pid"
  echo "$run,$(awk -v kb="$peak" 'BEGIN { printf "%.2f", kb * 1024 / 1e9 }')" \
    >> "$peaks"
}
# drop the pages of the backing file from the page cache (no root needed)
drop_cache() {
  sync "$bk"
  dd if="$bk" iflag=nocache count=0 status=none
}

# 1. in memory: the backing file is first read into the page cache
cat "$bk" > /dev/null
with_peak memory "${fit[@]}" memory $k $seed $cores "$runs"

# 2. the memory capped at 8 GiB, for R and the page cache together, from a
# cold page cache, on the disk as it is (an NVMe SSD)
drop_cache
with_peak capped_ssd \
  systemd-run --user --scope --quiet \
  -p MemoryMax=$memory_cap -p MemorySwapMax=0 \
  "${fit[@]}" capped_ssd $k $seed $cores "$runs"

# 3. the same, with the reads from the disk throttled to 500 MB/s: a scope
# of the system's systemd (as root), in which R runs as the user, with the
# environment it needs
drop_cache
with_peak capped_500mbps \
  sudo systemd-run --scope --quiet \
  -p MemoryMax=$memory_cap -p MemorySwapMax=0 \
  -p "IOReadBandwidthMax=$bk ${read_mbps}M" \
  sudo -u "$USER" env "PATH=$PATH" "R_LIBS=${R_LIBS:-}" \
  "LARGE_DIR=$large_dir" "OMP_PROC_BIND=${OMP_PROC_BIND:-false}" \
  "OMP_PLACES=${OMP_PLACES:-threads}" \
  "${fit[@]}" capped_500mbps $k $seed $cores "$runs"

## ---- large-session ----
# the machine, the disk that holds the data, and the size of the data
disk=$(lsblk -no pkname "$(df --output=source "$bk" | tail -1)")
{
  echo "item,value"
  echo "date,$(date +%F)"
  echo "CPU,$(grep -m1 '^model name' /proc/cpuinfo | sed 's/^.*: //')"
  echo "cores used,$cores"
  echo "RAM,$(awk '/^MemTotal/ {printf "%.0f GB", $2 / 1024^2}' /proc/meminfo)"
  echo "memory cap,$memory_cap"
  echo "read throttle,$read_mbps MB/s"
  echo "disk,$(sed 's/^ *//; s/ *$//' "/sys/block/$disk/device/model")"
  echo "backing file bytes,$(stat -c %s "$bk")"
  echo "tidygenclust,$(Rscript -e 'cat(format(packageVersion("tidygenclust")))')"
} > "$out_dir/large_session.csv"

#!/bin/bash
# The part of the benchmark of gt_radmix() on data larger than the memory,
# shown in the article admixture_benchmark.Rmd (in the parent directory of
# this script), which also shows the code of the chunks of this script (the
# sections that start with a "## ---- label" line, read with
# knitr::read_chunk()).
#
# gt_radmix() is run twice on the simulated data of
# admixture_benchmark_large_data.R (a gen_tibble with a 44 GB backing file),
# with admixture_benchmark_large_fit.R: once with the backing file in memory
# (in the page cache), and once with the memory of the process, page cache
# included, capped at 8 GiB, a fifth of the data. The cap is a transient
# scope of systemd (Linux, with cgroups v2), which needs no root:
# `systemd-run --user --scope -p MemoryMax=8G`. Run it from the directory of
# the article (vignettes/articles), on an otherwise idle machine with
# enough memory to hold the data (about 50 GB), with tidygenclust installed
# with R CMD INSTALL:
#
#   cd vignettes/articles
#   bash code/admixture_benchmark_large.sh
#
# It saves its results in data/admixture_benchmark (large_*.csv).

## ---- large-runs ----
set -euo pipefail
large_dir=${LARGE_DIR:-data/large}
out_dir=data/admixture_benchmark
mkdir -p "$out_dir"
k=5
seed=1
cores=8
memory_cap=8G
bk=$large_dir/sim_large.bk
runs=$out_dir/large_runs.csv
peaks=$out_dir/large_memory.csv
rm -f "$runs"
echo "run,peak_anon_mb" > "$peaks"

# run a command, and record the peak of the memory of its processes that is
# not the page cache (RssAnon, sampled every 0.2 s), in MB
with_peak() {
  local run=$1
  shift
  "$@" &
  local pid=$! peak=0 kb p
  while kill -0 "$pid" 2> /dev/null; do
    for p in $pid $(pgrep -P "$pid"); do
      kb=$(awk '/^RssAnon:/ {print $2}' "/proc/$p/status" 2> /dev/null || true)
      if [ -n "$kb" ] && [ "$kb" -gt "$peak" ]; then peak=$kb; fi
    done
    sleep 0.2
  done
  wait "$pid"
  echo "$run,$((peak / 1024))" >> "$peaks"
}

# in memory: the backing file is first read into the page cache
cat "$bk" > /dev/null
with_peak memory \
  Rscript code/admixture_benchmark_large_fit.R memory $k $seed $cores "$runs"

# memory capped at 8 GiB, for R and the page cache together, from a cold
# page cache (the pages of the backing file are dropped from it first)
sync "$bk"
dd if="$bk" iflag=nocache count=0 status=none
with_peak capped \
  systemd-run --user --scope --quiet \
  -p MemoryMax=$memory_cap -p MemorySwapMax=0 \
  Rscript code/admixture_benchmark_large_fit.R capped $k $seed $cores "$runs"

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
  echo "disk,$(sed 's/^ *//; s/ *$//' "/sys/block/$disk/device/model")"
  echo "backing file bytes,$(stat -c %s "$bk")"
  echo "bed file bytes,$(stat -c %s "$large_dir/sim_large.bed")"
  echo "tidygenclust,$(Rscript -e 'cat(format(packageVersion("tidygenclust")))')"
} > "$out_dir/large_session.csv"

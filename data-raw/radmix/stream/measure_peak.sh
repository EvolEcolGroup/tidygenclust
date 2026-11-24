#!/bin/sh
# The peak anonymous memory (RssAnon, sampled every 0.2 s) of an R script,
# in MB, written to <out_file>: the memory caps of the benchmarks are this
# plus the part of the file that the page cache may hold. Works on any
# Linux (no cgroups needed).
# Usage: measure_peak.sh <out_file> Rscript <script> [args...]
out=$1; shift
"$@" > /dev/null &
pid=$!
peak=0
while kill -0 $pid 2> /dev/null; do
  # (the R process is the child that Rscript execs, or Rscript itself)
  for p in $pid $(pgrep -P $pid); do
    kb=$(awk '/^RssAnon:/ {print $2}' /proc/$p/status 2> /dev/null)
    [ -n "$kb" ] && [ "$kb" -gt "$peak" ] && peak=$kb
  done
  sleep 0.2
done
echo $((peak / 1024)) > "$out"
echo "peak anonymous memory: $((peak / 1024)) MB"

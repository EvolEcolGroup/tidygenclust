#!/bin/bash
# Unattended covr timing of tidygenclust on the laptop (adapted from the
# benchmarking skill's desktop_job.sh). Sleeps until START, waits until the
# machine is idle, then runs the timing variants in balanced rounds, each run
# once (markers in $W/done), starting none after NO_START_AFTER. Results: one
# CSV per run in $W/results, the log in $W/results/status.md, the summary in
# $W/results/summary.md. Nothing is committed or posted.
set -uo pipefail
W=/media/andrea/annex/tmp/claude/tidygenclust/covrtime
START=${START:-2026-11-16 00:00}
NO_START_AFTER=${NO_START_AFTER:-2026-11-16 07:30}
IDLE_LOAD=${IDLE_LOAD:-1.0}
export NOT_CRAN=true TMPDIR=$W/rtmp
STATUS=$W/results/status.md; DONE=$W/done; LOGS=$W/logs
mkdir -p "$DONE" "$LOGS" "$TMPDIR"

say() { echo "- $(date '+%a %d %b %H:%M'): $*" | tee -a "$STATUS"; }
load() { cut -d' ' -f1-3 /proc/loadavg; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a < b) }'; }
too_late() { [ "$(date +%s)" -ge "$(date -d "$NO_START_AFTER" +%s)" ]; }
wait_idle() {  # three checks in a row, 5 minutes apart
  local ok=0
  while [ $ok -lt 3 ]; do
    too_late && return 1
    if lt "$(cut -d' ' -f2 /proc/loadavg)" "$IDLE_LOAD"; then ok=$((ok + 1)); else ok=0; fi
    [ $ok -lt 3 ] && sleep 300
  done
  return 0
}
settle() {  # the 1-min load below 1; after 30 min, wait for idle again
  local n=0
  while ! lt "$(cut -d' ' -f1 /proc/loadavg)" 1.0; do
    sleep 60; n=$((n + 1))
    if [ $n -ge 30 ]; then
      say "load still $(load) after 30 min: waiting for idle"
      wait_idle || return 1; n=0
    fi
  done
}
# run <variant> <omp threads> <round>
run() {
  local v=$1 omp=$2 r=$3 label="${1}_omp${2}_r${3}"
  [ -e "$DONE/$label" ] && return 0
  if too_late; then say "not started: $label (after $NO_START_AFTER)"; return 1; fi
  settle || { say "not started: $label (machine busy)"; return 1; }
  say "start $label (load $(load))"
  local t0=$SECONDS rc
  case $v in
    covrO0) OMP_NUM_THREADS=$omp Rscript "$W/run_covr.R" "$label" O0 ;;
    covrO1) OMP_NUM_THREADS=$omp Rscript "$W/run_covr.R" "$label" O1 ;;
    *) OMP_NUM_THREADS=$omp TGC_LIB=$W/lib_$v TGC_TESTS=$W/src/tests/testthat \
         TGC_OUT=$W/results/$label.csv Rscript "$W/run_files.R" ;;
  esac > "$LOGS/$label.log" 2>&1
  rc=$?
  if [ $rc -eq 0 ] && [ -s "$W/results/$label.csv" ]; then
    touch "$DONE/$label"
    say "done $label ($(( (SECONDS - t0) / 60 )) min; load after $(load))"
  else
    say "FAILED $label (exit $rc, $(( (SECONDS - t0) / 60 )) min), see logs/$label.log"
  fi
}

finish() {
  trap - EXIT
  Rscript "$W/summarise.R" >> "$LOGS/summarise.log" 2>&1 ||
    say "summary failed, see logs/summarise.log"
  say "finished"
}
trap finish EXIT

say "job started (commit $(cat "$W/results/commit.txt")); sleeping until $START"
sleep $(( $(date -d "$START" +%s) - $(date +%s) )) 2>/dev/null
say "waiting for the machine to be idle (5-min load < $IDLE_LOAD)"
if ! wait_idle; then say "the machine was never idle before $NO_START_AFTER: nothing run"; exit 0; fi
say "machine idle (load $(load))"

# the six set-ups, in balanced rounds (forward, reverse, forward); the most
# informative comparisons (covr at -O0 and -O1) come first in round 1
SETUPS=("covrO0 4" "covrO1 4" "covrO1 1" "O0 4" "O0cov 4" "O2 4")
for r in 1 2 3; do
  if [ $r -eq 2 ]; then order=(5 4 3 2 1 0); else order=(0 1 2 3 4 5); fi
  for i in "${order[@]}"; do
    run ${SETUPS[$i]} $r
    sleep 30
  done
done

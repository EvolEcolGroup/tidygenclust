#!/bin/bash
# Driver of the desktop benchmarks of streaming in gt_radmix() (local branch
# ada196_article, dev_notes/radmix_streaming.md), started by Andrea in a
# terminal after `sudo -v` (the throttled runs need root):
#
#   tmux new -s stream
#   sudo -v && bash ~/project_temp/claude/tidygenclust/stream/driver.sh
#
# - keeps sudo alive (sudo -n true every 2 minutes, from this terminal);
# - waits until the machine is idle, and lets the load settle before each
#   phase;
# - runs the phases in order, each the script $W/phases/<name>.sh, read when
#   the phase starts (a phase whose script is not there yet is waited for);
#   a marker in $W/done lets a rerun continue where it stopped;
# - stops before the next phase if $W/stop exists;
# - logs every step, with the load, to $W/status.md.
set -uo pipefail
W=${W:-$HOME/project_temp/claude/tidygenclust/stream}
PHASES=${PHASES:-"gate throttle_check ooc_third ooc passes large_data large_fits"}
IDLE_LOAD=${IDLE_LOAD:-1.0}
IDLE_CHECKS=${IDLE_CHECKS:-3}
unset OMP_NUM_THREADS OMP_THREAD_LIMIT
STATUS=$W/status.md; DONE=$W/done; LOGS=$W/logs
mkdir -p "$DONE" "$LOGS" "$W/phases" "$W/results"

say() { echo "- $(date '+%a %d %b %H:%M'): $*" | tee -a "$STATUS"; }
load() { cut -d' ' -f1-3 /proc/loadavg; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a < b) }'; }

# sudo, kept alive from this terminal (sudo's timestamps are per terminal)
if ! sudo -n true 2> /dev/null; then
  echo "run 'sudo -v' first, in this terminal"
  exit 1
fi
( while sleep 120; do sudo -n true || exit; done ) &
KEEP=$!
trap 'kill $KEEP 2> /dev/null; say "driver stopped"' EXIT

wait_idle() {  # the 5-min load below IDLE_LOAD for IDLE_CHECKS checks, 5 min apart
  local ok=0
  while [ $ok -lt $IDLE_CHECKS ]; do
    [ -e "$W/stop" ] && return 1
    if lt "$(cut -d' ' -f2 /proc/loadavg)" "$IDLE_LOAD"; then ok=$((ok + 1)); else ok=0; fi
    [ $ok -lt $IDLE_CHECKS ] && sleep 300
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

say "driver started; waiting for the machine to be idle (5-min load < $IDLE_LOAD)"
wait_idle || exit 0
say "machine idle (load $(load))"
for name in $PHASES; do
  [ -e "$DONE/$name" ] && continue
  # (the script of a phase may still be on its way)
  while [ ! -e "$W/phases/$name.sh" ]; do
    [ -e "$W/stop" ] && break 2
    sleep 60
  done
  [ -e "$W/stop" ] && { say "stop file: not started $name"; break; }
  sudo -n true 2> /dev/null || { say "sudo lapsed: not started $name"; break; }
  settle || { say "not started: $name (machine busy)"; break; }
  say "start $name (load $(load))"
  t0=$SECONDS
  if (cd "$W" && bash "$W/phases/$name.sh") >> "$LOGS/phase_$name.log" 2>&1; then
    touch "$DONE/$name"; say "done $name ($(( (SECONDS - t0) / 60 )) min)"
  else
    say "FAILED $name ($(( (SECONDS - t0) / 60 )) min), see logs/phase_$name.log"
  fi
done
say "all phases run"

#!/bin/bash
# Unattended desktop run (manica01) before the merge of tidygenclust PR #36
# (claude_jobs#3, tasks 1 and 2 of dev_notes/handover.md on newspeed_port):
#   1. port_timing: whole gt_radmix() fits of newspeed_port against dev
#      (a20eb08), with the kit of the tag archive/newspeed (worktree
#      ~/git/tidygenclust_newspeed): HGDP K = 5, 10 and EUR K = 3, 4, on 1 and
#      8 threads, 3 alternated rounds; about 30 min.
#   2. article_cli, article_r: the full benchmark of
#      vignettes/articles/admixture_benchmark.Rmd on hgdp_bench, with
#      newspeed_port (the -O2 build of phase 1) first on R's library path.
# Waits for the machine to be idle (another session's cisgem benchmark runs
# until about 19:45 on 16 Nov), starts no phase after NO_START_AFTER, logs to
# $WORK/status.md, and posts the log on claude_jobs#3 at the end. The results
# are not committed: a session reads them, then commits them to newspeed_port.
#
#   systemd-run --user --on-calendar='2026-11-16 19:45' --unit=tgc-port-bench \
#     bash ~/tgc_bench_article/run_desktop.sh
set -uo pipefail

# ---- CONFIG ------------------------------------------------------------------
JOB=tgc-port-bench
REPO=$HOME/git/tidygenclust                      # on newspeed_port
KITREPO=$HOME/git/tidygenclust_newspeed          # worktree with the kit
WORK=$HOME/tgc_bench_article
ISSUE=dramanica/claude_jobs/issues/3
IDLE_LOAD=${IDLE_LOAD:-1.0}
NO_START_AFTER=${NO_START_AFTER:-2026-11-17 08:00}
export TGC_WORK=$HOME/tgc_newspeed_work          # the kit's builds and data
PORTLIB=$TGC_WORK/libs/port
# ------------------------------------------------------------------------------

STATUS=$WORK/status.md; DONE=$WORK/done; LOGS=$WORK/logs
mkdir -p "$DONE" "$LOGS" "$WORK/results"
exec 9> "$WORK/lock"
flock -n 9 || { echo "already running"; exit 0; }

say() { echo "- $(date '+%a %d %b %H:%M'): $*" | tee -a "$STATUS"; }
load() { cut -d' ' -f1-3 /proc/loadavg; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a < b) }'; }
too_late() { [ "$(date +%s)" -ge "$(date -d "$NO_START_AFTER" +%s)" ]; }

wait_idle() {
  local ok=0
  while [ $ok -lt 3 ]; do
    too_late && return 1
    if lt "$(cut -d' ' -f2 /proc/loadavg)" "$IDLE_LOAD"; then ok=$((ok + 1)); else ok=0; fi
    [ $ok -lt 3 ] && sleep 300
  done
  return 0
}
settle() {
  local n=0
  while ! lt "$(cut -d' ' -f1 /proc/loadavg)" 1.0; do
    sleep 60; n=$((n + 1))
    if [ $n -ge 30 ]; then
      say "load still $(load) after 30 min: waiting for idle"
      wait_idle || return 1; n=0
    fi
  done
}
phase() {  # phase <name> <dir> <command...>
  local name=$1 dir=$2; shift 2
  [ -e "$DONE/$name" ] && return 0
  if too_late; then say "not started: $name (after $NO_START_AFTER)"; return 1; fi
  settle || { say "not started: $name (machine busy)"; return 1; }
  say "start $name (load $(load))"
  local t0=$SECONDS
  if (cd "$dir" && "$@") >> "$LOGS/phase_$name.log" 2>&1; then
    touch "$DONE/$name"; say "done $name ($(( (SECONDS - t0) / 60 )) min)"
  else
    say "FAILED $name ($(( (SECONDS - t0) / 60 )) min), see logs/phase_$name.log"
  fi
}

report() {
  trap - EXIT
  say "reporting"
  {
    echo "## $(date +%F): unattended desktop run of $JOB finished (manica01)"
    echo
    echo "Results on manica01 in \`$WORK\` (port timing: \`results/port_timing/summary.md\`) and \`vignettes/articles/data/admixture_benchmark/\`; not committed or interpreted yet."
    echo
    cat "$STATUS"
  } > "$WORK/comment.md"
  gh api "repos/$ISSUE/comments" -F body=@"$WORK/comment.md" >> "$LOGS/report.log" 2>&1 ||
    say "could not post the comment; it is in $WORK/comment.md"
}

# the article's programs: conda env ctidygenclust (ADMIXTURE, fastmixture),
# ADAMIXTURE from its own env, newspeed_port first on R's library path
article_env() {
  source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate ctidygenclust || return 1
  export R_LIBS=$PORTLIB
  export ADAMIXTURE=$HOME/miniconda3/envs/cadamixture/bin/adamixture
  export ADAMIXTURE_PYTHON=$HOME/miniconda3/envs/cadamixture/bin/python
  export OPENADMIXTURE_JULIA_PREFIX=$HOME/miniconda3/envs/openadmixturer_julia
}
article_cli() { article_env && git -C "$REPO" log -1 --format='newspeed_port %h' &&
  taskset -c 0-7 bash code/admixture_benchmark_cli.sh; }
article_r() { article_env && taskset -c 0-7 Rscript code/admixture_benchmark_radmix.R; }

trap report EXIT
# another session's timed cisgem benchmark: its last phase is single-threaded
# with pauses, so the load alone does not show it; wait for its unit to end
say "scheduled run started; waiting for cisgem_final.service to end"
while systemctl --user is-active --quiet cisgem_final.service ||
      [ "$(systemctl --user is-active cisgem_final.service)" = activating ]; do
  too_late && { say "cisgem_final.service still running at $NO_START_AFTER: nothing run"; exit 0; }
  sleep 120
done
say "cisgem_final.service not running; waiting for the machine to be idle (5-min load < $IDLE_LOAD)"
if ! wait_idle; then say "the machine was never idle before $NO_START_AFTER: nothing run"; exit 0; fi
say "machine idle (load $(load))"

phase port_timing "$KITREPO" env CASES=hgdp:5,hgdp:10,eur:3,eur:4 THREADS="1 8" \
  ROUNDS=3 MAX_ITER=1000 bash data-raw/radmix/newspeed/run_bench.sh \
  "$WORK/results/port_timing" dev port
phase article_cli "$REPO/vignettes/articles" article_cli
phase article_r "$REPO/vignettes/articles" article_r
say "all phases run"

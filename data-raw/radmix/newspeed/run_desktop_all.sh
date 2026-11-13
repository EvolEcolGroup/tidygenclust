#!/bin/bash
# Runs every timed desktop check of HANDOVER.md, unattended, in the order of
# priority, then commits the summaries to the branches `newspeed` and
# `supervised`, pushes them, and reports on the claude_jobs issue.
#
# Scheduled for Sunday 2026-11-15 09:00 on manica01 (systemd user timer
# tgc-desktop-checks). The builds (build_variants.sh, labels below) and the
# runs that do not depend on timing (CV trajectories, test suites, thread
# invariance) were done beforehand on the E-cores.
#
# - Waits until the machine is idle (5-min load average below IDLE_LOAD for
#   three checks in a row, 5 minutes apart) before starting, and before each
#   phase waits for the load left by the previous one to settle.
# - Each phase leaves a marker in $TGC_WORK/done, so the script can be run
#   again after an interruption and goes on where it stopped.
# - Does not start a phase after NO_START_AFTER (the user is back on Monday
#   morning; timings with the machine in use are worthless).
#
#   bash data-raw/radmix/newspeed/run_desktop_all.sh
set -uo pipefail
export TGC_WORK=${TGC_WORK:-$HOME/tgc_newspeed_work}
NS=${NS:-$HOME/git/tidygenclust_newspeed}         # worktree of newspeed
SUP=${SUP:-$HOME/git/tidygenclust_supervised}     # worktree of supervised
export TGC_HGDP_DIR=${TGC_HGDP_DIR:-$HOME/git/tidygenclust/vignettes/articles/data/hgdp}
export HGDP_DIR=$TGC_HGDP_DIR
export ALGO_SCRATCH=$TGC_WORK/algo
SUP_SCRATCH=$TGC_WORK/sup
IDLE_LOAD=${IDLE_LOAD:-1.0}
NO_START_AFTER=${NO_START_AFTER:-2026-11-16 07:00}
ISSUE=${ISSUE:-dramanica/claude_jobs/issues/3}
STATUS=$TGC_WORK/status.md
DONE=$TGC_WORK/done
K=data-raw/radmix/newspeed
mkdir -p "$DONE" "$ALGO_SCRATCH" "$SUP_SCRATCH"

say() { echo "- $(date '+%a %d %b %H:%M'): $*" | tee -a "$STATUS"; }
load5() { cut -d' ' -f2 /proc/loadavg; }
load1() { cut -d' ' -f1 /proc/loadavg; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a < b) }'; }
too_late() { [ "$(date +%s)" -ge "$(date -d "$NO_START_AFTER" +%s)" ]; }

# the machine idle: three checks in a row, 5 minutes apart
wait_idle() {
  local ok=0
  while [ $ok -lt 3 ]; do
    too_late && return 1
    if lt "$(load5)" "$IDLE_LOAD"; then ok=$((ok + 1)); else ok=0; fi
    [ $ok -lt 3 ] && sleep 300
  done
  return 0
}
# between phases: the load of the previous phase gone (1-min average)
settle() {
  local n=0
  while ! lt "$(load1)" 1.0; do
    sleep 60; n=$((n + 1))
    if [ $n -ge 30 ]; then
      say "load still $(cat /proc/loadavg | cut -d' ' -f1-3) after 30 min: waiting for idle"
      wait_idle || return 1
      n=0
    fi
  done
}

# phase <name> <command...>: run once (marker), from the newspeed worktree
phase() {
  local name=$1; shift
  [ -e "$DONE/$name" ] && return 0
  if too_late; then say "not started: $name (after $NO_START_AFTER)"; return 1; fi
  settle || { say "not started: $name (machine busy)"; return 1; }
  say "start $name (load $(cut -d' ' -f1-3 /proc/loadavg))"
  local t0=$SECONDS
  if (cd "$NS" && "$@") >> "$TGC_WORK/logs/phase_$name.log" 2>&1; then
    touch "$DONE/$name"
    say "done $name ($(( (SECONDS - t0) / 60 )) min)"
  else
    say "FAILED $name ($(( (SECONDS - t0) / 60 )) min), see logs/phase_$name.log"
  fi
}

bench() { local out=$1; shift; bash "$K/run_bench.sh" "$TGC_WORK/results/$out" "$@"; }
algo() { taskset -c 0-7 env LIB="$TGC_WORK/libs/head" NC=8 bash data-raw/radmix/newspeed_algo/run_desktop.sh "$@"; }
sup_timings() {
  cd "$SUP" || return 1
  bash data-raw/radmix/supervised/run_hgdp_timings.sh "$TGC_WORK/libs/sup" "$SUP_SCRATCH" 2 3 1 &&
  bash data-raw/radmix/supervised/run_hgdp_timings.sh "$TGC_WORK/libs/sup" "$SUP_SCRATCH" 0-7 3 8
}
sup_kernels() {
  cd "$SUP" && OMP_NUM_THREADS=1 taskset -c 2 Rscript data-raw/radmix/supervised/time_kernels.R \
    "$TGC_WORK/libs/sup" "$SUP_SCRATCH" 5
}
# ADMIXTURE --supervised -j8 against gt_radmix(n_cores = 8) (the fits of
# sup_timings), timed
sup_admixture() {
  cd "$SUP_SCRATCH/admix" || return 1
  for case in regions pops; do
    for ext in bed bim fam; do ln -sf "$TGC_HGDP_DIR/hgdp_bench.$ext" "hgdp_$case.$ext"; done
    if [ ! -s "hgdp_$case.5.Q" ]; then
      uptime
      /usr/bin/time -v taskset -c 0-7 admixture --supervised -s 1 "hgdp_$case.bed" 5 -j8 \
        > "hgdp_$case.5.log" 2> "hgdp_$case.5.time" || return 1
      grep -E "Elapsed|Maximum resident" "hgdp_$case.5.time"
    fi
    (cd "$SUP" && Rscript data-raw/radmix/supervised/compare_hgdp_admixture.R "$SUP_SCRATCH" "$case" \
      "$SUP_SCRATCH/hgdp_${case}_random_1_sup.rds") || return 1
  done
}

# copies the summaries into the branches, commits, pushes, and reports
report() {
  trap - EXIT
  say "reporting"
  local d=$NS/$K/desktop_results
  mkdir -p "$d"
  for r in "$TGC_WORK"/results/*/; do
    [ -d "$r" ] || continue
    local p; p=$(basename "$r"); mkdir -p "$d/$p"
    cp "$r"/summary.md "$r"/*.csv "$r"/run.log "$d/$p/" 2>/dev/null
  done
  local a=$NS/data-raw/radmix/newspeed_algo/results_desktop
  mkdir -p "$a"
  cp "$ALGO_SCRATCH"/out/*.csv "$ALGO_SCRATCH"/out/*.txt "$a/" 2>/dev/null
  cp -r "$ALGO_SCRATCH"/out/cv_policies* "$a/" 2>/dev/null
  cp "$ALGO_SCRATCH"/logs/*.log "$a/" 2>/dev/null
  cp "$STATUS" "$TGC_WORK"/logs/phase_*.log "$d/" 2>/dev/null
  local s=$SUP/data-raw/radmix/supervised/results_desktop
  mkdir -p "$s"
  cp "$SUP_SCRATCH"/*.csv "$s/" 2>/dev/null
  cp "$SUP_SCRATCH"/admix/*.log "$SUP_SCRATCH"/admix/*.time "$s/" 2>/dev/null
  cp "$TGC_WORK"/logs/phase_sup_*.log "$s/" 2>/dev/null
  local msg="Desktop checks of HANDOVER.md (unattended run of $(date +%F))

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
  local pushed=""
  for wt in "$NS" "$SUP"; do
    if (cd "$wt" && git add -A data-raw/radmix && git commit -qm "$msg" && git push -q origin HEAD) \
      >> "$TGC_WORK/logs/report.log" 2>&1; then
      pushed="$pushed $(cd "$wt" && git rev-parse --abbrev-ref HEAD)@$(cd "$wt" && git rev-parse --short HEAD)"
    fi
  done
  local body=$TGC_WORK/comment.md
  {
    echo "**Unattended desktop run (manica01), finished $(date '+%a %d %b %Y %H:%M')**"
    echo
    echo "Pushed:${pushed:- nothing (see ~/tgc_newspeed_work/logs/report.log)}. The tables are in \`data-raw/radmix/newspeed/desktop_results/<phase>/summary.md\` and \`data-raw/radmix/newspeed_algo/results_desktop/\` on \`newspeed\`, and in \`data-raw/radmix/supervised/results_desktop/\` on \`supervised\`. Everything, with the raw outputs, is in \`~/tgc_newspeed_work\` on the desktop."
    echo
    echo "Log of the run:"
    echo
    cat "$STATUS"
    echo
    echo "Not interpreted yet: the decisions of HANDOVER.md still need a read of the summaries."
  } > "$body"
  gh api "repos/$ISSUE/comments" -F body=@"$body" >> "$TGC_WORK/logs/report.log" 2>&1 \
    || say "could not post the comment; it is in $body"
}

mkdir -p "$TGC_WORK/logs"
trap report EXIT
say "scheduled run started; waiting for the machine to be idle (5-min load < $IDLE_LOAD)"
# the CV trajectories (machine-independent) may still run on the E-cores
if systemctl --user is-active -q tgc-cvtraj; then
  systemctl --user stop tgc-cvtraj
  say "stopped the CV trajectories still running on the E-cores (restartable)"
fi
if ! wait_idle; then
  say "the machine was never idle before $NO_START_AFTER: nothing run"
  exit 0
fi
say "machine idle (load $(cut -d' ' -f1-3 /proc/loadavg))"

# 1. per-iteration kit (newspeed/README.md, "Commands")
phase kernels env MODE=kernels CASES=hgdp:5,hgdp:10,eur:4 THREADS=1 ROUNDS=3 \
  bash "$K/run_bench.sh" "$TGC_WORK/results/kernels" base fused inline final final_llnoinl
phase iter8 env MAX_ITER=8 CASES=hgdp:5,hgdp:10,hgdp_full:5 THREADS="1 4 8" ROUNDS=3 \
  bash "$K/run_bench.sh" "$TGC_WORK/results/iter8" base fused final head
phase whole env CASES=hgdp:5,hgdp:10,hgdp_full:5,eur:3,eur:4 THREADS="1 4 8" ROUNDS=3 \
  bash "$K/run_bench.sh" "$TGC_WORK/results/whole" base final head
phase seeds env SEEDS=1,2,3 CASES=hgdp:5,hgdp:10,eur:3,eur:4 THREADS=8 ROUNDS=1 \
  bash "$K/run_bench.sh" "$TGC_WORK/results/seeds" base final
phase cvkit env MODE=cv CASES=hgdp:5,hgdp_full:5 THREADS="1 8" ROUNDS=3 \
  bash "$K/run_bench.sh" "$TGC_WORK/results/cv" base cvmask final
# 3. CV policy, wall-clock; 4. warm-ups; 5. supervised
phase cv_timing algo cv_timing
phase warmup_hgdp algo warmup hgdp
phase warmup_eur algo warmup eur
phase sup_timings sup_timings
phase sup_kernels sup_kernels
phase sup_admixture sup_admixture
# 2. tiles and schedule ("only if step 1 leaves time"): after the rest
phase tiles env MAX_ITER=8 CASES=hgdp:5,hgdp:10,hgdp_full:5 THREADS="1 8" ROUNDS=3 \
  bash "$K/run_bench.sh" "$TGC_WORK/results/tiles" final final_q64 final_q128 final_guided final_tj32
# the rest of the CV trajectories, if the E-core run did not finish them
# (machine-independent; the script skips the runs already done), then the
# summary tables of newspeed_algo
phase cv_traj algo cv
phase algo_summary algo summary
say "all phases run"

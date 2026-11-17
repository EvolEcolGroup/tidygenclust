#!/bin/bash
# Focused single runs, after the overnight job: covr with the comparison with
# ADMIXTURE skipped under covr (abe2d2a) and the fold tests on 2 threads
# (99eef49), at covr's default -O0 (as in CI on dev) and at -O1 with bounds
# checks (PR #37). One run each: the effects looked for are 5 vs 50 min.
W=/media/andrea/annex/tmp/claude/tidygenclust/covrtime
export NOT_CRAN=true TMPDIR=$W/rtmp OMP_NUM_THREADS=4
STATUS=$W/results/status.md
say() { echo "- $(date '+%a %d %b %H:%M'): $*" | tee -a "$STATUS"; }
load() { cut -d' ' -f1-3 /proc/loadavg; }
# wait for the overnight job to end
while systemctl --user is-active --quiet tgc-covr-timing; do sleep 60; done
# the two test changes, in the source copy
T=$W/src/tests/testthat
grep -q "skip_on_covr" $T/test_radmix_vs_admixture.R ||
  sed -i 's/^skip_on_cran()$/skip_on_cran()\nskip_on_covr()/' $T/test_radmix_vs_admixture.R
sed -i 's/n_cores = 3)/n_cores = 2)/' $T/test_radmix_cv.R
grep -n "skip_on" $T/test_radmix_vs_admixture.R | tee -a $W/logs/focus.log
grep -n "n_cores = [23])" $T/test_radmix_cv.R | tee -a $W/logs/focus.log
for v in "covrO0_skip O0" "covrO1_skip O1"; do
  set -- $v
  sleep 30
  say "start $1 (load $(load))"
  t0=$SECONDS
  Rscript $W/run_covr.R $1 $2 > $W/logs/$1.log 2>&1
  say "done $1 (exit $?, $(( (SECONDS - t0) / 60 )) min; $(grep -h 'coverage ' $W/logs/$1.log); load after $(load))"
done

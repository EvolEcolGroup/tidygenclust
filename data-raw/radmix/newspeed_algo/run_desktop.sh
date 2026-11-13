#!/usr/bin/env bash
# The runs of the newspeed_algo experiments for the desktop (i9-13900, 8
# P-cores): wall-clock confirmation of the laptop's pass counts, and the
# runs too heavy for the laptop. See README.md for what each step answers.
#
# Usage (from any folder):
#   LIB=<library with the tidygenclust build> \
#   ALGO_SCRATCH=<scratch folder, several GB free> \
#   HGDP_DIR=<folder with hgdp_bench.bed, hgdp_full.bed, hgdp_bench_quick.bed> \
#   [EUR_BED=<path to EUR_subset.bed>] [NC=8] \
#   bash run_desktop.sh <step> [<dataset>]
# Steps: smoke (HGDP quick, a few minutes), cv, cv_timing, warmup, summary.
# Every script skips the runs whose output already exists, so a step can be
# interrupted and restarted.
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
: "${LIB:?set LIB to the library holding the tidygenclust build}"
: "${ALGO_SCRATCH:?set ALGO_SCRATCH to a scratch folder}"
: "${HGDP_DIR:?set HGDP_DIR to the folder of the HGDP .bed files}"
export ALGO_SCRATCH HGDP_DIR
NC=${NC:-8}
# threads bound to the P-cores (OpenMP threads drift onto the E-cores
# otherwise), R's temporary files in the scratch folder
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=$NC
export TMPDIR=$ALGO_SCRATCH/rtmp
mkdir -p "$TMPDIR" "$ALGO_SCRATCH/out" "$ALGO_SCRATCH/logs"
out=$ALGO_SCRATCH/out
log=$ALGO_SCRATCH/logs
run() { echo "== $*"; Rscript "$@"; }

case "${1:-}" in
  smoke)
    # everything once on HGDP quick (toy data: checks that the scripts run,
    # not their conclusions)
    run "$here/check_replica.R" "$LIB" hgdp_quick 3 "$NC"
    run "$here/cv_trajectories.R" "$LIB" hgdp_quick 2:3 1 "$NC" "$out/smoke_cv"
    run "$here/cv_policies.R" "$out/smoke_cv"
    run "$here/cv_policy_timing.R" "$LIB" hgdp_quick 3 1 conv,it5,rtol1e-6 \
      "$NC" "$out/smoke_cv_timing.csv"
    run "$here/warmup_compare.R" "$LIB" hgdp_quick 4 1:2 \
      random,svd,sub10,mbtol3 "$NC" "$out/smoke_warmup"
    run "$here/warmup_summary.R" "$out/smoke_warmup"
    ;;
  cv)
    # Part A trajectories: HGDP K = 2..10 (seeds 1-3), EUR K = 1..6 (seeds
    # 1-5); the laptop did HGDP seed 1 for some K and EUR seed 1 only
    run "$here/cv_trajectories.R" "$LIB" hgdp 2:10 1 "$NC" "$out/cv" \
      2>&1 | tee -a "$log/cv_hgdp.log"
    run "$here/cv_trajectories.R" "$LIB" eur 1:6 1:5 "$NC" "$out/cv" \
      2>&1 | tee -a "$log/cv_eur.log"
    run "$here/cv_trajectories.R" "$LIB" hgdp 2:10 2:3 "$NC" "$out/cv" \
      2>&1 | tee -a "$log/cv_hgdp.log"
    ;;
  cv_timing)
    # Part A wall-clock: the package's radmix_cv() with each policy
    run "$here/cv_policy_timing.R" "$LIB" eur 2:6 1 \
      conv,it5,it10,rtol1e-6,rtol1e-5 "$NC" "$out/cv_timing_eur.csv"
    run "$here/cv_policy_timing.R" "$LIB" hgdp 2:8 1 \
      conv,it5,it10,rtol1e-6,rtol1e-5 "$NC" "$out/cv_timing_hgdp.csv"
    ;;
  warmup)
    # Part B: starts and warm-ups, seeds 1-5 (K = 10 and the unpruned HGDP
    # were too heavy for the laptop)
    ds=${2:-all}
    if [[ $ds == all || $ds == eur ]]; then
      for k in 3 4; do
        run "$here/warmup_compare.R" "$LIB" eur $k 1:5 \
          random,svd,sub10,sub25,mbtol3 "$NC" "$out/warmup"
      done
    fi
    if [[ $ds == all || $ds == hgdp ]]; then
      for k in 5 10; do
        run "$here/warmup_compare.R" "$LIB" hgdp $k 1:5 \
          random,svd,sub5,sub10,sub25,mbtol3 "$NC" "$out/warmup"
      done
    fi
    if [[ $ds == all || $ds == hgdp_full ]]; then
      for k in 5 10; do
        run "$here/warmup_compare.R" "$LIB" hgdp_full $k 1:5 \
          random,svd,sub2,sub5,sub10,mbtol3 "$NC" "$out/warmup"
      done
    fi
    ;;
  summary)
    run "$here/cv_policies.R" "$out/cv" "$out/cv_policies" \
      "$here/../../../tests/testthat/fixtures/admixture_eur_reference.rds"
    run "$here/warmup_summary.R" "$out/warmup" "$out/warmup.csv"
    ;;
  *)
    echo "usage: bash run_desktop.sh smoke|cv|cv_timing|warmup|summary [dataset]" >&2
    exit 1
    ;;
esac

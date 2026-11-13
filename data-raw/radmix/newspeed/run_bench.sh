#!/bin/bash
# Times gt_radmix() fits of several builds (made by build_variants.sh),
# alternated (A B A B ...) for each case and number of threads, with pauses,
# threads pinned, and the load recorded with each fit; then summarises them.
#
#   bash data-raw/radmix/newspeed/run_bench.sh <out dir> <label> <label> ...
#
# The first label is the reference of the summary. Environment variables
# (defaults in brackets):
#   MODE       fits: whole fits or fixed iterations (bench_fits.R, the
#              default); kernels: single passes (bench_kernels.R, KREPS
#              timings each [5]); cv: fold fits of MAX_ITER iterations, plain
#              view against the training view with and without the bit mask
#              (bench_cv.R, KREPS timings each [3]; MAX_ITER defaults to 4)
#   TGC_WORK   the work directory of build_variants.sh [../tgc_newspeed_work]
#   CASES      dataset:K pairs [hgdp:5,hgdp:10,hgdp_full:5,eur:3,eur:4]
#   THREADS    numbers of threads [1 4 8]
#   ROUNDS     repeats of each build, alternated [3]
#   SEEDS      seeds of the random starts, e.g. 1,2,3 [1]
#   MAX_ITER   iteration limit [1000: whole fits; e.g. 8 gives the cost of a
#              fixed number of iterations, the same work for every build]
#   RTOL       convergence threshold [1e-8]
#   PIN1       CPU of the 1-thread runs [2]
#   PCORES     CPUs of the multi-thread runs, the P-cores [0-7, those of the
#              i9-13900 with SMT off]; threads are also bound with
#              OMP_PROC_BIND=close OMP_PLACES=cores
#   PAUSE      seconds between fits [5]
# Results: <out dir>/fits.csv (one row per fit), <out dir>/fits/*.rds
# (log-likelihood histories and Q, for the comparisons), <out dir>/run.log
# (uptime before each fit, and errors), <out dir>/summary.md.
set -uo pipefail
REPO=$(git rev-parse --show-toplevel)
KIT=$REPO/data-raw/radmix/newspeed
export TGC_WORK=${TGC_WORK:-$(dirname "$REPO")/tgc_newspeed_work}
[ $# -ge 2 ] || { echo "usage: run_bench.sh <out dir> <label> ..."; exit 1; }
OUT=$1; shift
LABELS="$*"
CASES=${CASES:-hgdp:5,hgdp:10,hgdp_full:5,eur:3,eur:4}
THREADS=${THREADS:-1 4 8}
ROUNDS=${ROUNDS:-3}
SEEDS=${SEEDS:-1}
MODE=${MODE:-fits}
if [ "$MODE" = cv ]; then MAX_ITER=${MAX_ITER:-4}; else MAX_ITER=${MAX_ITER:-1000}; fi
RTOL=${RTOL:-1e-8}
PIN1=${PIN1:-2}
PCORES=${PCORES:-0-7}
export PAUSE=${PAUSE:-5}
mkdir -p "$OUT"
cd "$REPO"
for label in $LABELS; do
  [ -d "$TGC_WORK/libs/$label/tidygenclust" ] || { echo "no build $label in $TGC_WORK/libs"; exit 1; }
done
{ echo "== $(date) run_bench.sh $OUT $LABELS"; echo "MODE=$MODE CASES=$CASES THREADS=$THREADS ROUNDS=$ROUNDS SEEDS=$SEEDS MAX_ITER=$MAX_ITER RTOL=$RTOL PIN1=$PIN1 PCORES=$PCORES"; } >> "$OUT/run.log"
# the datasets are converted once (common.R caches them in $TGC_WORK/data)
Rscript -e "source('$KIT/common.R'); for (cs in tgc_cases('$CASES')) invisible(tgc_read(cs\$data))" >> "$OUT/run.log" 2>&1
for n in $THREADS; do
  if [ "$n" -eq 1 ]; then pre="taskset -c $PIN1"; else pre="taskset -c $PCORES"; fi
  for cs in ${CASES//,/ }; do
    for r in $(seq "$ROUNDS"); do
      for label in $LABELS; do
        echo "$(date +%T) $label $cs n=$n round $r: $(uptime)" >> "$OUT/run.log"
        case $MODE in
          fits) cmd=("$KIT/bench_fits.R" "$TGC_WORK/libs/$label" "$label" "$OUT" "$cs" "$SEEDS" "$n" "$MAX_ITER" "$RTOL") ;;
          kernels) cmd=("$KIT/bench_kernels.R" "$TGC_WORK/libs/$label" "$label" "$OUT" "$cs" "$n" "${KREPS:-5}") ;;
          cv) cmd=("$KIT/bench_cv.R" "$TGC_WORK/libs/$label" "$label" "$OUT" "$cs" "$n" "$MAX_ITER" "${KREPS:-3}") ;;
          *) echo "unknown MODE $MODE"; exit 1 ;;
        esac
        OMP_NUM_THREADS=$n OMP_PROC_BIND=close OMP_PLACES=cores OPENBLAS_NUM_THREADS=1 \
          $pre Rscript "${cmd[@]}" >> "$OUT/run.log" 2>&1 \
          || echo "FAILED: $label $cs n=$n (see $OUT/run.log)" | tee -a "$OUT/run.log"
      done
    done
  done
done
Rscript "$KIT/summarise.R" "$OUT" > "$OUT/summary.md" 2>> "$OUT/run.log"
echo "finished: $OUT/summary.md"

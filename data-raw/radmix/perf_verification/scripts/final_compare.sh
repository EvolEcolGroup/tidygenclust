#!/bin/bash
# The final comparison of the laptop, repeated: the baseline, the combination
# and the combination without multiversioning (or other builds), alternated,
# with pauses: kernels (1 core), fits on 1 and more cores, with
# cross-validation, and at rtol = 1e-8. Summary in <out dir>/summary.txt.
#
#   bash final_compare.sh <out dir>
#
# Environment variables: TGC_WORK, PAUSE (default 10), LABELS (default
# "base comb comb_nomv"; the first is the reference), NCORES (core counts of
# the multi-core fits, default "4 16"), ROUNDS (fit rounds, default 2),
# KERNEL_CASES and FIT_CASES (dataset:K cases; defaults: those of the
# laptop; FIT_CASES also applies to the CV fits), PIN (a CPU to pin the
# 1-core runs to with taskset, e.g. PIN=2) and PIN_MULTI (CPUs to confine
# the multi-core runs to, e.g. PIN_MULTI=0-7 for the P-cores of a hybrid
# processor); by default nothing is pinned.
set -uo pipefail
REPO=$(git rev-parse --show-toplevel)
SCRIPTS=$REPO/data-raw/radmix/perf_verification/scripts
WORK=${TGC_WORK:-$(dirname "$REPO")/tgc_perf_work}
PAUSE=${PAUSE:-10}
OUT=$1
B=${LABELS:-base comb comb_nomv}
read -r REF COMB _ <<< "$B"
NCORES=${NCORES:-4 16}
ROUNDS=${ROUNDS:-2}
RUN1=""; [ -n "${PIN:-}" ] && RUN1="taskset -c $PIN"
RUNN=""; [ -n "${PIN_MULTI:-}" ] && RUNN="taskset -c $PIN_MULTI"
mkdir -p "$OUT"
cd "$REPO"
run() {
  local pre=$RUN1
  [ "$2" -gt 1 ] && pre=$RUNN
  $pre Rscript "$SCRIPTS/bench_fits.R" "$WORK/libs/$1" "$1" "$OUT" "$2" "$3" "$4" "${FIT_CASES:-}" \
    >> "$OUT/run.log" 2>&1 || echo "fits $* failed (see $OUT/run.log)"
  sleep "$PAUSE"
}
for r in 1 2 3; do for v in $B; do
  $RUN1 Rscript "$SCRIPTS/bench_kernels.R" "$WORK/libs/$v" "$v" "$OUT" 1 "${KERNEL_CASES:-}" \
    >> "$OUT/run.log" 2>&1 || echo "kernels $v failed (see $OUT/run.log)"
  sleep "$PAUSE"
done; done
for r in $(seq "$ROUNDS"); do for v in $B; do run "$v" 1 nocv 1e-7; done; done
for n in $NCORES; do
  for r in $(seq "$ROUNDS"); do for v in $B; do run "$v" "$n" nocv 1e-7; done; done
done
for r in $(seq "$ROUNDS"); do for v in $B; do run "$v" 1 cv 1e-7; done; done
for r in $(seq "$ROUNDS"); do for v in $REF $COMB; do run "$v" 1 nocv 1e-8; done; done
Rscript "$SCRIPTS/final_summary.R" "$OUT" "$REF" > "$OUT/summary.txt" 2>&1
echo "finished $OUT"

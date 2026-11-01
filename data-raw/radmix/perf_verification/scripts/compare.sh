#!/bin/bash
# Compares builds of tidygenclust: runs the kernel and fit benchmarks of each
# build in turn (A B C A B C ...), with pauses, then compares each build
# with the first one (summary_<label>.txt in the output directory).
#
#   bash compare.sh <out dir> <kernel rounds> <fit rounds> <label A> <label B> [...]
#
# Run from the root of the repository. Environment variables:
# - TGC_WORK: where the builds are (libs/<label>, see build_variants.sh)
# - PAUSE: seconds between runs (default 10)
# - PIN: a CPU to pin these 1-core runs to with taskset (e.g. PIN=2;
#   default: not pinned)
# - KERNEL_CASES, FIT_CASES: the dataset:K cases (defaults: those of the
#   scripts, as on the laptop)
set -uo pipefail
REPO=$(git rev-parse --show-toplevel)
SCRIPTS=$REPO/data-raw/radmix/perf_verification/scripts
WORK=${TGC_WORK:-$(dirname "$REPO")/tgc_perf_work}
PAUSE=${PAUSE:-10}
OUT=$1; KR=$2; FR=$3; shift 3; LABELS="$*"; A=$1
RUN=""
[ -n "${PIN:-}" ] && RUN="taskset -c $PIN"
mkdir -p "$OUT"
cd "$REPO"
for r in $(seq "$KR"); do for v in $LABELS; do
  $RUN Rscript "$SCRIPTS/bench_kernels.R" "$WORK/libs/$v" "$v" "$OUT" 1 \
    "${KERNEL_CASES:-eur:4,hgdp_quick:5,hgdp_quick:10}" \
    >> "$OUT/run.log" 2>&1 || echo "kernels $v failed (see $OUT/run.log)"
  sleep "$PAUSE"
done; done
for r in $(seq "$FR"); do for v in $LABELS; do
  $RUN Rscript "$SCRIPTS/bench_fits.R" "$WORK/libs/$v" "$v" "$OUT" 1 nocv 1e-7 \
    "${FIT_CASES:-eur:3,eur:4,hgdp_quick:5,hgdp_quick:7}" \
    >> "$OUT/run.log" 2>&1 || echo "fits $v failed (see $OUT/run.log)"
  sleep "$PAUSE"
done; done
shift
for B in "$@"; do
  Rscript "$SCRIPTS/summarise.R" "$OUT" "$A" "$B" > "$OUT/summary_$B.txt" 2>&1
done
echo "finished $OUT"

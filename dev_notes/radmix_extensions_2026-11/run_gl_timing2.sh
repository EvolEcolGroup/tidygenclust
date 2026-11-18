#!/bin/bash
# The definitive timing of the sources of genotype likelihoods (quiet
# machine): passes and whole fits of gt_radmix() on the simulated gl_tibble
# (400 x 100,000, K = 3) for
#   tidygenclust: ext (radmix_extensions 7561833) and inl (+ slice_gl() and
#     loglik_gl_slice() forced inline into the AVX2 clones, tgc_inline.patch)
#   tidypopgl source: base (a87e576), tab (runs read directly, one lookup in
#     a 64K table per code, branchless indicator of data), br (tab with the
#     two-table decode), tabidx (tab without the run path), simd (tab with
#     an AVX2 clone and omp simd: GCC does not vectorise it), avx (tab with
#     AVX2 gathers for the runs, gl_radmix_source_avx2.cpp)
# in a rotating order (reversed every other round); the in-memory example
# source (dense) with ext + base. 1 thread on CPU 2 (3 rounds), then 8 on
# the P-cores (5 rounds).
#   bash run_gl_timing2.sh <out csv> <rounds 1 thread> <rounds 8 threads>
set -u
W=/home/andrea/tgc_supfused_work
OUT=$1; ROUNDS=${2:-3}; ROUNDS8=${3:-5}
combos=("ext base" "inl base" "inl tab" "inl br" "inl tabidx" "inl simd" "inl avx" "ext avx")
one() {  # one <tgc> <gl> <threads> [extra]
  local t=$3
  if [ "$t" = 1 ]; then
    OMP_NUM_THREADS=1 taskset -c 2 Rscript $W/time_gl_passes.R $W/lib_$1 $W/lib_gl_$2 "$1_$2" 1 5 "$OUT" fit ${4:-} > /dev/null 2>> $W/logs/gl_timing_err.log
  else
    OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores taskset -c 0-7 \
      Rscript $W/time_gl_passes.R $W/lib_$1 $W/lib_gl_$2 "$1_$2" $t 5 "$OUT" fit ${4:-} > /dev/null 2>> $W/logs/gl_timing_err.log
  fi
  sleep 2
}
echo "start $(date) load $(cut -d' ' -f1-3 /proc/loadavg)"
for t in 1 8; do
  R=$ROUNDS; [ "$t" = 8 ] && R=$ROUNDS8
  for r in $(seq 1 $R); do
    # rotate the order of the combinations, reversed every other round
    order=()
    for i in 0 1 2 3 4 5 6 7; do order+=("${combos[$(( (i + 3 * (r - 1)) % 8 ))]}"); done
    if [ $((r % 2)) = 0 ]; then
      rev=(); for (( i=${#order[@]}-1; i>=0; i-- )); do rev+=("${order[$i]}"); done; order=("${rev[@]}")
    fi
    for c in "${order[@]}"; do
      set -- $c
      extra=""; [ "$1 $2" = "ext base" ] && extra=dense
      one $1 $2 $t $extra
    done
    echo "threads $t round $r done $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)"
  done
done
echo "end $(date)"

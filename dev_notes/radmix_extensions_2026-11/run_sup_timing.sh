#!/bin/bash
# Whole supervised gt_radmix() fits of the pruned HGDP (K = 5), builds of
# tidygenclust alternated in balanced order: dev (00deb4f), new (the fused
# pass of supervised fits, branch supervised_fused) and allB (the
# alternative: dev's fused pass over every individual, then the references
# restored; allB.patch). Cases pops (139 references) and regions (660), and
# unsupervised as a control (dev and new run the same code there).
# 1 thread pinned to CPU 2; 8 threads on the P-cores 0-7.
#   bash run_sup_timing.sh <out csv> <rounds_1> <rounds_8>
set -u
W=/home/andrea/tgc_supfused_work
OUT=$1; R1=${2:-3}; R8=${3:-5}
mkdir -p $W/rds
orders=("dev new allB" "allB new dev" "new dev allB" "allB dev new" "dev allB new")
fit() {  # fit <build> <case> <threads>
  local b=$1 c=$2 t=$3
  if [ "$t" = 1 ]; then
    OMP_NUM_THREADS=1 taskset -c 2 Rscript $W/fit_sup.R $W/lib_$b $c 1 1 "$OUT" $W/rds > /dev/null 2>> $W/logs/timing_err.log
  else
    OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores taskset -c 0-7 \
      Rscript $W/fit_sup.R $W/lib_$b $c $t 1 "$OUT" $W/rds > /dev/null 2>> $W/logs/timing_err.log
  fi
  sleep 2
}
echo "start $(date) load $(cut -d' ' -f1-3 /proc/loadavg)"
for t in 8 1; do
  rounds=$R8; [ "$t" = 1 ] && rounds=$R1
  cases="pops regions unsupervised"; [ "$t" = 1 ] && cases="pops regions"
  for r in $(seq 1 $rounds); do
    for c in $cases; do
      for b in ${orders[$(( (r - 1) % 5 ))]}; do
        [ "$c" = unsupervised ] && [ "$b" = allB ] && continue
        fit $b $c $t
      done
    done
    echo "threads $t round $r done $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)"
  done
done
echo "end $(date)"

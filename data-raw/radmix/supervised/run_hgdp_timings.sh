#!/bin/bash
# Timings of supervised and unsupervised fits on HGDP pruned, K = 5, one
# core, alternating the cases (U R P U R P), with hgdp_supervised.R.
# Usage, from the root of the package:
#   bash data-raw/radmix/supervised/run_hgdp_timings.sh <library> \
#     <scratch dir> <cores> [repeats] [n_cores]
# <cores> is a taskset CPU list (e.g. 2, or 0-7 for the P-cores of the
# desktop); n_cores (default 1) is the number of threads of each fit.
lib=$1
scratch=$2
core=$3
reps=${4:-2}
nc=${5:-1}
export OMP_NUM_THREADS=$nc OMP_PROC_BIND=close OMP_PLACES=cores
export TMPDIR=$scratch/rtmp
mkdir -p "$TMPDIR"
for rep in $(seq 1 "$reps"); do
  for case in unsupervised regions pops; do
    uptime
    taskset -c "$core" Rscript data-raw/radmix/supervised/hgdp_supervised.R \
      "$lib" "$scratch" "$case" random "$nc" 1
  done
done

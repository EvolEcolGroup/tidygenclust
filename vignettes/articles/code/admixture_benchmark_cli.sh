#!/usr/bin/env bash
# Benchmark of ADMIXTURE and fastmixture, run from the command line, on the
# HGDP data prepared by admixture_benchmark_data.R (in the same directory as
# this script). Together with admixture_benchmark_radmix.R, which runs
# gt_radmix() with the same settings, it gives the results shown in
# the article admixture_benchmark.Rmd (in the parent directory of this
# script), which also shows the code of the chunks of this script (the
# sections that start with a "## ---- label" line, read with
# knitr::read_chunk()).
#
# Run from the directory of the article (vignettes/articles), with ADMIXTURE
# and fastmixture available, e.g. in the conda environment created by
# tidygenclust::tgc_tools_install() (the programs can also be given with the
# environment variables ADMIXTURE and FASTMIXTURE):
#
#   cd vignettes/articles
#   conda activate ctidygenclust
#   bash code/admixture_benchmark_cli.sh
#
# Run it on an otherwise idle machine with at least 8 cores. It saves the Q
# matrix, log and timing of each run in data/admixture_benchmark, which the
# article reads (the other outputs, such as the P matrices, are left in
# data/hgdp/runs). `bash code/admixture_benchmark_cli.sh quick` runs
# everything on a small subset of the data, to check the script, and saves
# its results in data/admixture_benchmark_quick instead.

## ---- cli-setup ----
set -euo pipefail
if [[ "${1:-}" == "quick" ]]; then
  bed=data/hgdp/hgdp_bench_quick.bed
  out_dir=data/admixture_benchmark_quick
  runs_dir=data/hgdp/runs_quick
else
  bed=data/hgdp/hgdp_bench.bed
  out_dir=data/admixture_benchmark
  runs_dir=data/hgdp/runs
fi
admixture=${ADMIXTURE:-admixture}
fastmixture=${FASTMIXTURE:-fastmixture}

## ---- cli-settings ----
# the same settings as in admixture_benchmark_radmix.R
k=5
seed=1
n_folds=5
cores=(1 4 8)

## ---- cli-functions ----
# the wall-clock time, in seconds
now() { date +%s.%N; }

# the first match of a regular expression (with one group, which is
# returned) at the start of a line of a log, or NA if there is none
extract() {
  local value
  value=$(sed -n -E "s/$1/\1/p" "$2" | head -n 1)
  echo "${value:-NA}"
}

# ADMIXTURE writes its output in the current directory, named after the bed
# file (e.g. hgdp_bench.5.Q), so it is run from the directory of the run;
# its log is what it prints
run_admixture() {
  local n_cores=$1 crossval=$2 run_dir=$3
  local bed_path=$PWD/$bed
  local args=(-s "$seed" -j"$n_cores")
  if [[ $crossval == TRUE ]]; then
    args+=(--cv="$n_folds")
  fi
  (cd "$run_dir" && "$admixture" "${args[@]}" "$bed_path" "$k" > admixture.log 2>&1)
}

# fastmixture writes its output (Q, P and a log with its settings and
# results) with the prefix given by --out, and prints its progress
run_fastmixture() {
  local n_cores=$1 crossval=$2 run_dir=$3
  local args=(
    --bfile "${bed%.bed}" --K "$k" --threads "$n_cores" --seed "$seed"
    --out "$run_dir/fastmixture"
  )
  if [[ $crossval == TRUE ]]; then
    args+=(--cv "$n_folds")
  fi
  "$fastmixture" "${args[@]}" > "$run_dir/console.log" 2>&1
}

## ---- cli-runs ----
# each combination of program, number of cores and cross-validation (or
# not) is run once, in its own directory; the time includes reading the
# data. The Q matrix and log of each run are then copied to the results,
# and its log-likelihood, cross-validation error and number of iterations
# are extracted from its log
mkdir -p "$out_dir/Q" "$out_dir/logs"
timings=$out_dir/timings_cli.csv
echo "tool,n_cores,crossval,seconds,loglik,cv_error,iterations" > "$timings"
for tool in admixture fastmixture; do
  for n_cores in "${cores[@]}"; do
    for crossval in FALSE TRUE; do
      if [[ $crossval == TRUE ]]; then
        run=${tool}_t${n_cores}_cv
      else
        run=${tool}_t${n_cores}_nocv
      fi
      run_dir=$runs_dir/$run
      rm -rf "$run_dir"
      mkdir -p "$run_dir"
      echo "Running $run"
      start=$(now)
      "run_$tool" "$n_cores" "$crossval" "$run_dir"
      seconds=$(awk -v start="$start" -v end="$(now)" \
        'BEGIN { printf "%.3f", end - start }')
      case $tool in
        admixture)
          q_file=$run_dir/$(basename "$bed" .bed).$k.Q
          log_file=$run_dir/admixture.log
          loglik=$(extract '^Loglikelihood: ([^[:space:]]+).*' "$log_file")
          cv_error=$(extract \
            '^CV error \(K=[0-9]+\): ([^[:space:]]+).*' "$log_file")
          ;;
        fastmixture)
          q_file=$run_dir/fastmixture.K$k.s$seed.Q
          log_file=$run_dir/fastmixture.K$k.s$seed.log
          loglik=$(extract \
            '^Final log-likelihood: ([^[:space:]]+).*' "$log_file")
          cv_error=$(extract \
            '^Cross-validation error \(SD\): ([^[:space:]]+).*' "$log_file")
          ;;
      esac
      iterations=$(extract '^Converged in ([0-9]+) iterations.*' "$log_file")
      cp "$q_file" "$out_dir/Q/$run.Q"
      cp "$log_file" "$out_dir/logs/$run.log"
      echo "$tool,$n_cores,$crossval,$seconds,$loglik,$cv_error,$iterations" \
        >> "$timings"
    done
  done
done

## ---- cli-session ----
# the machine and the versions of the programs (the CPU and memory are read
# from /proc, so they are only available on Linux)
row() { printf '"%s","%s"\n' "$1" "$2"; }
cpu=$(sed -n -E 's/^model name[[:space:]]*: //p' /proc/cpuinfo 2> /dev/null |
  head -n 1 || true)
ram=$(awk '/^MemTotal/ { printf "%.0f GB", $2 / 1024^2 }' /proc/meminfo \
  2> /dev/null || true)
admixture_version=$("$admixture" 2>&1 |
  sed -n -E 's/.*ADMIXTURE Version ([^[:space:]]+).*/\1/p' | head -n 1 || true)
fastmixture_version=$("$fastmixture" --version 2>&1 | sed 's/^v//' || true)
{
  row item value
  row date "$(date +%F)"
  row CPU "${cpu:-NA}"
  row cores "$(getconf _NPROCESSORS_ONLN)"
  row RAM "${ram:-NA}"
  row ADMIXTURE "${admixture_version:-NA}"
  row fastmixture "${fastmixture_version:-NA}"
  row individuals "$(wc -l < "${bed%.bed}.fam" | tr -d ' ')"
  row loci "$(wc -l < "${bed%.bed}.bim" | tr -d ' ')"
} > "$out_dir/session_cli.csv"
cat "$timings"
echo "Results saved in $out_dir"

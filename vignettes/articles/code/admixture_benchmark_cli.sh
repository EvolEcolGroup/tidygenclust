#!/usr/bin/env bash
# Benchmark of the programs run from the command line, on the HGDP data
# prepared by admixture_benchmark_data.R (in the same directory as this
# script): ADMIXTURE, fastmixture, ADAMIXTURE, gt_radmix() (with Rscript,
# running admixture_benchmark_radmix_cli.R) and OpenADMIXTURE.jl (with julia,
# running admixture_benchmark_openadmixture.jl). Together with
# admixture_benchmark_radmix.R, which runs gt_radmix() in an R session with
# the same settings, it gives the results shown in the article
# admixture_benchmark.Rmd (in the parent directory of this script), which
# also shows the code of the chunks of this script (the sections that start
# with a "## ---- label" line, read with knitr::read_chunk()).
#
# Run from the directory of the article (vignettes/articles), with ADMIXTURE,
# fastmixture and ADAMIXTURE available, e.g. in the conda environment created
# by tidygenclust::tgc_tools_install() (the programs can also be given with
# the environment variables ADMIXTURE, FASTMIXTURE and ADAMIXTURE, and the
# python of ADAMIXTURE, for its version, with ADAMIXTURE_PYTHON), tidygenclust
# installed
# (Rscript can be given with RSCRIPT), and OpenADMIXTURE.jl installed in the
# conda environment created by data-raw/radmix/setup_julia_env.sh (from the
# root of the package; another environment can be given with
# OPENADMIXTURE_JULIA_ENV, or its path with OPENADMIXTURE_JULIA_PREFIX):
#
#   cd vignettes/articles
#   conda activate ctidygenclust
#   bash code/admixture_benchmark_cli.sh
#
# Run it on an otherwise idle machine with at least as many cores as the
# largest number in cores (below). It saves the Q matrix, log and timing of
# each run in data/admixture_benchmark, which the article reads (the other
# outputs, such as the P matrices, are left in data/hgdp/runs). `bash code/admixture_benchmark_cli.sh quick` runs
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
adamixture=${ADAMIXTURE:-adamixture}
adamixture_python=${ADAMIXTURE_PYTHON:-python}
rscript=${RSCRIPT:-Rscript}
julia_env=${OPENADMIXTURE_JULIA_ENV:-openadmixturer_julia}
julia_prefix=${OPENADMIXTURE_JULIA_PREFIX:-$(conda info --base)/envs/$julia_env}
# the data on which OpenADMIXTURE.jl is run first, so that Julia compiles its
# code before the fit is timed (the small subset of the data)
warmup_bed=data/hgdp/hgdp_bench_quick.bed

## ---- cli-settings ----
# the same settings as in admixture_benchmark_radmix.R
k=5
seed=1
n_folds=5
cores=(1 2 4)

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

# ADAMIXTURE writes its output (Q, P and a plot, which cannot be turned off)
# in --save_dir, named after --name and K (e.g. adamixture.5.Q), and prints
# its progress. Its command line drops the chromosomes that are not
# autosomes by default (hence --chrom_mode all), and -t alone does not limit
# its number of threads (it sets OMP_NUM_THREADS after loading its OpenMP
# library), so OMP_NUM_THREADS is set before starting it
run_adamixture() {
  local n_cores=$1 crossval=$2 run_dir=$3
  local args=(
    --data_path "$bed" --save_dir "$run_dir" --name adamixture -k "$k"
    -t "$n_cores" -s "$seed" --chrom_mode all
  )
  if [[ $crossval == TRUE ]]; then
    args+=(--cv "$n_folds")
  fi
  OMP_NUM_THREADS=$n_cores "$adamixture" "${args[@]}" \
    > "$run_dir/console.log" 2>&1
}

# gt_radmix() is run by a new R process, which loads tidygenclust, reads the
# bed file into a gen_tibble, runs gt_radmix(), saves the results, and
# prints the log-likelihood and cross-validation error
run_radmix_cli() {
  local n_cores=$1 crossval=$2 run_dir=$3
  "$rscript" code/admixture_benchmark_radmix_cli.R "$bed" "$k" "$seed" \
    "$n_cores" "$crossval" "$n_folds" "$run_dir" > "$run_dir/console.log" 2>&1
}

# julia, with the variables that activating its conda environment would
# set, so that it finds the packages installed in it (without the time
# taken by conda run)
julia_run() {
  env JULIA_DEPOT_PATH="$julia_prefix/share/julia:" \
    JULIA_PROJECT="@$julia_env" JULIA_LOAD_PATH="@:@$julia_env:@stdlib" \
    "$julia_prefix/bin/julia" "$@"
}

# OpenADMIXTURE.jl is run by a new Julia process, which writes Q and P with
# the prefix openadmixture, and prints its progress, the log-likelihood and
# the time of the fit (Julia is started without an interactive thread, see
# admixture_benchmark_openadmixture.jl); it has no cross-validation
run_openadmixture_cli() {
  local n_cores=$1 run_dir=$3
  julia_run --threads="$n_cores,0" code/admixture_benchmark_openadmixture.jl \
    "$bed" "$k" "$seed" "$run_dir/openadmixture" > "$run_dir/console.log" 2>&1
}

# the same, after a warm-up on the small subset of the data, so that the
# fit, timed in Julia, does not include compiling the code
run_openadmixture() {
  local n_cores=$1 run_dir=$3
  julia_run --threads="$n_cores,0" code/admixture_benchmark_openadmixture.jl \
    "$bed" "$k" "$seed" "$run_dir/openadmixture" "$warmup_bed" \
    > "$run_dir/console.log" 2>&1
}

## ---- cli-runs ----
# each combination of program, number of cores and cross-validation (or
# not) is run once, in its own directory; the time includes starting the
# program, reading the data and saving the results. OpenADMIXTURE.jl, which
# has no cross-validation, is run without it, both from the command line
# (openadmixture_cli) and after a warm-up (openadmixture, timed in Julia).
# The Q matrix and log of each run are then copied to the results, and its
# log-likelihood, cross-validation error and number of iterations are
# extracted from its log
mkdir -p "$out_dir/Q" "$out_dir/logs"
timings=$out_dir/timings_cli.csv
echo "tool,n_cores,crossval,seconds,loglik,cv_error,iterations" > "$timings"
for tool in admixture fastmixture adamixture radmix_cli openadmixture_cli \
  openadmixture; do
  for n_cores in "${cores[@]}"; do
    for crossval in FALSE TRUE; do
      if [[ $crossval == TRUE && $tool == openadmixture* ]]; then
        continue
      fi
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
      iterations_pattern='^Converged in ([0-9]+) iterations.*'
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
        adamixture)
          q_file=$run_dir/adamixture.$k.Q
          log_file=$run_dir/console.log
          loglik=$(extract \
            '^ *Final log-likelihood: ([^[:space:]]+).*' "$log_file")
          cv_error=$(extract \
            '^ *K=[0-9]+: CV index = ([^[:space:]]+).*' "$log_file")
          # (the lines indented further are those of the initialisation)
          iterations_pattern='^ {4}Converged at iteration ([0-9]+).*'
          ;;
        radmix_cli)
          q_file=$run_dir/radmix.Q
          log_file=$run_dir/radmix.log
          loglik=$(extract '^Loglikelihood: ([^[:space:]]+).*' \
            "$run_dir/console.log")
          cv_error=$(extract '^CV error \(K=[0-9]+\): ([^[:space:]]+).*' \
            "$run_dir/console.log")
          ;;
        openadmixture*)
          q_file=$run_dir/openadmixture.Q
          log_file=$run_dir/console.log
          loglik=$(extract '^Loglikelihood: ([^[:space:]]+).*' "$log_file")
          cv_error=NA
          iterations_pattern='^Main algorithm converged in ([0-9]+) iterations.*'
          if [[ $tool == openadmixture ]]; then
            # after the warm-up, the time of the fit alone
            seconds=$(extract '^Fit time \(s\): ([^[:space:]]+).*' "$log_file")
          fi
          ;;
      esac
      iterations=$(extract "$iterations_pattern" "$log_file")
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
# (the command line of ADAMIXTURE has no --version)
adamixture_version=$("$adamixture_python" -c \
  'import adamixture; print(adamixture.__version__)' 2>&1 || true)
julia_version=$(julia_run --startup-file=no -e 'println(VERSION)' || true)
openadmixture_version=$(julia_run --startup-file=no -e '
  import Pkg
  for p in values(Pkg.dependencies())
    p.name == "OpenADMIXTURE" &&
      println(p.version, " (", something(p.git_revision, "?"), ")")
  end' || true)
{
  row item value
  row date "$(date +%F)"
  row CPU "${cpu:-NA}"
  row cores "$(getconf _NPROCESSORS_ONLN)"
  row RAM "${ram:-NA}"
  row ADMIXTURE "${admixture_version:-NA}"
  row fastmixture "${fastmixture_version:-NA}"
  row ADAMIXTURE "${adamixture_version:-NA}"
  row Julia "${julia_version:-NA}"
  row OpenADMIXTURE.jl "${openadmixture_version:-NA}"
  row individuals "$(wc -l < "${bed%.bed}.fam" | tr -d ' ')"
  row loci "$(wc -l < "${bed%.bed}.bim" | tr -d ' ')"
} > "$out_dir/session_cli.csv"
cat "$timings"
echo "Results saved in $out_dir"

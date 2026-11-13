# Fewer passes over the genotypes: CV refit policies and warm-ups

Experiments of the night of 2026-11-13 on the algorithm of `gt_radmix()`
(not its kernels): how many passes over the genotypes the cross-validation
refits and the main fit need, on large data. This file says what each
script does and how to run them on the desktop; the laptop's results are
summarised in `results_laptop/` (CSV files written by `cv_policies.R` and
`warmup_summary.R`).

Measures are machine-independent where possible: iterations, passes over
the genotypes (a pass over a subset of the loci counts as its fraction of a
full pass), log-likelihoods and CV errors, which are bit-identical for any
number of threads. Times on the laptop are indicative only: the desktop
runs below give the wall-clock times.

## Scripts

| script | what it does |
|---|---|
| `algo_lib.R` | shared helpers: loading the datasets (gen_tibbles created once in `$ALGO_SCRATCH/data`), pass counting from the log of a fit, an R replica of the fitting loop that records every iterate, the subset and mini-batch warm-ups |
| `check_replica.R` | checks that the replica gives bit-identical iterates to `radmix_fit()` (main fit and a CV refit) |
| `cv_trajectories.R` | Part A: for each K, the main fit as `gt_radmix(seed)` runs it, then each fold refitted to convergence from the full-data estimates, scored on its held-out genotypes after every iteration |
| `cv_policies.R` | Part A summary: the CV error, the K chosen and the cost of each refit policy (converged, at most m iterations, a looser rtol), read from the trajectories |
| `cv_policy_timing.R` | Part A wall-clock: the package's `radmix_cv()` run with each policy, timed |
| `warmup_compare.R` | Part B: the default random start, the SVD start, the subset warm-ups and the mini-batch warm-up, each followed by the full fit at rtol 1e-8 |
| `warmup_summary.R` | Part B summary: passes, iterations and log-likelihoods (distance from the best mode found) per method |
| `run_desktop.sh` | the desktop runs, by step |

All scripts take the library holding the tidygenclust build as their first
argument (load an installed `-O2` build: never `load_all()`), and skip the
runs whose output exists, so a step can be interrupted and restarted.

## Running on the desktop

Install the build to test (here the branch merged into `newspeed`, or
`dev`) into its own library, with optimisation, and check the flags:

```bash
git archive <ref> | tar -x -C /path/to/src_algo
R CMD INSTALL --preclean --library=/path/to/lib_algo /path/to/src_algo
# the log must show -O2 for src/radmix_*.cpp
```

Then, from the repository root (paths are examples):

```bash
export LIB=/path/to/lib_algo
export ALGO_SCRATCH=/path/to/scratch/algo      # several GB free
export HGDP_DIR=/path/to/vignettes/articles/data/hgdp
export EUR_BED=$HOME/.cache/R/tidygenclust/EUR_subset.bed   # if elsewhere
S=data-raw/radmix/newspeed_algo/run_desktop.sh
bash $S smoke          # HGDP quick, a few minutes: the scripts run
bash $S cv_timing      # Part A wall-clock (rough guess: ~1 h)
bash $S warmup hgdp    # Part B on the pruned HGDP, K = 5 and 10 (~1 h)
bash $S warmup hgdp_full   # Part B on the unpruned HGDP (~2-3 h)
bash $S warmup eur     # Part B on EUR_subset, 10 seeds: the modes (~20 min)
bash $S cv             # Part A, more K and seeds (several hours)
bash $S summary        # the summary tables, from all the outputs
```

`smoke` was run on the laptop (HGDP quick, 1 core) on 2026-11-13 and
completed. The durations above are guesses from the laptop's pass counts
and the desktop's earlier timings.

## The laptop's results

`results_laptop/` holds the summaries of the laptop runs of 2026-11-13
(one thread per job, on a machine loaded by other jobs, so the times in
them are indicative only):

- `cv_policies_errors.csv` and `cv_policies.txt`: Part A, HGDP K = 2, 3,
  5, 8 (seed 1) and EUR_subset K = 1-5 (seed 1), every policy;
- `warmup.csv` and `warmup.txt`: Part B, HGDP K = 5 (seeds 1-3) and
  K = 10 (seed 1), unpruned HGDP K = 5 (seed 1), EUR_subset K = 3 and 4
  (seeds 1-5);
- `smoke_cv_timing.csv`: the smoke run of `cv_policy_timing.R` (HGDP
  quick, K = 3).

`NC` (default 8) sets the number of threads; they are bound to the
P-cores with `OMP_PROC_BIND=close OMP_PLACES=cores`. Keep the machine
otherwise idle during `cv_timing` and `warmup` (their times are the
point); `cv` and `summary` only produce machine-independent numbers.
Outputs go to `$ALGO_SCRATCH/out`: one `.rds` per run, and the CSV
summaries (`cv_policies_errors.csv`, `warmup.csv`, `cv_timing_*.csv`).
The laptop's outputs can be copied into the same folders: the desktop runs
then skip them (their passes and log-likelihoods are the same on both
machines, but their times are not comparable).

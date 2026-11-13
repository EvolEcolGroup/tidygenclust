# newspeed: per-iteration speed of `gt_radmix()` on large data

Scripts to build variants of tidygenclust from commits and to time them on
large data (HGDP pruned and unpruned, EUR_subset), alternated, pinned,
with the load recorded, and summarised as markdown tables. They were
written for the desktop (i9-13900, 8 P-cores = CPUs 0-7, SMT off) and
tested on the laptop with a tiny configuration. What was done and measured
on the laptop is in `dev_notes/radmix_performance.md` (the changes of
2026-11-13), with the background.

## Files

| file | what it does |
|---|---|
| `common.R` | shared helpers: datasets, loading a build, counts from the log |
| `build_variants.sh` | `label=commit` arguments: `git archive` of each commit, `R CMD INSTALL --preclean` into `$TGC_WORK/libs/<label>`, check of `-O2` and of the AVX2 clones |
| `bench_fits.R` | fits (random start, seeds given), whole or of a fixed number of iterations; `fits.csv` + log-likelihood histories and Q in `fits/*.rds` |
| `bench_kernels.R` | single passes (loglik, update_q, update_p, update_q_loglik) from the same (q, p); `kernels.csv` |
| `bench_cv.R` | fold fits of a fixed number of iterations: plain view, training view hashing every genotype, training view with the bit mask; `cv.csv` |
| `run_bench.sh` | runs one of the three for several builds, alternated, pinned, with the load; then `summarise.R` |
| `summarise.R` | markdown tables (`summary.md`): medians, ratios to the reference build, iterations, passes, log-likelihood differences, thread invariance |

## Data

`common.R` reads `vignettes/articles/data/hgdp/hgdp_bench.bed` (`hgdp`, 1043
x 130,943), `hgdp_full.bed` (`hgdp_full`, 1043 x 642,794),
`hgdp_bench_quick.bed` (`hgdp_quick`, toy data) and the EUR_subset of the
tests (`eur`, 379 x 54,051). Set `TGC_HGDP_DIR` if the HGDP files are
elsewhere, `TGC_EUR_BED` for EUR. To make them: the data script of the
benchmark article (`vignettes/articles/code/admixture_benchmark_data.R`) for
`hgdp` and `hgdp_quick`, `git show
archive/radmix-perf-verify:data-raw/radmix/perf_verification/scripts/prepare_data.R`
for `hgdp_full` and `eur`. The first use of a dataset converts it into a
backing file under `$TGC_WORK/data` (about 0.7 GB for `hgdp_full`); later
runs reload it.

## Protocol

- Only builds made by `build_variants.sh` (fresh sources of a commit,
  `--preclean`, `-O2` checked); never a `load_all()` build.
- An otherwise idle machine; the 1-thread runs on one P-core (`PIN1`, CPU 2
  by default), the multi-thread runs confined to the P-cores (`PCORES`, 0-7)
  with `OMP_PROC_BIND=close OMP_PLACES=cores` (without confinement, OpenMP
  threads drift onto the E-cores).
- Builds alternated (A B C A B C) for each case and number of threads,
  medians of at least 3; check in `run.log` and in the load columns that
  nothing else ran; check with the repeats that the noise is below the
  effect.
- Machine-independent checks come with the times: identical iteration
  counts, log-likelihoods equal up to rounding (relative differences of
  ~1e-12; the fused pass adds the same terms in another order), identical
  results for any number of threads, passes over the genotypes per fit.

## Commands (desktop)

From the root of the repository, on branch `newspeed` (or any branch, after
`git fetch`, as long as the commits are known). The work directory holds
the builds, the cached data and nothing else:

```bash
export TGC_WORK=$HOME/tgc_newspeed_work
K=data-raw/radmix/newspeed
# 1. the builds: the baseline (origin/dev), each change of newspeed, and the
#    tip of newspeed, `head`, which adds the options of the algorithmic work
#    (cv_max_iter/cv_rtol, init = "subset", all opt-in, so its default fits
#    should match `final`); add as many label=commit as needed (the commits
#    of newspeed: git log --oneline b5999cc..newspeed)
MAKEFLAGS=-j8 PIN=0-7 bash $K/build_variants.sh base=b5999cc \
  fused=378d830 cvmask=e2e3661 inline=d358c02 final=13f89be \
  head=newspeed
# variants of the final code with compiler flags (a separate call each, as
# EXTRA_FLAGS applies to all the builds of a call): the log-likelihood of the
# fused pass called in baseline code; more Q-side tasks; other schedules
EXTRA_FLAGS=-DRADMIX_FUSED_LL_NOINLINE MAKEFLAGS=-j8 PIN=0-7 \
  bash $K/build_variants.sh final_llnoinl=13f89be
EXTRA_FLAGS=-DRADMIX_MIN_Q_TASKS=64 MAKEFLAGS=-j8 PIN=0-7 \
  bash $K/build_variants.sh final_q64=13f89be
EXTRA_FLAGS=-DRADMIX_MIN_Q_TASKS=128 MAKEFLAGS=-j8 PIN=0-7 \
  bash $K/build_variants.sh final_q128=13f89be
EXTRA_FLAGS=-DRADMIX_OMP_SCHEDULE=guided MAKEFLAGS=-j8 PIN=0-7 \
  bash $K/build_variants.sh final_guided=13f89be
EXTRA_FLAGS=-DRADMIX_TILE_J=32 MAKEFLAGS=-j8 PIN=0-7 \
  bash $K/build_variants.sh final_tj32=13f89be

# 2. single passes, 1 thread (a few minutes): the cost of each pass and of
#    an iteration made of them (does the laptop's 0.78 hold?)
MODE=kernels CASES=hgdp:5,hgdp:10,eur:4 THREADS=1 ROUNDS=3 \
  bash $K/run_bench.sh $TGC_WORK/results/kernels base fused inline final final_llnoinl

# 3. a fixed number of iterations (the same work for every build), then
#    whole fits; 1, 4 and 8 threads
MAX_ITER=8 CASES=hgdp:5,hgdp:10,hgdp_full:5 THREADS="1 4 8" ROUNDS=3 \
  bash $K/run_bench.sh $TGC_WORK/results/iter8 base fused final head
CASES=hgdp:5,hgdp:10,hgdp_full:5,eur:3,eur:4 THREADS="1 4 8" ROUNDS=3 \
  bash $K/run_bench.sh $TGC_WORK/results/whole base final head
# correctness at several seeds (iterations, log-likelihoods), 8 threads
SEEDS=1,2,3 CASES=hgdp:5,hgdp:10,eur:3,eur:4 THREADS=8 ROUNDS=1 \
  bash $K/run_bench.sh $TGC_WORK/results/seeds base final

# 4. cross-validation passes: fold fits of 4 iterations, plain view against
#    the training view with and without the bit mask (in builds without the
#    mask, both training views hash)
MODE=cv CASES=hgdp:5,hgdp_full:5 THREADS="1 8" ROUNDS=3 \
  bash $K/run_bench.sh $TGC_WORK/results/cv base cvmask final

# 5. scaling of the Q side (decisions for many cores): 8 threads (and 1, to
#    see the cost of the extra tasks); results differ by rounding with other
#    MIN_Q_TASKS, so compare times and iterations, not log-likelihoods
MAX_ITER=8 CASES=hgdp:5,hgdp:10,hgdp_full:5 THREADS="1 8" ROUNDS=3 \
  bash $K/run_bench.sh $TGC_WORK/results/tiles final final_q64 final_q128 final_guided final_tj32
```

Each `$TGC_WORK/results/<phase>/summary.md` has the tables; `fits.csv`,
`kernels.csv` and `cv.csv` the raw rows (with the load). Rough durations on
the desktop (per round and build): kernels ~1 min; iter8 at 1 thread ~1.5
min for hgdp K = 5, ~4 min for K = 10, ~6 min for hgdp_full; whole fits at 1
thread ~0.5 min (hgdp K = 5), ~2 min (K = 10), ~3 min (hgdp_full), ~5 s
(eur). Run the 1-thread phases overnight if needed (`THREADS=1`), then the
multi-thread ones.

## Tested on the laptop

The kit was used for all the laptop measurements of 2026-11-13 (builds,
`MODE=kernels`, `MODE=cv`, whole fits), and run with a tiny configuration
(toy data, 1 iteration, 1 and 2 threads), with the builds already made:

```bash
export TGC_WORK=/media/andrea/annex/tmp/claude/tidygenclust/newspeed
export TGC_HGDP_DIR=/media/andrea/annex/git/tidygenclust/vignettes/articles/data/hgdp
K=data-raw/radmix/newspeed
PIN=0 MAKEFLAGS=-j1 bash $K/build_variants.sh cvmask=<commit>
MAX_ITER=1 CASES=hgdp_quick:3 THREADS="1 2" ROUNDS=2 PIN1=0 PCORES=0 PAUSE=0 \
  bash $K/run_bench.sh $TGC_WORK/out/kit_tiny base sg5ll
MODE=kernels CASES=hgdp:5 THREADS=1 ROUNDS=2 PIN1=0 KREPS=3 \
  bash $K/run_bench.sh $TGC_WORK/out/kernels_sg5 base fused sg5 sg5ll
MODE=cv MAX_ITER=2 CASES=hgdp:5 THREADS=1 ROUNDS=1 PIN1=0 KREPS=2 \
  bash $K/run_bench.sh $TGC_WORK/out/cv base cvmask
```

## Laptop results

`laptop_results/` holds the tables of the laptop runs of 2026-11-13 (1 core,
the other three busy with other jobs, so the times are indicative). The
labels are those of builds made from the work-in-progress commits, before
the history of the branch was tidied; the code of each is that of:

| label | code |
|---|---|
| `base` | b5999cc (origin/dev) |
| `fused` | 378d830 (the fused pass) |
| `cvmask` | e2e3661 (the bit mask of cross-validation) |
| `ll2`, `ll2inline` | a first version of the cheaper log-likelihood, whose loop did not vectorise (1.5x slower; `ll2inline` was meant to inline it into the clone but did not) |
| `sg5` | 13f89be built with `-DRADMIX_FUSED_LL_NOINLINE` |
| `sg5ll` | 13f89be (the final code) |

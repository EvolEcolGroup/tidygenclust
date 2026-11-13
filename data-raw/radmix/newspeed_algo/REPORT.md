# CV refit policy and warm-ups for gt_radmix(): laptop report

Night of 2026-11-13, on the laptop (i7-10510U). The work was done by an
agent, and this report was saved by the session that coordinated it. The
summary tables of the scripts are in `results_laptop/`. The commands for the
desktop are in `README.md` and `run_desktop.sh`.

## Commits (on the branch `newspeed`, where they were copied from the temporary branch `newspeed-algo`)

- `05b5e5a` The experiment scripts: `algo_lib.R`, `check_replica.R`,
  `cv_trajectories.R`, `cv_policies.R`, `cv_policy_timing.R`,
  `warmup_compare.R`, `warmup_summary.R`, `run_desktop.sh` and `README.md`.
- `b092e3d` `gt_radmix(cv_max_iter = max_iter, cv_rtol = rtol)`: optional
  shorter refits of the CV folds. The default is unchanged. The commit
  includes the docs with the measurements, NEWS, dev_notes and tests in
  `test_gt_radmix.R`.
- `3cf356f` `gt_radmix(init = "subset", init_subset = 0.05)`: a start from
  a fit to a subset of the loci, in the new file `R/radmix_subset.R`. It is
  an option, and the default is unchanged. The commit includes the docs,
  NEWS, dev_notes and `tests/testthat/test_radmix_subset.R`. **This start
  moves runs to lower modes (see Part B), so the commit may be dropped.**
- `07ea43e` The laptop results: `results_laptop/*.csv` and `*.txt`.

## How the work was checked

- **An R replica of the fitting loop.** `fit_replica()` calls the kernels
  of the package in the same order as `radmix_fit()`, and
  `check_replica.R` shows that the two are bit-identical: the same
  log-likelihood histories, Q and pass counts, for a main fit and for a CV
  refit.
- **One trajectory gives every CV policy.** Each fold is refitted once, to
  convergence, and the held-out genotypes are scored after every
  iteration. Every policy (converged, at most m iterations, or a looser
  rtol) stops somewhere on that same deterministic path, so all the
  policies are read from one run.
- **The converged policy matches the package.** The main fit and the seeds
  of the folds are drawn exactly as in
  `gt_radmix(seed = 1, crossval = TRUE)`, so the converged policy
  reproduces its CV error.
- **Tests.** An -O2 build (checked in the install log and with `readelf`)
  passes `test_gt_radmix.R`, `test_radmix_cv.R` and `test_radmix_subset.R`.
  The new tests fail without the change (unknown argument). lintr is clean,
  apart from `object_usage` artefacts of linting single files.
- **The desktop script.** `run_desktop.sh smoke` ran to completion on HGDP
  quick with 1 core.

## Measurements

Laptop, one thread per job, on cores 2 and 3. Other agents kept the load
average at 4-14, so the times are indicative only. A "pass" is one pass over
the genotypes; a pass over a subset of the loci counts as its fraction of a
full pass. A pass with a CV fold hidden costs about 1.5x a plain pass in the
code of b5999cc (about 1x with the bit mask of the branch `newspeed`).

### Part A: the CV refit policy

Seed 1 only. The biases are the differences from the converged refits, in
units of 1e-4.

HGDP pruned, K = 2, 3, 5 and 8:

- The refits converge in 6-7 iterations, and CV costs 1.76x the passes of
  the main fit.
- The CV errors are 0.4588, 0.4373, 0.4224 and 0.4175 for K = 2, 3, 5 and 8.

| policy | largest bias | passes saved |
|---|---|---|
| 3 iterations | 0.54 | |
| 5 iterations | 0.04 | 24% |
| rtol 1e-6 | 0.013 | 20% |
| rtol 1e-5 | -0.45 (K = 8) | 52% |

Every policy chooses K = 8, as the curve keeps falling.

EUR_subset, K = 1-5:

- The refits take 7-34 iterations, and CV costs 3.2x the passes of the main
  fit.

| policy | bias at K = 1, 2, 3, 4, 5 | passes saved |
|---|---|---|
| 5 iterations | +0.18, -0.05, -0.38, -1.8, -4.3 | 61% |
| rtol 1e-6 | -0.03, -0.04, -0.24, -0.87, -2.65 | 50% |
| rtol 1e-7 | 0, -0.04, -0.02, -0.37, -0.66 | 30% |

- The bias grows with K. The gaps between the errors of successive K are 5
  to 32.
- Every policy chooses K = 2, except no refit at all, which would choose
  K = 5.
- The differences from ADMIXTURE's errors (K = 2-4) are -0.39 to +1.15
  for converged refits, -0.44 to -0.76 for 5 iterations, and -0.61 to
  +0.29 for rtol 1e-6: all within the noise from different folds.

Wall clock, HGDP quick, K = 3, 1 core (the smoke run): converged refits
22.9 s, 5 iterations 13.7 s, rtol 1e-6 19.6 s.

### Part B: fewer full-data iterations

Passes, with the iterations of the final full fit in brackets.

HGDP pruned, K = 5:

| start | seeds | passes (iterations) |
|---|---|---|
| random | 1, 2, 3 | 111 (21), 141 (27), 111 (21) |
| subset 5% | 1, 2, 3 | 58.7 (9), 67.3 (11), 55.0 (9) |
| subset 2% | 2 | 61.1 (11) |
| subset 10% | 1 | 73.8 |
| subset 25% | 1 | 99.8 |
| admixer's mini-batch warm-up (`mbtol3`) | 1 | 86 |

- All these runs reached the same mode.
- Seed 1 took 258 s from a random start and 150 s from a 5% subset.
- admixer's mini-batch warm-up needs a halving tolerance: with admixer's own
  rule the batches never halved, and it ran 100 epochs (306 passes on HGDP
  quick).

Unpruned HGDP, K = 5, seed 1:

| start | passes (iterations) | time |
|---|---|---|
| random | 130 (24) | 1227 s |
| subset 2% | 53.8 (9) | 658 s |
| subset 5% | 69.6 (11) | |

The runs reached the same mode.

Scaling with the number of loci: the gain does not grow. The subset start
needs 0.43-0.53 of the passes on the pruned HGDP and 0.41-0.54 on the
unpruned one, because the full fit after the warm-up takes 9-11 iterations
either way, and that full fit dominates the cost.

Modes, the problem with the subset start:

- HGDP, K = 10, seed 1: a random start took 206 passes (37 iterations). A
  5% subset took 64.3 passes (11 iterations), but **ended 67,600
  log-likelihood units below** the random run.
- EUR, K = 3, seeds 1-5: the best mode was reached by 4 of 5 random runs,
  1 of 5 with 10% subsets, 2 of 5 with 25% and 3 of 5 with 50%. The other
  runs ended 61-64 units below. The median passes were 108 (random), 89
  (10%), 88.5 (25%) and 99.5 (50%).
- EUR, K = 4, seeds 1-5: all the runs reached the same mode. The median
  passes were 129 (random) and 78 (10%).

The SVD start was only run on HGDP quick: 44 passes, against 86-99 from
random starts.

## Decisions for the user

1. **The CV default:** keep refitting folds to convergence, or switch to
   rtol 1e-7, rtol 1e-6 or 5 iterations. On the HGDP every one of these is
   harmless. On EUR the bias grows with K, and with 5 iterations it reaches
   10-15% of the gap between successive K at K = 5. rtol 1e-7 is the
   cautious middle (on EUR: 30% of the passes saved, bias at most 0.66e-4).
2. **The subset start:** keep commit `3cf356f` (an option with warnings, like
   `init = "svd"`) or drop it. It needs about half the passes, but it moves
   runs to lower modes (HGDP K = 10, EUR K = 3).
3. **The default of `init_subset`,** if the option is kept: 2% was as fast as
   5% on the HGDP, but smaller subsets moved more EUR runs to lower modes.

## Open problems

- **The coverage of Part A.** Only seed 1 was run: HGDP K = 2, 3, 5 and 8
  (not 4, 6 or 7), and EUR K = 1-5 (K = 6 and seeds 2-3 were stopped to free
  a core).
- **Statistics of the modes.** HGDP K = 10 had a single seed. Random starts
  also reach different modes at K = 10 (the desktop record of November has
  2 of 5 runs ending 68,000-118,000 units below the best), so this needs
  seeds 1-5 or more.
- **The wall clock of subset passes** was not established: the time ratios
  (0.47-0.58) were a little higher than the pass ratios, under uneven load.
- **The bias of converged refits.** Converged refits also start warm, from
  estimates that have seen the masked genotypes. A reference with cold
  starts was not measured.

## Desktop tests

From the root of the repository on the desktop:

```bash
LIB=<lib> ALGO_SCRATCH=<dir> HGDP_DIR=<hgdp dir> [EUR_BED=...] NC=8 \
  bash data-raw/radmix/newspeed_algo/run_desktop.sh <step>
```

The steps, in order of priority:

1. `warmup hgdp`: K = 5 and 10, seeds 1-5; the methods random, svd, sub2,
   sub5, sub10 and mbtol3. Wall clock, and the question of the modes at
   K = 10.
2. `warmup hgdp_full`: K = 5 and 10, seeds 1-5.
3. `warmup eur`: K = 3 and 4, seeds 1-10, subsets of 10, 25 and 50%.
4. `cv_timing`: the package's `radmix_cv()` for each policy, EUR K = 2-6 and
   HGDP K = 2-8.
5. `cv`: HGDP K = 2-10 with seeds 1-3, and EUR K = 1-6 with seeds 1-5.
6. `summary`: rebuilds the tables.

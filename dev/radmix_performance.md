# Performance record of `gt_radmix()`

What made `gt_radmix()` faster, what was tried and rejected, and what is
left to try, with the measurements behind each decision. The design of the
code is described in `dev/admixture_implementation.md`; this file is its
history. The full reports, scripts and results of the work summarised here
are archived under the git tag `archive/radmix-perf-verify`; "Full
records", at the end, explains how to read them.

Machines: a laptop (Intel i7-10510U, 4 cores, 15 W, AVX2; it throttles by
about 60% after 10-20 s with 4 busy cores) and a desktop (Intel i9-13900,
8 P-cores used, AVX2, no AVX-512). Data: "HGDP quick" (209 individuals x
25,000 loci, the quick benchmark of the article), EUR_subset (379 x
54,051), HGDP (1043 x 130,943) and the unpruned HGDP (1043 x 642,794).
Times are as measured then; the code has only become faster since.

## The changes, in order

Before October 2026, while the code was being ported and restructured
(see "Origins of the code" in `dev/admixture_implementation.md`): the
Q side was split into segments of loci when there are few blocks of
individuals, so that small samples use all the threads; and the fitting
loop moved to C++ with reused buffers, bringing memory down from about 34
copies of the parameters to 10.

| date | change | effect | results |
|---|---|---|---|
| 2026-10-31 (#31) | cross-validation views: fold numbers by a `switch` over 2-10 folds instead of a 64-bit `%`, and a branch-free hide test | passes with a fold hidden 19-25% faster (1.9x to 1.5x the cost of a plain pass); 5-fold CV on HGDP quick 60.5 to 41.2 s | bit-identical |
| 2026-11-02 (#32) | SVD start, `init = "svd"` (an option) | see below: 1.5-2.4x faster first runs at K = 5-10 on the HGDP | different starts |
| 2026-11-02 (#32) | the ALS of the SVD start runs with one BLAS thread | HGDP K = 10, 8 cores: start 5.1 to 2.2 s, and 117 to 10 s of CPU time | rounding |
| 2026-11-02 (#32) | four kernel changes (next table) | fits take 0.57-0.63 of the time on 1 core, 0.52-0.56 on the HGDP at any core count | rounding |
| 2026-11-02 (#32) | default `rtol` 1e-8 instead of 1e-7 | +7-17% time, closer to the optimum (below) | more iterations |
| 2026-11-03 | consolidation of the kernels; plain Q-side derivatives | 2-5% faster (desktop, 1 and 4 cores); 9 copies of the parameters instead of 10 | rounding |

### The kernel changes of #32

Fit times as a fraction of those before (EUR K = 3 and 4, HGDP quick K = 5
and 7, 1 core, medians of 3), each change alone:

| change | laptop | desktop | results |
|---|---|---|---|
| two reciprocals instead of four divisions per genotype (Newton coefficients) | 0.94-0.98 | 0.86-0.87 | rounding |
| decoding fast path: bytes of contiguous individuals read directly, decoded with two 256-entry tables (also with a fold hidden) | 0.86-0.95 | 0.89-0.91 | bit-identical |
| two loci per pass in the Q- and P-side accumulations (half the loads and stores of the accumulators) | 0.86-0.94 | 0.83-0.95 | rounding |
| `target_clones("arch=x86-64-v3", "default")` (AVX2 + FMA) on the Newton, EM and gradient kernels, x86-64 Linux with GCC >= 12 | 0.71-0.86 | 0.70-0.83 | rounding |
| all four | 0.61-0.63 | 0.57-0.60 | |
| all four, built with `-DRADMIX_NO_MULTIVERSION` | 0.70-0.76 | 0.73-0.75 | |

On 4 and 8 cores the combination takes 0.56-0.63 of the time, and 0.62 with
5-fold CV. Iterations and CV errors are unchanged, log-likelihoods agree to
about 1e-8 and Q to about 1e-11. The AVX2 and baseline versions differ by
the same amount (Q 2e-12, log-likelihood 4e-9), so results can differ by
rounding between machines with and without AVX2; on one machine they are
identical for any number of threads. The laptop and the desktop gave the
same iterations and log-likelihoods (to the 15 digits stored) for every
fit.

### The SVD start

Why: ADAMIXTURE needed 6 iterations where `gt_radmix()` needed 16-21, at
the same cost per iteration (HGDP quick, laptop, 1 core: `gt_radmix()` 8.2
s and 16 iterations, ADAMIXTURE 4.2 s and 6, ADMIXTURE 31.9 s and 22,
OpenADMIXTURE.jl 24.7 s and 22, fastmixture 68.8 s and 315 EM-type
iterations, all at the same log-likelihood). The difference was its start.

HGDP, desktop, 8 cores, rtol 1e-7, seeds 1-5, medians. "First run"
includes the SVD; further runs of the same call share it:

| K | random: time (iterations) | SVD: SVD + start + fit (iterations) | first run | further runs |
|---|---|---|---|---|
| 3 | 3.6 s (13) | 1.6 + 0.6 + 1.5 s (6) | 0.97x | 1.7x |
| 5 | 7.4 s (20) | 1.6 + 1.1 + 2.2 s (6) | 1.5x | 2.3x |
| 7 | 21.3 s (37) | 1.6 + 2.1 + 5.5 s (11) | 2.3x | 2.8x |
| 10 | 26.4 s (32) | 1.6 + 3.2 + 6.1 s (8) | 2.4x | 2.8x |

All SVD runs ended at the best mode (at K = 7 only with rtol 1e-8: at 1e-7
they stopped 1.5-3.7 units short); 2 of 5 random runs at K = 10 ended in
modes 68,000 and 118,000 units below it. But SVD runs explore few modes
(the ALS solutions are nearly deterministic), and on the HGDP quick subset
at K = 7 they all ended in a mode 7,800 units below the best (as did
ADAMIXTURE), so random starts stay the default and the recommended choice
for analyses of modes with `gt_clumppling()`. On tiny data (the Anolis) the
SVD costs more than a fit. The choices (partial SVD, K + 5 components,
refined ALS, best of 3) come from `data-raw/radmix/explore_svd_start_refined.R`:

- `gt_pca_partialSVD()` gives the same subspace as `gt_pca_randomSVD()`
  about 5x faster (HGDP, rank 7: 8.5 against 45.5 s);
- fastmixture's ALS (unweighted, clipped) converged slowly to poor modes on
  EUR K = 3 (36 iterations, 5-44 units below the best); ADAMIXTURE's
  refinements (bounded least squares, Fisher-weighted Q step, momentum 0.8)
  fixed that;
- K components missed weak structure (EUR K = 3: 2 of 5 runs more than 5
  units below the best, up to 61); K + 5 put all runs within 1 unit of it,
  in 7 iterations instead of 12; 2K was no better and once led to a poor
  mode;
- the best of 3 ALS solutions (by log-likelihood) removed the poor modes
  of single solutions (HGDP quick K = 5, Africa split in two, 19,800 units
  below); the best of 5 added nothing.

### The convergence threshold

Desktop, 1 core, seed 1, before the kernel changes. "Short" is the
log-likelihood short of a fit run to ADMIXTURE's threshold (an absolute
1e-4), "Q diff" the largest difference of Q from that fit:

| data, K | rtol 1e-7: time, iterations, short, Q diff | rtol 1e-8 | absolute 1e-4 |
|---|---|---|---|
| EUR 3 | 10.2 s, 19, 2.70, 0.050 | +17%, 23, 0.08, 0.009 | +43%, 29 |
| EUR 4 | 14.1 s, 22, 0.64, 0.026 | +7%, 24, 0.17, 0.015 | +32%, 31 |
| HGDP quick 5 | 3.3 s, 16, 0.07, 1e-4 | +10%, 18, 1e-5, 1e-5 | +13%, 19 |

With the kernel changes, 1e-8 takes 0.59-0.70 of the time 1e-7 took before
them. ADAMIXTURE stops at an absolute 0.1; ours is still the loosest rule
of the programs compared.

## Where the time goes

`perf` profile of a fit before the kernel changes (EUR K = 4, 1 core):
Hessian accumulation 23%, the Newton coefficients (four divisions per
genotype) 21%, genotype decoding 17%, gradients 11%, `qp` 9%, the
log-likelihood 7%; EM, the QP solver, the quasi-Newton step and R together
under 5%. `perf` charges stalls to the instruction that waits, so a share
is not a gain: the lines with the divisions held 21% of the samples, but
two reciprocals instead of four divisions saved 2-6% of a fit on the
laptop (13-14% on the desktop).

Scaling, desktop, HGDP, whole fits (seconds, speed-up over 1 core; threads
bound to the P-cores):

| case | 1 core | 2 | 4 | 6 | 8 |
|---|---|---|---|---|---|
| K = 5 (20 iterations) | 53.5 | 27.1 (2.0) | 13.7 (3.9) | 9.6 (5.6) | 7.6 (7.1) |
| K = 10 (35 iterations) | 197.0 | 99.5 (2.0) | 50.8 (3.9) | 36.8 (5.4) | 28.4 (6.9) |
| unpruned HGDP, K = 5 (23 iterations) | 309.1 | | 78.9 (3.9) | | 43.2 (7.2) |

The kernels alone reach 7.4-8.1x on 8 cores, with the same clock and
slightly higher instructions per cycle than on 1 core, and a pass reads
the genotypes at 2-4 GB/s, far below the memory bandwidth. What is lost
comes from the granularity of the Q side (17 blocks of 64 individuals x 2
segments = 34 tasks, at most ~95% efficiency on 8 threads; the P side has
~2,000 tasks and reaches 98-101%) and from 1-2% of serial work per fit (QP
solves, quasi-Newton step, R). Beyond 16 threads the 34 Q-side tasks
(`MIN_Q_TASKS` = 32 in `src/radmix_kernels.cpp`) would limit the Q side.
Results are bit-identical for every number of threads and binding. Without
binding, OpenMP threads drift onto the slower E-cores of this CPU.

## Tried and rejected

Listed so that they are not retried without new evidence.

- One reciprocal per genotype instead of two: no faster (0.98-1.04), and
  `qp = 0` would give NaN.
- ADAMIXTURE's expansion of the Q-side Hessian, (1 - p)(1 - p)' = 11' -
  1p' - p1' + pp': slower on the laptop (1.01-1.06; the loop is bound by the
  loads and stores of the accumulators, not by arithmetic), faster on the
  desktop (0.87-0.97) but partly from code layout (the same instructions on
  the P side in 10% fewer cycles). Moot since 2026-11-03: the Q side now
  uses the plain derivatives, whose Hessian is a single sum of d2 p p'.
- Longer P-side tiles (128 or 256 loci): 1.02-1.09 (the TILE_I x K corner
  of Q no longer stays in L1).
- `-fno-trapping-math`: no effect (the loops have no branches).
- AVX2 clones of the log-likelihood and the deviance: 7-9% slower on the
  laptop (scalar `log()` calls inside the loop), so they are not cloned.
- BLAS products (GEMMs) for the Q-side sums: faster only for K >= ~8 and
  only with an optimised BLAS; 3-7x slower with the reference BLAS (the
  default of the CRAN builds of R for Windows and macOS), and R's BLAS
  threads compete with OpenMP's.
- Avoiding the per-genotype `log()` (as admixer does): the logarithms are
  already batched, 8 per 64 genotypes.
- Reading the genotypes through bigstatsr's `SubBMCode256Acc`: 19% slower
  than the custom view.
- For the SVD start: perturbing it towards a random start (iterations back
  to those of random starts, no better exploration of modes), and mixing
  SVD and random starts in one analysis (it would bias the modes found by
  `gt_clumppling()` towards those the ALS favours).

## Leads not pursued yet

- `-march=native` is 11-16% faster than the multiversioned build (desktop),
  mostly in the log-likelihood (0.59-0.64 of the baseline with
  `-march=native`, 0.86-0.89 multiversioned) and in code outside the cloned
  kernels. Compiling the whole file for AVX2 makes the log-likelihood faster
  too, so the lead is how it is cloned (or a vectorisable `log()`).
- `-O3` saves 2-5% of a fit with identical results, from inlining or loop
  unswitching rather than vectorisation; which loops is not known, and a
  package should not override R's optimisation flags.
- Cross-validation refits each fold to convergence (7-11 iterations on
  HGDP quick); ADAMIXTURE polishes each fold with 5 warm-started
  iterations. The biggest lever for CV, but it changes the CV errors, and
  the warm start has seen the held-out genotypes: it needs a statistical
  check (same K chosen, agreement with ADMIXTURE).
- A decoded in-memory copy of the genotypes (1 byte or 2 bits each), with
  the fold mask applied once: passes with a fold hidden would cost the same
  as plain ones (now 1.5x).
- The log-likelihood of the quasi-Newton point folded into the next P
  update, which reads the same genotypes: 4 passes instead of 5 per
  accepted iteration (admixer measured ~4%).
- admixer's damped retries of a rejected quasi-Newton step (1.2-4.2x on
  hard data there, nothing elsewhere; ours accepted 16 of 16 steps on HGDP
  quick) and its mini-batch warm-up between EM and the main loop (there,
  K = 5 28.8 to ~20 s, K = 20 126.7 to 54.7 s).
- More Q-side tasks for many cores (see "Where the time goes"), and tile
  sizes and the OpenMP schedule, never tuned.
- For K >= 8, single-precision Hessians or BLAS products (admixer: 6-35%).

## How to measure

- Time installed builds only (`R CMD INSTALL --preclean`, and check that
  the log shows `-O2`): `devtools::load_all()` compiles at `-O0`, and its
  leftover `.o` files once made `gt_radmix()` 6x slower.
- One core first; medians of at least 3 runs; alternate the builds (A B A
  B) with pauses; pin the threads (`taskset`, or `OMP_PROC_BIND=close
  OMP_PLACES=cores`); check with repeats that the noise is below the
  effect. Judge multi-core scaling on a desktop, not on a laptop that
  throttles.
- Effects of up to ~10% on one kernel can be code-layout artefacts on one
  CPU: trust changes that help on two different CPUs.
- Correctness gates: `test_radmix_vs_julia.R` (identical iteration
  counts), `test_radmix_vs_admixture.R`, the kernels against the R
  reference (`helper-radmix-reference.R`), and the `expect_identical()`
  checks of thread invariance; any sum whose order depends on the number
  of threads breaks them.
- Tools in `data-raw/radmix/`: `profile_eur.R` (per-kernel times),
  `benchmark_cv_access.R` (passes with and without a fold hidden),
  `benchmark_admixture.R`; and the benchmark of the article
  (`vignettes/articles/admixture_benchmark.Rmd`), which compares the
  current code with ADMIXTURE, fastmixture, ADAMIXTURE and
  OpenADMIXTURE.jl.

## Full records

The full reports of this work, with the scripts, patches and results
behind them, are not in the working tree: they are archived under the git
tag `archive/radmix-perf-verify` (commit f423896), pushed to GitHub. The
tag keeps them reachable after the branches of that work
(`radmix-perf-verify`, `kernel-speedups`, `svd-start`, `radmix-speedups`)
are deleted; it also keeps the commit the patches apply to (7baa46c, the
SVD start before the kernel changes). Everything is under
`data-raw/radmix/perf_verification/` in the tagged commit:

| path | content |
|---|---|
| `notes/speedup_plan.md` | the plan: the comparison with ADAMIXTURE and admixer, the `perf` profile, the experiments ranked, the exploration of the SVD start, the BLAS measurements |
| `notes/kernel_experiments_laptop.md` | the log of the four kernel experiments on the laptop |
| `HANDOVER.md` | the summary of the laptop work, the decisions left open, and the protocol of the verification |
| `desktop_results/REPORT.md` | the report of the verification on the desktop, laptop and desktop figures side by side, and the evidence for the decisions |
| `desktop_results/E_scaling/scaling.png` | the plot of the multi-core scaling |
| `patches/` | one patch per experiment (on 7baa46c), including the rejected `hessian_expansion_REJECTED.patch`, and the BLAS-thread fix |
| `scripts/` | the builds of the variants, the benchmarks, the summaries |
| `laptop_results/`, `desktop_results/*/` | the CSV results and summaries of each phase (the `.rds` outputs were never committed) |

To read them in a clone of the repository:

```bash
# fetch the tag, if `git tag -l 'archive/*'` does not list it (a fresh
# clone has it; an older clone does not fetch it once the branch is gone)
git fetch origin tag archive/radmix-perf-verify
# list the files, or read one
git ls-tree -r --name-only archive/radmix-perf-verify -- data-raw/radmix/perf_verification
git show archive/radmix-perf-verify:data-raw/radmix/perf_verification/desktop_results/REPORT.md
# or check out the whole archive next to the repository, to run the scripts
git worktree add ../radmix-perf-archive archive/radmix-perf-verify
```

On GitHub:
<https://github.com/EvolEcolGroup/tidygenclust/tree/archive/radmix-perf-verify/data-raw/radmix/perf_verification>.

Copies of the four reports were also kept in `data-raw/radmix/speedups/`
on `dev` (commit 3c3fe1b, #32) until 2026-11-03, when this file replaced
them.

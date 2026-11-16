# Performance record of `gt_radmix()`

What made `gt_radmix()` faster, what was tried and rejected, and what is
left to try, with the measurements behind each decision. The design of the
code is described in `dev_notes/admixture_implementation.md`; this file is its
history. The full reports, scripts and results of the work summarised here
are archived under the git tags `archive/radmix-perf-verify` (October and
early November 2026) and `archive/newspeed` and `archive/supervised` (the
work of 13-15 November 2026); "Full records", at the end, explains how to
read them.

Machines: a laptop (Intel i7-10510U, 4 cores, 15 W, AVX2; it throttles by
about 60% after 10-20 s with 4 busy cores) and a desktop (Intel i9-13900,
8 P-cores used, AVX2, no AVX-512). Data: "HGDP quick" (209 individuals x
25,000 loci, the quick benchmark of the article), EUR_subset (379 x
54,051), HGDP (1043 x 130,943) and the unpruned HGDP (1043 x 642,794).
Times are as measured then; the code has only become faster since.

## The changes, in order

Before October 2026, while the code was being ported and restructured
(see "Origins of the code" in `dev_notes/admixture_implementation.md`): the
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
| 2026-11-13 | fused pass: the log-likelihood of a point computed with the next update of Q (`update_q_loglik()`) | 4 passes per accepted iteration instead of 5 (a fit of 8 iterations: 33 passes instead of 41, the start included); fits 0.93-0.95 of the time alone | log-likelihoods by rounding |
| 2026-11-13 | cross-validation: the hidden genotypes of a fold as a bit mask, built once per fold fit | passes with a fold hidden cost 1.08-1.12x a plain one (were 1.7-1.8x on the desktop); mask 17 MB for the HGDP, 83 MB unpruned | bit-identical |
| 2026-11-13 | `slice_geno()` and `slice_coef()` forced inline into the AVX2 clones (they were calls to the baseline code) | with the next change, below | rounding (FMA) |
| 2026-11-13 | cheaper log-likelihood of a slice (scalar fix-up loop skipped for hard calls, arithmetic masks), inlined into the AVX2 clone in the fused pass | all four changes: fits take 0.75-0.85 of the time (desktop, 1-8 threads, EUR K = 3-4, HGDP K = 5-10, unpruned HGDP K = 5) | rounding, at most 3e-15 relative in the log-likelihood, same iterations |
| 2026-11-15 | `MIN_Q_TASKS` 128 instead of 32 (more Q-side tasks for samples under 8,192 individuals) | fits 3-6% faster on 8 threads, 0-3% on 1 (desktop) | rounding |
| 2026-11-15 | refits of cross-validation to `cv_rtol = 1e-7` instead of `rtol` | CV passes 0.46 of before on EUR_subset, 0.89 on the HGDP (see "Cross-validation refits") | CV errors biased by at most 2e-4, same K chosen |

The changes of 13 and 15 November were made on the branches `newspeed` and
`supervised` (from b5999cc), measured on the desktop on 15 November, and
ported onto `dev` after #35 (genotype likelihoods), which made the kind of
data a template parameter of the kernels; the timings below are those of
the branches, before the port.

### The changes of 13 November 2026 (branch `newspeed`)

Laptop first (1 core, the other three busy with other jobs, load 4-11, so
indicative only), HGDP 1043 x 130,943, K = 5, single passes from the same
estimates (medians of 2 runs of medians of 3, builds alternated), seconds:

| build | loglik | update_q | update_p | update_q_loglik | accepted iteration |
|---|---|---|---|---|---|
| before (b5999cc) | 2.13 | 1.84 | 2.14 | | 10.09 (5 passes) |
| fused pass | 2.20 | 2.08 | 2.30 | 3.31 | 9.99 (0.99) |
| + inlining, cheaper log-likelihood, called in baseline code | 1.83 | 1.65 | 1.89 | 2.77 | 8.20 (0.81) |
| + log-likelihood inlined into the AVX2 clone (the default) | 1.75 | 1.68 | 1.94 | 2.27 | 7.83 (0.78) |

Desktop (i9-13900, idle machine, threads pinned to the P-cores, installed
`-O2` builds; unattended run of 15 November 2026), builds `base` (b5999cc),
`fused` (the fused pass alone), `inline` (fused + inlining), `final` (all,
13f89be) and `final_llnoinl` (the log-likelihood of the fused pass called
in the baseline code). Single passes, 1 thread, seconds; an iteration is
2 `update_q` + 2 `update_p` + `loglik` before, `update_q` +
`update_q_loglik` + 2 `update_p` with the fused pass:

| data, K | base: iteration | fused | inline | final | final_llnoinl |
|---|---|---|---|---|---|
| EUR, 4 | 0.304 | 0.277 (0.91) | 0.261 (0.86) | 0.222 (0.73) | 0.241 (0.79) |
| HGDP, 5 | 2.355 | 2.170 (0.92) | 1.938 (0.82) | 1.772 (0.75) | 1.847 (0.78) |
| HGDP, 10 | 4.942 | 4.591 (0.93) | 4.091 (0.83) | 3.980 (0.81) | 3.971 (0.80) |

Inlining the log-likelihood into the AVX2 clone makes the fused pass
10-22% cheaper at K <= 5 (`update_q_loglik` 0.521 against 0.646 s on the
HGDP) and changes nothing at K = 10, so it is the default
(`RADMIX_FUSED_LL_NOINLINE` restores the call). Fits of 8 iterations
(medians of 3, 1, 4 and 8 threads) take 0.93-0.95 of the time with the
fused pass alone and 0.76-0.80 with all the changes, on the HGDP (K = 5
and 10) and the unpruned HGDP (K = 5); whole fits to convergence (EUR
K = 3-4, HGDP K = 5-10, 1-8 threads) 0.75-0.85, with the same iterations
and log-likelihoods within 3e-15 (relative). Every build gives identical
results on 1, 4 and 8 threads. Fold fits of cross-validation (HGDP and
unpruned HGDP, K = 5, 4 iterations): with the hash, 1.7-2.0x the time of
the same fit without a fold hidden; with the bit mask, 1.08-1.12x.

Decisions (Andrea, 15 November 2026): keep all four changes, accept the
rounding differences that FMA brings into the clones (at most 3.5e-15 in
the log-likelihood, the same iterations), inline the log-likelihood into
the clone, build the bit mask automatically (when it takes less than
1 GiB), and raise `MIN_Q_TASKS` to 128 (next section).

Lessons: GCC did not inline `slice_geno()` (qp and the decoding of every
slice) into the AVX2 clones, so a share of every pass ran in SSE2 code
(check with `objdump -d` that the clones only call `memset()` and `log()`);
and in the baseline code, `(double) (x > y)`, `std::min()`, `std::max()` and
`omp simd reduction(max:)` keep a loop from vectorising (a first version of
the cheaper log-likelihood was 1.5x slower for that), while sums of products
do not. The fused pass alone made `update_q` and `update_p`, whose code it
does not change, 5-10% slower on the laptop, a layout effect that
disappeared with the inlining.

#### The port onto `dev` (PR #36), against `dev`

Desktop (i9-13900, idle machine after a 15-min idle check, threads pinned:
1 thread on CPU 2, 8 on the P-cores 0-7; installed `-O2` builds with the
AVX2 clones), 16 November 2026: `dev` (a20eb08) against `newspeed_port`
(2e2f5b4), whole fits to convergence (`rtol` 1e-8, seed 1), 3 rounds
alternated, with the kit of the tag `archive/newspeed`
(`run_bench.sh`, `CASES=hgdp:5,hgdp:10,eur:3,eur:4 THREADS="1 8"`).
Medians in seconds (min-max of the repeats within 2%), ratio port / dev:

| data, K | 1 thread: dev | port | ratio | 8 threads: dev | port | ratio |
|---|---|---|---|---|---|---|
| EUR, 3 | 6.26 | 5.02 | 0.80 | 0.89 | 0.69 | 0.78 |
| EUR, 4 | 7.52 | 6.22 | 0.83 | 1.06 | 0.84 | 0.80 |
| HGDP, 5 | 50.31 | 40.02 | 0.80 | 6.68 | 5.28 | 0.79 |
| HGDP, 10 | 183.32 | 168.01 | 0.92 | 25.79 | 21.89 | 0.85 |

The same iterations, log-likelihoods within 7.4e-16 (relative), and
identical results on 1 and 8 threads for both builds. The gain is that of
`newspeed` (0.75-0.85) except at HGDP K = 10 on 1 thread: there `dev`
became faster with #35 (183 s against 198 s for b5999cc on 15 Nov) and the
port is a little slower than `newspeed` (168 against 161 s), so the gain
shrinks from 0.81 to 0.92. Comparisons across the two days are indicative
only (separate runs). Summary, fits and the
driver script: `dev_notes/radmix_port_timing_2026-11/`.

### Tiles and the OpenMP schedule

Desktop, fits of 8 iterations (medians of 3), against `final`: with
`MIN_Q_TASKS` = 128, 0.94-0.97 of the time on 8 threads and 0.97-1.00 on 1
(HGDP K = 5 and 10, unpruned HGDP K = 5); 64 gained less (0.96-0.98 on 8
threads); `TILE_J` = 32 changed nothing (0.99-1.02); a guided schedule was
slower (1.05-1.07 on 8 threads). With 128, the HGDP has 8 segments of loci
(136 Q-side tasks) instead of 2 (34). The schedule never changes the
results; `MIN_Q_TASKS` changes them by rounding (4e-16 to 3e-15).

### Supervised analyses (branch `supervised`)

A supervised fit skips the Q-side work of the reference individuals (their
sums and quadratic programs), so it is faster than an unsupervised one
with the same data, and it agrees with ADMIXTURE `--supervised` (HGDP
K = 5, 660 or 139 references: Q of the individuals to estimate within
6e-5, log-likelihoods 1-4 units higher). On the desktop, at 8 threads, it
ran about 50 times as fast as ADMIXTURE `--supervised -j8` (HGDP K = 5:
3.3 s against 158 s with 660 references, 4.9 s against 264 s with 139;
an unsupervised fit of the same data takes 7.1 s). Its fits do not use the
fused pass (their Q side does not read the references, so it cannot give
the log-likelihood).

### Cross-validation refits

Each fold is refitted from the estimates for the full data; ADAMIXTURE
polishes each fold with 5 warm-started iterations. A refit stopped early
predicts the held-out genotypes too well (the start has seen them), so the
error is biased downwards. Desktop (E-cores, the measures are counts, not
times; `newspeed_algo`, 13 November 2026), the trajectory of each refit
scored after every iteration, HGDP K = 2-10 with seeds 1-3 and EUR_subset
K = 1-6 with seeds 1-5; bias = the error minus that of refits converged to
`rtol` = 1e-8, x 1e-4; passes = the CV passes relative to converged
refits:

| policy | HGDP bias | HGDP passes | EUR bias, K = 1-4 | EUR bias, K = 5 / 6 | EUR passes | K chosen |
|---|---|---|---|---|---|---|
| `cv_rtol` 1e-7 | <= 0.007 | 0.89 | <= 0.15 | -1.2 / -2.0 | 0.46 | unchanged |
| `cv_rtol` 1e-6 | <= 0.009 | 0.83 | <= 0.62 | -3.4 / -6.7 | 0.26 | unchanged |
| at most 5 iterations | <= 0.06 | 0.74 | <= 1.4 | -5.7 / -10.3 | 0.20 | unchanged |
| `cv_rtol` 1e-5 | <= 0.41 | 0.55 | <= 1.4 | -5.7 / -10.2 | 0.19 | unchanged |

On the HGDP, converged refits take about 7 iterations, so every policy is
harmless but saves little; on EUR_subset they take 24 on average (CV costs
2.9 times the main fits), which is where the savings are, and where the
bias grows with K (successive K differ by 5e-4 to 3e-3). With no refit at
all, EUR would choose K = 6 instead of 2. Decision (Andrea, 15 November
2026): `cv_rtol = 1e-7` by default (or `rtol`, if larger).

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
(Since 15 November 2026, `MIN_Q_TASKS` = 128: 136 Q-side tasks on the
HGDP, see "Tiles and the OpenMP schedule".)
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
- A subset start (`init = "subset"`, 13 November 2026, never merged): a
  random start fitted to a random subset of the loci (rtol 1e-6), P
  extended to all the loci given its Q, then the fit of all the loci. It
  took 0.4-0.5 of the passes of a random start on the HGDP (pruned or not:
  the full fit needs 9-11 iterations after it, whatever the size of the
  subset), but it moved runs to lower modes (HGDP K = 10: 67,600
  log-likelihood units below; EUR K = 3: the smaller the subset, the more
  often a lower mode), and the gain does not grow with the number of loci.
  Dropped by Andrea on 15 November 2026. admixer's mini-batch warm-up was
  less effective (HGDP K = 5, 86 passes). Records: `archive/newspeed`,
  `data-raw/radmix/newspeed_algo/`.

## Leads not pursued yet

- `-march=native` is 11-16% faster than the multiversioned build (desktop),
  mostly in the log-likelihood (0.59-0.64 of the baseline with
  `-march=native`, 0.86-0.89 multiversioned) and in code outside the cloned
  kernels. Partly addressed on 2026-11-13: `slice_geno()` and
  `slice_coef()` were outside the clones (now inlined into them), and the
  log-likelihood now mostly runs inside the clone of the fused pass. To
  check on the desktop: what is left between the multiversioned build and
  `-march=native`; a vectorisable `log()` would only replace the 8 batched
  logarithms per 64 genotypes.
- `-O3` saves 2-5% of a fit with identical results, from inlining or loop
  unswitching rather than vectorisation; which loops is not known, and a
  package should not override R's optimisation flags.
- (Done 2026-11-15: the refits of cross-validation stop at
  `cv_rtol = 1e-7` by default, see "Cross-validation refits".)
- (Done 2026-11-13, as a bit mask rather than a decoded copy of the
  genotypes, 1/8 of the memory: passes with a fold hidden cost as plain
  ones.) A decoded copy (1 byte per genotype) would also make the rows of
  any subset contiguous (the fast path of the decoding), at the memory of
  the genotypes (0.7 GB for the unpruned HGDP); not tried.
- (Done 2026-11-13, folded into the next Q update rather than P, as that
  pass reads the same (q, p): 4 passes instead of 5 per accepted
  iteration.)
- admixer's damped retries of a rejected quasi-Newton step (1.2-4.2x on
  hard data there, nothing elsewhere; ours accepted 16 of 16 steps on HGDP
  quick) and its mini-batch warm-up between EM and the main loop (there,
  K = 5 28.8 to ~20 s, K = 20 126.7 to 54.7 s).
- (Done 2026-11-15: `MIN_Q_TASKS` = 128, see "Tiles and the OpenMP
  schedule"; `TILE_J` and the schedule were measured and kept. They can
  still be set at compile time, `-DRADMIX_MIN_Q_TASKS`, `-DRADMIX_TILE_J`
  and `-DRADMIX_OMP_SCHEDULE`, for benchmarks.)
- A fused pass for supervised fits (update every column of Q, then restore
  the references), which would give them the 4 passes per iteration of
  unsupervised fits.
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

The work of 13-15 November 2026 is archived in the same way, under two
tags that keep the branches `newspeed` and `supervised` as they were at
the end of the desktop run (the package code of both was ported onto
`dev`, without the subset start):

| tag (commit) | path | content |
|---|---|---|
| `archive/newspeed` (ec96315) | `data-raw/radmix/newspeed/` | the kit of the per-iteration speed-ups (builds of the variants, benchmarks, summaries), `REPORT.md` (laptop), `HANDOVER.md` (the desktop protocol and the decisions), `desktop_results/*/summary.md` and the raw CSV files of the desktop run |
| `archive/newspeed` | `data-raw/radmix/newspeed_algo/` | the CV refit policies and the warm-ups (subset start): scripts, `REPORT.md`, `results_laptop/`, `results_desktop/cv_policies.txt` |
| `archive/supervised` (fc43f9b) | `data-raw/radmix/supervised/` | the validation against ADMIXTURE `--supervised`, the timings, `REPORT.md`, `results_desktop/` |

```bash
git fetch origin tag archive/newspeed tag archive/supervised
git show archive/newspeed:data-raw/radmix/newspeed/HANDOVER.md
```

Copies of the four reports were also kept in `data-raw/radmix/speedups/`
on `dev` (commit 3c3fe1b, #32) until 2026-11-03, when this file replaced
them.

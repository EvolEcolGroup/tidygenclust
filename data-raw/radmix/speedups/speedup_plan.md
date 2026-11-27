> Copied from `data-raw/radmix/perf_verification/notes/speedup_plan.md` on the branch `radmix-perf-verify` (commit f423896),
> which also holds the scripts, patches and results (CSV) it refers to.

> Copied on 2026-11-01 from the laptop notes (Dropbox, radmix_speedup_plan.md): the history of the speed-up work (comparison with ADAMIXTURE and admixer, perf profile, SVD start exploration, BLAS measurements). Scratchpad paths refer to the laptop.

# Speed-ups for gt_radmix(): what to learn from ADAMIXTURE and admixer

Written on 2026-10-31 (branch `adamixture`), after adding ADAMIXTURE to the
benchmark. Not implemented yet: this is a ranked list of experiments, to be
picked up later. The aim is not to change our approach, only to find simple
tweaks worth trying.

## 1. Where we stand

`gt_radmix()` (C++/OpenMP, a modified port of OpenADMIXTURE.jl): random
start, 5 EM priming steps, then block updates of Q and P (SQP: one
Newton/QP step per column), with quasi-Newton (QN) acceleration using 3
secant pairs; stops when the relative change of the log-likelihood is
< 1e-7; cross-validation as in ADMIXTURE (each fold refitted from the full
estimates).

Quick HGDP benchmark (`vignettes/articles/data/admixture_benchmark_quick`,
209 individuals x 25,000 loci, K = 5, seed 1, i7-10510U laptop), all programs
reaching the same log-likelihood (-3872985.2):

| 1 core | no CV | 5-fold CV | iterations |
|---|---|---|---|
| gt_radmix (R session) | 8.2 s | 49.2 s | 16 |
| gt_radmix (command line, incl. R start-up) | 10.8 s | 48.8 s | 16 |
| ADAMIXTURE 1.8.9 (command line, incl. python start-up and a plot) | 4.2 s | 16.0 s | 6 |
| ADMIXTURE 1.3.0 | 31.9 s | 55.7 s | 22 |
| fastmixture 1.3.0 | 68.8 s | 121.4 s | 315 |
| OpenADMIXTURE.jl (Julia session) | 24.7 s | - | 22 |

gt_radmix in an R session on 1 / 2 / 4 cores: 8.2 / 6.5 / 5.0 s. The poor
scaling is mostly thermal on this 15 W laptop (all programs scale as badly;
`data-raw/radmix/Notes_to_do`: "slows down by about 60% after 10-20 s with 4
busy cores"), so measure scaling on a desktop before optimising it.

## 2. What the comparison found

### Per iteration we are already as fast as ADAMIXTURE

~0.5 s per iteration for gt_radmix (8.2 s / 16), ~0.47 s for ADAMIXTURE (its
CV: 11.8 s for 5 folds x 5 iterations; its log prints `Iteration i,
Log-likelihood, Time:`). Its default CPU solver, BR-QN
(`adamixture/model/br_qn.py::optimize_original`), is the same loop as ours:
Q sweep, P sweep, Q sweep, P sweep, ZAL quasi-Newton with 3 secant pairs,
log-likelihood of the QN point (5 passes over the genotypes; 6 when the QN
step is rejected), one log per 8 genotypes. Its wheel is built with
`-O3 -ffast-math` but SSE2 only (no AVX/FMA), like ours (R's `-O2`).

### Main fit: the gap is the number of iterations (16 vs 6)

ADAMIXTURE starts from a randomised SVD + ALS (also fastmixture's default,
and admixer's `--svd`):
- `model/svd.py::_RSVD_once`: rank max(K+10, 20), genotypes centred and
  scaled with Patterson weights 1/sqrt(2 max(f(1-f), 1e-4)), missing set to 0,
  float32, chunks of 8192 loci, at most 5 power iterations (stop when the
  top singular values change by < 0.1); 4-7 cheap passes in all.
- `model/als.py::ALS`: alternating bounded least squares on the SVD factors
  only (never touches G again): reconstruction `[0.5 U S, f] [V, 1]'`, P step
  unweighted, Q step Fisher-weighted (w = 1/max(f(1-f), 1e-3), mean 1),
  momentum 0.8 on Q, stop when the RMSE of Q changes by < 1e-4 (max 1000);
  bounded LS by `bvls.pyx` (unconstrained pinv first, block principal
  pivoting only for the rows that violate the bounds).
- Its stopping rule (|delta loglik| < 0.1 absolute, or 3 iterations without
  improvement) is stricter than ours here (1e-7 x 3.9e6 = 0.39), so it is not
  the reason for its fewer iterations.

### Cross-validation: two causes

1. **Refit policy.** ADAMIXTURE (`src/cv.py`) polishes each fold with exactly
   5 BR-QN iterations warm-started from the full estimates
   (`polish_sqp_qn(n_iters=5, tol=0.0)`), keeping the best; we refit each
   fold to convergence (`R/radmix_cv.R:44-47`, full `max_iter` and `rtol`):
   on HGDP quick our folds took 8, 11, 7, 11 and 8 iterations. Its error is
   the same as ours (mean deviance of the held-out genotypes). Caveat: the
   warm start contains the held-out genotypes, so a short polish may bias
   the CV error downwards (untested).
2. **Cost of a CV pass.** Our CV passes cost ~2x a normal pass (~175 ms vs
   ~85 ms per pass). Every genotype read in CV mode
   (`src/radmix_geno_access.h:207-235`) runs a splitmix64 hash, a 64-bit
   `z % _n_folds` by a run-time value and a data-dependent branch that hides
   ~20% of the genotypes (hard to predict). ADAMIXTURE masks the fold once,
   in memory (set to missing), so its CV passes cost like normal ones.

Passes on HGDP quick: main fit 10 (EM) + 1 + 16 x 5 = 91; CV 5 x (1 + ~9 x 5)
+ 5 deviance passes = ~235.

### admixer (github.com/aalbrechtsen/admixer, v0.2.13, MIT, C++17)

Same scheme as ours (random start, 5 EM steps, P/Q Newton/QP block updates,
3-pair QN; absolute tolerance 1e-4; no CV). Its speed-ups, with its own
measurements (`bench/perf_results.md`):
- **Mini-batch warm-up** between EM and QN (`src/fit.hpp`,
  `Fitter::warmup`): loci shuffled once; B = min(32, M/500) batches; each
  batch does a P step on its rows and a Q step from its loci only; B halves
  whenever an epoch does not raise the full log-likelihood, until 1 (max 100
  epochs). K = 5: 28.8 s / 34 QN iterations -> ~20 s / 15; K = 20: 126.7 s /
  54 -> 54.7 s / 10. Off for genotype likelihoods (found the best optimum
  less often).
- **Branch-free, vectorised per-genotype kernel**
  (`Genotypes::row_pass<MODE>`, `src/io.hpp:116-164`): genotypes widened to
  double first, selects compiled to blends, one division per genotype
  (`inv = 1/(h(1-h)); ih = (1-h) inv; ia = h inv`); vectorised by GCC only
  with `-fno-trapping-math` (1.6-2.1x on this step alone on AVX2); 12.65 ->
  3.02 ns per genotype; 1.6-2.7x end to end at K = 5.
- **Log-likelihood without per-entry log()**: 8 running products, exponent
  bits moved to a counter every 16 factors, one log per lane per row.
- **Log-likelihood folded into the P step**: `sqp_P` returns log L of its
  input, and the QN point's P step becomes the next iteration's first half,
  so an accepted iteration costs 4 passes (K = 5: 22.6 -> 21.8 s).
- **Damped retries** of a rejected QN step (scaled 0.5, then 0.25, before
  falling back): 1.2-4.2x on hard data with many rejections, nothing
  elsewhere.
- Hessians as BLAS products on tiles (single-threaded BLAS inside OpenMP),
  float Hessians in the QN phase for K >= 8 (6-35%), register-blocked
  kernels for K = 2-7 (9-27%), tile shapes (P: 128 loci x 1024 individuals,
  dynamic; Q: 1024 individuals x 64 loci, up to 8 splits with private
  accumulators: Q pass 2-2.6x faster), a faster QP (primal-dual active set,
  Schur complement), 2-bit genotypes with a 256-entry table (memory only),
  `target_clones` builds for AVX2/AVX-512.
- Tried and not adopted: reusing Hessians within a QN iteration, sparse
  Hessians, warm-started QPs, K as a template parameter in the QP, capping
  warm-up epochs.

### Our kernels (`src/`)

- Genotypes are decoded from the memory-mapped FBM in every pass
  (`radmix_geno_access.h`: `_pMat[_rows[i] + _col_off[j]]`, missing table,
  `_lut[b] * _scale[i]`); `slice_geno` (`radmix_kernels.cpp:243-252`) is
  scalar and branchy. Q (K x I) and P (K x J) are double; each kernel makes a
  transposed copy of Q (L188-193) and `em`, `grad_q`, `update_q` allocate
  `one_minus(p)` per call (L175-179).
- Tiles 64 x 64 (L152-153; `TILE_I == LOG_BATCH^2` enforced by a
  `static_assert`, L615).
- `slice_coef<true>` (L288-294): 4 divisions per genotype (`a = g/x; b = h/y;
  wa = a/x; wb = b/y`); EM 2.
- `accum_q_block` (L314-378): per genotype K gradient updates and K(K+1)/2
  Hessian updates of `wa*c1 + wb*c2` (~40 multiply-adds at K = 5; the
  heaviest kernel). Its accumulators exceed L1 from K ~ 10.
- `accum_p_block` (L514-591): 2K + K(K+1)/2 short SIMD dot products per
  (locus, tile), each ending with a horizontal reduction.
- Segmented Q path when there are fewer than 32 blocks of individuals
  (L411-432; HGDP quick: 4 blocks x 8 segments, the last block only 17
  individuals); `grad_part`/`hess_part` allocated and zeroed per call
  (L452-453).
- Per iteration (`radmix_fit.cpp:229-236, 273, 302`): `update_q`,
  `update_p`, `update_q`, `update_p`, `loglik(x_qn)`: 5 passes, 6 on
  rejection; QN accepted in 16/16 iterations on HGDP quick.
- `radmix_quasi_newton.cpp`: `qn_push` serial (L63-80); `QN_CHUNK = 65536`
  (L56), so only 2 chunks on HGDP quick; stage 1 accumulates
  `sums[...] +=` through a pointer in the inner loop (L114-129) with false
  sharing of `partial`. All of this is far less than one pass.
- `src/Makevars`: only the OpenMP flags; R's `CXXFLAGS` here are `-g -O2`
  (GCC 13.3), so `omp simd` loops get SSE2 only, no FMA.

### perf profile of a fit (2026-10-31)

`perf record -F 1999 --call-graph fp` of `gt_radmix(k = 4, seed = 1,
n_cores = 1)` on EUR_subset (random start, 22 iterations, 30 s), installed
with `-g -O2`; `perf report --sort srcline --symbol-filter=<kernel>`.
Shares of all the samples of the fit (scripts in the session scratchpad:
`perf/perf_fit.R`):

| What | Where (`src/radmix_kernels.cpp` unless noted) | Share |
|---|---|---|
| Hessian accumulation, Q and P sides | L357 (`hk += wa*c1 + wb*c2`), L580 (`sh += wq*qk2`) | 23% |
| Likelihood coefficients (4 divisions per genotype) | `slice_coef`, L284-293 | 21% |
| Genotype decoding | `slice_geno` L248-251, `radmix_geno_access.h` L205-217 | 17% |
| Gradients | L343, L566 (+ `a -= b`, L558-559) | 11% |
| `qp = sum_k q_k p_k`, recomputed in every pass | `slice_geno` L230-233 | 9% |
| Log-likelihood (`loglik`, incl. its own `qp`) | L599-669 | 7% |
| `wq = wa*qk` (P side) | L574 | 2% |
| EM, QP solver and sweeps, quasi-Newton, R | | < 5% |

HGDP quick (K = 5, with missing data) has the same profile at the function
level (`accum_p_block` 24%, `accum_q_block` 19%, `slice_geno` 19%, `loglik`
4% of a process that also loads the data). What it suggests:
- the divisions (experiment 3a: one or two reciprocals instead of four
  divisions) are the largest share that is cheap to cut: up to ~15% of a
  fit;
- decoding is larger than expected (17%): a faster path for the common case
  (rows of the gen_tibble in the order of the backing file, no CV mask, all
  diploid) reading the bytes of a column directly, with two 256-entry
  tables giving g and h (0 for missing) instead of the `_rows`
  indirection, the missing table, the LUT and the scaling; or the
  pre-decoded copy of experiment 7;
- the Hessian accumulation (23%) is where 3b (ADAMIXTURE's expansion) and
  FMA (experiment 2, flags) would act;
- `qp` (9%) is recomputed in each of the 5 passes of an iteration; it is a
  small part of each pass, so only flags or FMA would help.

**Follow-up (2026-11-01, unsupervised overnight run):** experiments 2, 3a,
3b and a new decoding fast path were run, with two-loci accumulation and
x86-64-v3 multiversioning; the combination makes fits 1.6x faster (1.3-1.4x
without multiversioning). Results, patches and decisions in
`kernel_experiments/README.md` (next to this file).

### Suggestions already in `dev/`

- `dev/admixture_implementation.md` ("Possible further improvements", from
  `data-raw/radmix/Notes_to_do`): secant pairs in single precision (memory);
  `1 - P` kept in a reused buffer; tune the tile sizes and OpenMP schedule
  (never tuned since the port); for very large data, fit on a subset of loci
  and estimate P for the others in chunks (memory). `Notes_to_do` also: the
  optional in-memory copy of the subset was never built.
- `dev/simplification_plan.md`: tier 1 (prototyped, identical results, same
  speed, less memory): `1 - p` in place; the QN point written into the x0
  buffers (keeping the 5cf546f fix). Tier 3: drop `mu` from the Q tableau
  (rounding-level changes); R's LAPACK instead of RcppArmadillo. Not
  recommended: `SubBMCode256Acc` (19% slower); merging Q- and P-side
  accumulation.

## 3. Experiments, in order of priority

Each on its own branch/commit, measured alone (see section 4).

1. **CV genotype access fast path** (bit-identical; trivial; biggest CV win).
   In `src/radmix_geno_access.h` (L207-235): replace `z % _n_folds` with a
   `switch` over common fold counts (2-10, compile-time constants, so the
   compiler uses a multiplication instead of a 64-bit division; same folds,
   same results), and make the hide test branch-free (multiply the observed
   weight by 0/1 instead of branching). Expected: CV passes close to normal
   passes, CV time ~49 -> ~30 s on HGDP quick. Optionally (memory cost of 1
   byte per genotype) precompute the fold map once.

   **Done (2026-10-31, branch `radmix-cv-speed` from `dev`, not committed
   yet).** `mod_n_folds()` (switch over 2-10 folds) and a branch-free hide
   test in `GenoView::get()`; the decode loop of `slice_geno()` multiplies g
   and h by a 0/1 weight instead of branching. No memory added (the fold map
   was not built, so it scales to any size). Results bit-identical (all
   kernel outputs and a full CV fit, on quick and full HGDP); all 2443
   radmix tests pass. Benchmark: `data-raw/radmix/benchmark_cv_access.R`
   (kernels on 1 core in plain, training and held-out views, K = 5):
   training-view passes 19-25% faster on both quick (209 x 25,000) and full
   (1043 x 130,943) HGDP, i.e. from ~1.9x to ~1.5x a plain pass (the rest is
   the hash itself); a full fit with 5-fold CV on quick HGDP 60.5 -> 41.2 s.
   On the full data the plain log-likelihood pass also got faster (1.77 ->
   1.42 s; its missing genotypes no longer mispredict a branch).
2. **Compiler flags, measure only** (trivial). Rebuild with a personal
   `~/.R/Makevars` (`CXXFLAGS = -O3`, then add `-fno-trapping-math`, then
   `-march=native`) and time the kernels, to see what the existing `omp simd`
   loops gain. If large, look for a portable route (in-source
   `#pragma GCC optimize` or `__attribute__((target_clones(...)))` behind
   `__GNUC__`/ELF guards); `-march=native` cannot ship on CRAN.
3. **Cheaper Q-side kernel** (rounding-level changes; easy; 5-30% of a pass,
   to measure).
   a. Two reciprocals instead of four divisions in `slice_coef<true>`
      (`rx = 1/x; a = g*rx; wa = a*rx`, and the same for y). `a = g*(1/x)`
      alone is bit-identical (g, h are in {0, 1/2, 1, 2}), `wa` is not.
   b. ADAMIXTURE's expansion in `accum_q_block`
      (`adamixture/src/utils_c/cython/sqp.pyx::_grad_hess_q_block`):
      (1-p)(1-p)' = 11' - 1p' - p1' + pp', so one FMA per Hessian pair
      (`h[k,k2] += p_k p_k2 * ts`), then `XtX[k,k2] = H - A[k] - A[k2] + C`;
      ~27 instead of ~40 multiply-adds per genotype at K = 5. Check
      cancellation against the R reference kernels
      (`tests/testthat/helper-radmix-reference.R`).
4. **Fold refit policy** (trivial code; changes CV values; needs a
   statistical check). Try a looser fold `rtol` (e.g. 1e-5) and a cap of 5
   warm-started iterations (ADAMIXTURE's rule), and compare the CV errors and
   the chosen K with converged CV and with ADMIXTURE on the Anolis, HGDP quick
   and EUR_subset data, at several K (watch the leakage risk of short
   warm-started refits). Adopt only a rule that keeps the K ranking and the
   1e-3 agreement with ADMIXTURE of `test_radmix_vs_admixture.R`.
5. **SVD+ALS start, prototype in R first** (medium; fewer iterations).
   `radmix_fit()` already accepts `init_q`/`init_p` (and then skips EM), so
   the optimiser is unchanged: build a start from a randomised SVD (rank
   max(K+10, 20), Patterson-scaled, missing set to 0; e.g.
   `bigstatsr::big_randomSVD` on an imputed copy of the FBM) and a few ALS
   iterations (as in `adamixture/model/als.py`), and count iterations and
   time to the same optimum on HGDP quick and EUR_subset, at several K and
   seeds (also checking how often the best optimum is reached). Port to C++
   as an option (e.g. `init = "svd"`) only if it saves most of the 16 -> 6
   iterations.

   The start of ADAMIXTURE is a refinement of that of fastmixture (also our
   `gt_fastmixture()` default since `random_init = FALSE`); both are a
   randomised SVD of the centred genotypes (missing set to 0) followed by
   ALS on the factors, reconstructing the frequencies as 0.5 U S V' + f. The
   SVD routine is the same algorithm (oversampling to max(K+10, 20) columns,
   dynamically shifted power iterations, K-1 components kept; fastmixture
   `utils.py:82`, `randomSVD`). Differences:

   | | fastmixture 1.3.0 | ADAMIXTURE 1.8.9 |
   |---|---|---|
   | genotypes in the SVD | centred | centred and Patterson-scaled |
   | power iterations | 11, fixed | up to 5, early stop |
   | loci | a subset if M > 50,000 (70%, at most 500,000), the others projected | all |
   | ALS | unconstrained (`pinv`), then clipping P / normalising Q | bounded LS (BVLS), Fisher-weighted Q step, momentum 0.8 |
   | ALS stop | RMSE of Q < 1e-4, max 1000 | the same |

   Prototype fastmixture's version first (its `factorALS`, `utils.py:133`,
   is ~15 lines of matrix algebra, easy in R on top of
   `bigstatsr::big_randomSVD`), and try ADAMIXTURE's refinements only if it
   still needs many more iterations than ADAMIXTURE's 6. (fastmixture's 315
   iterations say nothing about its start: its iterations are EM-based and
   not comparable.)

   **First results (2026-10-31).** Script:
   `data-raw/radmix/explore_svd_start.R` in the repository (not committed):
   tidypopgen's `gt_pca_randomSVD()` (Patterson scaling, K - 1 components,
   default tol 1e-4) on genotypes imputed with `gt_impute_simple("mean0")`,
   then fastmixture's ALS (unweighted, clipping; 3 ALS seeds), then
   `gt_radmix()` from that start (no EM), against 5 random starts; 1 core,
   rtol 1e-7. Medians (`start` = SVD + ALS time; `at best` = runs within 0.1
   of the best log-likelihood of all the runs at that K):

   | data, K | random: iter / fit s / at best | SVD: iter / start s / fit s / at best |
   |---|---|---|
   | HGDP quick, 3 | 14 / 4.5 / 5 of 5 | 7 / 1.3 / 2.0 / 3 of 3 |
   | HGDP quick, 5 | 18 / 11.6 / 5 of 5 | 12 / 1.9 / 6.9 / 3 of 3 |
   | HGDP quick, 7 | 21 / 15.7 / 1 of 5 | 13 / 2.1 / 10.3 / 0 of 3 (all 7796 below) |
   | EUR_subset, 2 | 13 / 13.7 / 2 of 5 | 10 / 3.9 / 9.7 / 3 of 3 |
   | EUR_subset, 3 | 19 / 26.3 / 1 of 5 | 36 / 4.6 / 48.2 / 0 of 3 (5, 44, 44 below) |
   | EUR_subset, 4 | 21 / 36.1 / 1 of 5 | 11 / 7.8 / 22.2 / 0 of 3 (0.9-2.6 below) |
   | Anolis, 2-4 | 11-12 / 0.1-0.2 / 1 of 15 | 7-9 / 0.1 / 0.1 / 5 of 9 |

   Findings:
   - The start is much closer to the optimum (HGDP K = 5: log-likelihood
     -3.92e6 against -4.29e6 after a random start and EM) and usually cuts
     the iterations by 30-50% (total time, start included, ~25-40% less on
     HGDP).
   - It is nearly deterministic: the ALS seed hardly changes it (EUR K = 2:
     identical starts), so repeated runs explore few modes. At HGDP K = 7 all
     SVD runs found a mode 7796 below the best (found by 1 random run of 5);
     ADAMIXTURE (seeds 1-3) found that same mode once and another one, 645
     below, twice: SVD-type starts are steered away from the best mode there.
   - At EUR K = 3, our simple start converged slowly (36 iterations) to modes
     5-44 below the best, whereas ADAMIXTURE's refined start reached the best
     with all 3 seeds (in 24-46 s, slower than our random starts): its
     refinements (Fisher-weighted Q step, momentum 0.8, bounded least
     squares) seem to matter for the quality of the start.
   - The SVD itself is costly on EUR_subset (3.8-10 s, about a third of a
     fit), so its tolerance (gt_pca_randomSVD's `tol`, 1e-4) or
     gt_pca_partialSVD() should be tried.

   Next steps: (a) add ADAMIXTURE's ALS refinements and check EUR K = 3 and
   HGDP K = 7; (b) a cheaper SVD (looser tol, or the partial SVD); (c) more
   seeds, and a perturbation of the start (or random ALS initial P on more
   seeds) to see whether runs from SVD starts can be made diverse enough for
   the multi-run analyses we recommend with gt_clumppling(); (d) if adopted,
   probably as an option (e.g. `init = "svd"`, or the first run of several
   SVD-started and the others random), not as the default.

   **Refined starts (2026-10-31, branch `svd-start`).** Script:
   `data-raw/radmix/explore_svd_start_refined.R` (untracked). Changes:
   - SVD: `gt_pca_partialSVD()` gives the same subspace as
     `gt_pca_randomSVD()` (cosines of the principal angles 1.0000) ~5x
     faster (EUR rank 4: 1.3 vs 6.5 s; HGDP quick rank 7: 0.9 vs 2.6 s; full
     HGDP 1043 x 130,943 rank 7: 8.5 vs 45.5 s; loosening the tolerance of
     the random SVD only saves ~20-35%). Computed once per dataset at the
     largest K (as ADAMIXTURE does) and reused for all K and runs.
   - ALS "refined" as ADAMIXTURE (`model/als.py`): K components (not K - 1),
     least squares bounded in [0, 1] (here by cyclic coordinate descent over
     all rows at once), Fisher-weighted Q step (w = 1/max(f(1-f), 1e-3)),
     momentum 0.8 on Q for the P step. In R it costs 1-11 s (ADAMIXTURE's,
     compiled: ~0.2 s); on a random 25% of the loci (P of all loci then from
     Q in one P step) 0.6-2.4 s, with the same quality.
   Results (5 runs per strategy; iterations are medians; "at best" within
   0.1, "near" within 5 of the best of all runs at that K):

   | data, K | random: iter / at best / near | refined: iter / at best / near | refined, 25% of loci: iter / at best / near |
   |---|---|---|---|
   | HGDP quick, 5 | 18 / 5 / 5 | 6 / 4 / 4 | 6 / 3 / 3 |
   | HGDP quick, 7 | 21 / 1 / 1 | 9 / 0 / 0 | 9 / 0 / 0 |
   | EUR_subset, 3 | 19 / 1 / 3 | 12 / 0 / 3 | 11 / 0 / 2 |
   | EUR_subset, 4 | 21 / 1 / 5 | 10 / 0 / 5 | 13 / 0 / 5 |
   | Anolis, 3-4 | 11-12 | 7-8 | 7-9 |

   - The refined start cuts the iterations by 50-70% (HGDP K = 5: 6, as
     ADAMIXTURE; EUR K = 3, where the simple start needed 36: 12), and the
     fit time by about half to two thirds.
   - Variability: different ALS seeds converge to different ALS solutions
     (2 modes in 5 runs in most cases); subsampling the loci adds little.
     Perturbing the start towards a random one (30% or 60%) brings the
     iterations back to those of random starts (12-23) with no better
     exploration of modes: dropped.
   - Weakness: some ALS solutions lead to a poor mode. At HGDP K = 5, 1-2 of
     5 SVD starts end 19,808 below the best (Africa split in two clusters,
     Oceania merged with East Asia; no empty cluster), which random starts
     never do; at HGDP K = 7 and EUR K = 3 no SVD start reaches the best mode
     (found by 1 random start in 5; ADAMIXTURE misses it too at HGDP K = 7).
     The poor HGDP K = 5 start has a much worse starting log-likelihood
     (-3,899,022 vs -3,878,322), so choosing the best of a few ALS solutions
     by their starting log-likelihood (cheap: ALS, plus one pass each) may
     avoid it; but at EUR K = 3 the start that led to the worse mode had the
     better starting log-likelihood, so this is not a general guarantee.

   **More than K components** (refined ALS with K + 2, K + 5 and 2K
   components of the SVD; 5 runs each; median iterations / runs within 5 of
   the best):

   | data, K | K | K + 2 | K + 5 | 2K |
   |---|---|---|---|---|
   | EUR_subset, 3 | 12 / 3 (max gap 61) | 8 / 5 (1.3) | 7 / 5 (0.6) | 8 / 5 (1.3) |
   | EUR_subset, 4 | 10 / 5 | 11 / 5 | 10 / 5 | 10 / 4 (one 901 below) |
   | HGDP quick, 5 | 6 / 4 | 6 / 4 | 6 / 4 | 6 / 4 |
   | HGDP quick, 7 | 9 / 0 | 9 / 0 | 11 / 0 | 11 / 0 |

   A few extra components fix the poor starts at EUR K = 3 (weak structure
   among Europeans, where the top K components miss what the fit needs):
   all runs near the best, in 7-8 iterations instead of 12, with a better
   starting log-likelihood (-15,606,583 vs -15,609,805); the ALS costs about
   the same. They do not change HGDP K = 5 (the poor mode, splitting
   Africa, comes from the ALS solution, not from the rank) or K = 7. 2K is
   no better than K + 5, and once led to a poor mode (EUR K = 4); the
   extra components also reduce the variability (EUR K = 3: one mode). So
   K + a few (e.g. K + 5) is a cheap improvement; it does not solve the
   exploration of modes.

   **Best of n ALS seeds** (K + 5 components; each run chooses, among n
   ALS starts from different seeds, the one with the best starting
   log-likelihood, one pass over the genotypes each; 5 runs; median
   iterations / runs at best / modes found):

   | data, K | random | 1 SVD start | best of 3 | best of 5 |
   |---|---|---|---|---|
   | HGDP quick, 5 | 18 / 5 / 1 | 6 / 4 / 2 | 6 / 5 / 1 | 6 / 5 / 1 |
   | HGDP quick, 7 | 21 / 1 / 3 | 11 / 0 / 2 | 8 / 0 / 2 | 8 / 0 / 2 |
   | EUR_subset, 3 | 19 / 1 (3 near) / 3 | 7 / 0 (5 near) / 1 | 7 / 0 (5 near) / 1 | same |
   | EUR_subset, 4 | 21 / 1 (5 near) / 1 | 10 / 0 (5 near) / 1 | same | same |
   | Anolis, 4 | 12 / 0 / 2 (worst 1958) | 8 / 0 / 2 (worst 1918) | 8 / 2 / 1 (worst 11) | same |

   - Best of 3 removes the poor modes (HGDP K = 5: the Africa split; Anolis
     K = 4: a mode 1918 below) at the cost of n - 1 more ALS and n passes
     (negligible in compiled code; 0.2-50 s in this R prototype, which is
     dominated by the R ALS). Best of 5 adds nothing over best of 3.
   - It does not reach the best mode at HGDP K = 7 (as ADAMIXTURE: SVD+ALS
     starts seem not to lead there), and it removes most of the variability:
     the runs converge to the same mode.
   - Overall: K + 5 components, refined ALS, best of 3 is a fast (6-10
     iterations instead of 18-21) and reliable single start, but not a way
     to explore modes; for multi-run analyses (gt_clumppling()), random
     starts are still needed. A natural design: `init = "svd"` as an option
     for a fast, good solution, or, with several runs, the first run from
     the SVD start and the others random.

   Next steps: (a) done (best of 3, above);
   (b) a multi-run policy mixing SVD and random starts (e.g. one SVD start
   and the others random), judged by the chance of reaching the best mode
   and the total time for n runs; (c) if adopted, port the ALS to C++ (it
   needs no pass over the genotypes, so it would be negligible) and offer it
   as an option (`init = "svd"`), with the SVD computed once for all K and
   runs in `gt_radmix()`.

   **Decision and implementation (2026-10-31, branch `svd-start`,
   uncommitted).** (b) dropped: mixing SVD and random starts would bias the
   mode exploration of gt_clumppling() towards the modes the ALS favours.
   (c) done: `gt_radmix(init = "svd")`, an option only, documented as a fast
   single fit for large data or a first exploration, not for multi-run mode
   analyses. `R/radmix_svd.R` (SVD factors with a scaling that leaves
   monomorphic loci unscaled, which tidypopgen's default rejects; best of 3
   ALS solutions by log-likelihood) and `src/radmix_als.cpp` (Armadillo
   products, per-row solves with OpenMP: the unconstrained solution A^-1 b
   of all the rows at once, used where it lies in [0, 1]^K, coordinate
   descent for the others). Needs imputed genotypes (error otherwise) and
   more individuals and loci than max(k) + 5. Same quality and iterations as
   the R prototype (seeds 1-5); best-of-3 start time:

   | data, K | R prototype | C++, coordinate descent only | C++, with the direct solution |
   |---|---|---|---|
   | HGDP quick 5 | 6.4 s | 1.1 s | 0.7 s |
   | HGDP quick 7 | 20 s | 2.4 s | 2.1 s |
   | EUR 3 | 23 s | 3.9 s | 1.4 s |
   | EUR 4 | 31 s | 5.7 s | 2.0 s |

   (coordinate descent alone spent 90% of the ALS in the small solves, 8-10
   sweeps per row). EUR K = 4 per run: SVD 2.4 s (once per call; 10 s with
   the reference BLAS), start 2.0 s, fit 13.5 s (10 iterations), against
   34 s (21 iterations) for a random start.
6. **Small QN-loop tweaks** (easy; small gains).
   a. admixer's damped retries of a rejected QN step (0.5x, then 0.25x):
      helps hard data (ours accepted 16/16 steps on HGDP quick).
   b. Fold `loglik(x_qn)` into the next pass (5 -> 4 passes per accepted
      iteration; admixer measured ~4%).
   c. The `dev/simplification_plan.md` tier-1 buffer changes, as cleanup
      (no speed change expected).
7. **Later, only if 1-5 leave room**: admixer's mini-batch warm-up (needs P
   and Q updates restricted to a range of loci; changes the trajectories, so
   check how often the best optimum is reached); the pre-decoded in-memory
   copy of the genotypes (the `dev/` idea; 1 byte per genotype; removes the
   `_rows` indirection and per-read decoding, and makes CV masking free);
   tile and OpenMP schedule tuning (decouple `TILE_I` from `LOG_BATCH^2`
   first); BLAS or float Hessians for K >= 8.

   **BLAS / Armadillo / Eigen for the kernels (measured 2026-10-31).** The
   Q-side gradient and Hessian can be written as GEMMs per tile (coefficient
   matrix I x J times the J x K(K+1)/2 products of P, as admixer does). HGDP
   quick (209 x 25,000), 1 thread; our `grad_q` minus `loglik` (~ the
   accumulation) against the four equivalent GEMMs on dense matrices
   (`scratchpad/kernel_blas.R`):

   | K | our accumulation | GEMMs, OpenBLAS | GEMMs, reference BLAS |
   |---|---|---|---|
   | 3 | 0.007 s | 0.048 s | 0.097 s |
   | 5 | 0.032 s | 0.051 s | 0.190 s |
   | 10 | 0.13 s | 0.075 s | 0.574 s |

   So: BLAS only wins at K >= ~8, and only with an optimised BLAS (it brings
   AVX2/FMA kernels that our `-O2` build cannot use); with the reference
   BLAS (the default of CRAN's Windows and macOS builds of R) it would be
   3-7x slower. R's BLAS inside our OpenMP threads also needs it
   single-threaded (OpenBLAS is multithreaded by default; bigstatsr handles
   this with bigparallelr/RhpcBLASctl), and results would differ slightly
   between BLAS libraries. Eigen (header-only, RcppEigen in LinkingTo) avoids
   the threading and BLAS-dependence issues, but compiles with our flags
   (SSE2 only), so its gain over our loops is unclear: try it only together
   with the K >= 8 experiment above. Not worth it for K <= 7. (Armadillo is
   used for the products of the SVD start's ALS, where it does not matter:
   90% of the ALS time was in the small bounded least squares, see 5.)
   The partial SVD itself (bigstatsr) is BLAS-bound: EUR, 9 components,
   2.4 s with OpenBLAS against 10.1 s with the reference BLAS.

Not worth trying: avoiding per-entry `log()` (we already batch 8 logs per 64
genotypes); a looser stopping rule (ours is already the loosest of the
programs).

## 4. Protocol

- Install with `R CMD INSTALL --preclean .` before every timing:
  `load_all()`, `document()` and `devtools::test()` leave unoptimised `.o`
  files in `src/` that a plain `R CMD INSTALL .` reuses (it made gt_radmix ~6x
  slower once). Check that the install log shows `-O2` for `src/radmix_*.cpp`.
- Time on 1 core first (the laptop throttles with 4 busy cores), median of 3
  runs, with pauses between runs; record seconds, iterations, passes and the
  log-likelihood.
- Tools:
  - `data-raw/radmix/profile_eur.R`: per-kernel times (`radmix_loglik_cpp`,
    `grad_q`, `update_q`, `grad_p`, `update_p`) on EUR_subset, K = 4, at 1
    and 4 cores. Extend it with the CV view
    (`ns$radmix_set_cv(dat, list(mode = "train", fold = 1, n_folds = 5,
    seed = 1))`) and `deviance`, to measure the CV overhead directly.
  - `data-raw/radmix/benchmark_admixture.R`: EUR_subset, K = 2-4, radmix at
    1 and 4 cores and with ADMIXTURE's threshold, and ADMIXTURE; writes
    `benchmark_admixture.csv`.
  - `vignettes/articles/code/admixture_benchmark_radmix.R quick`: HGDP quick,
    1/2/4 cores, with and without CV (its times include reading the
    `gen_tibble` and `saveRDS`). The data (`data/hgdp/hgdp_bench_quick.*`)
    are made by `vignettes/articles/code/admixture_benchmark_data.R`.
  - EUR_subset is downloaded by `data-raw/radmix/download_eur_subset.R`
    (cached in `tools::R_user_dir("tidygenclust", "cache")` by the tests).
- Correctness gates:
  - `tests/testthat/test_radmix_vs_julia.R`: log-likelihood and EM 1e-12,
    `update_q`/`update_p` 1e-10, `expect_identical(fit$n_iter, ...)`,
    log-likelihood history 1e-10, Q and P 1e-8.
  - `tests/testthat/test_radmix_vs_admixture.R`: log-likelihood 1e-6, Q/P
    mean absolute difference < 0.005 and correlation > 0.999, CV error 1e-3,
    warm start <= 7 iterations, same CV ranking of K.
  - Thread invariance (`expect_identical`): `test_radmix_kernels.R:154-169,
    316-321`, `test_gt_radmix.R:66-71`, `test_radmix_cv.R:147`,
    `test_radmix_quasi_newton.R:57`; any reduction whose order depends on
    the number of threads breaks them.
  - Kernels against the R reference at 1e-10
    (`helper-radmix-reference.R:102`); deviance 1e-12; fold masks
    complementary and of similar sizes (`test_radmix_cv.R:37-65`).
  - Experiments 1-2 must pass unchanged; the rounding-level ones (3, 6b) may
    need the Julia-fixture tolerances or the identical iteration counts
    reviewed, decided per experiment.
- Record the results and decisions in `dev/admixture_implementation.md`
  ("Possible further improvements").

## 5. Sources

- ADAMIXTURE 1.8.9 (BSD-3; github.com/AI-sandbox/ADAMIXTURE, tag 1.8.9; also
  installed in the `ctidygenclust` conda environment, with the `.pyx`
  files): `model/br_qn.py`, `model/als.py`, `model/svd.py`,
  `src/adamixture.py`, `src/cv.py`, `src/utils_c/cython/{sqp,tools,br_qn,
  bvls}.pyx`. CPU performance commits: `8085960ad6` (kernels vectorised
  across samples, the (1-p)(1-p)' identity, missing as mask weights; one
  BR-QN iteration 7.5 -> 2.9 s on an Apple M5), `bd4e93a28c`, `4f82efbdc9`,
  `28cf4c7d2c`, `9573fa3427`; fold polishing from 3 to 5 iterations in
  `79bc777618`.
- admixer 0.2.13 (MIT; github.com/aalbrechtsen/admixer): `src/fit.hpp`
  (`Fitter::warmup`, `Fitter::qn`), `src/model.hpp`, `src/io.hpp`
  (`row_pass`), `src/smallk.hpp`, `src/linalg.hpp`, `src/kernel.hpp`,
  `bench/perf_results.md`, `DETAILS.md`. Its benchmark figures are the
  author's own; the repository was a month old when read.

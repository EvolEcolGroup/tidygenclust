> Copied on 2026-11-01 from the laptop notes (Dropbox, kernel_experiments/README.md). Paths to /tmp worktrees and the scratchpad refer to the laptop and no longer exist: the patches are in ../patches, the results (CSV and summaries; the .rds outputs stayed on the laptop) in ../laptop_results, and the scripts, made portable, in ../scripts. Laptop label "hgdp" = dataset "hgdp_quick" here.

# Kernel speed-up experiments (overnight, 2026-10-31 to 2026-11-01)

Four experiments from the perf profile in `../radmix_speedup_plan.md`
("perf profile of a fit"), run without supervision. This file holds the plan,
the protocol and the results, updated after each step; the summary at the top
is written last.

## Summary

**`gt_radmix()` fits are 1.6x faster** on this laptop (i7-10510U, AVX2):
0.61-0.63 of the time on 1 core, 0.47-0.65 on 4 cores, 0.66 with
cross-validation, with the same iterations, log-likelihoods (to ~1e-8) and
CV errors. Where the AVX2 versions are not built (other platforms or
compilers, or the opt-out), 1.3-1.4x (0.71-0.76 on 1 core, 0.81 with CV).

| fits, seed 1, rtol 1e-7 (seconds) | baseline | combination | without multiversioning |
|---|---|---|---|
| EUR K = 3, 1 core | 19.8 | 12.2 (0.62) | 14.0 (0.70) |
| EUR K = 4, 1 core | 30.6 | 18.7 (0.61) | 22.4 (0.73) |
| HGDP K = 5, 1 core | 6.4 | 4.0 (0.63) | 4.9 (0.76) |
| HGDP K = 7, 1 core | 12.9 | 7.9 (0.61) | 9.8 (0.76) |
| EUR K = 3, 4 cores | 7.4 | 3.5 (0.47) | 4.4 (0.60) |
| EUR K = 4, 4 cores | 13.2 | 8.5 (0.65) | 9.5 (0.72) |
| HGDP K = 5, 4 cores | 2.3 | 1.2 (0.52) | 1.5 (0.63) |
| HGDP K = 7, 4 cores | 5.3 | 3.1 (0.59) | 4.1 (0.78) |
| HGDP K = 5, CV (5 folds), 1 core | 31.2 | 20.5 (0.66) | 25.3 (0.81) |

The combination (worktree/branch `kernel-speedups`, uncommitted;
`kernel_speedups_combined.patch`, only `src/radmix_kernels.cpp` and
`src/radmix_geno_access.h`):

| change | experiment | gain alone (fits) | results |
|---|---|---|---|
| two reciprocals instead of four divisions | 1a | 2-5% | rounding |
| genotype decoding fast path (also for CV views) | 2 | 5-14% | bit-identical |
| two loci per pass in the Q and P accumulations | 4c + 4d | 6-14% | rounding |
| x86-64-v3 (AVX2 + FMA) multiversioning, Linux + GCC >= 12 | 3 | 14-29% | rounding |

Rejected: one reciprocal (1b), ADAMIXTURE's Hessian expansion (4a: slower,
the loop is limited by memory traffic), longer P-side tiles (4b),
`-fno-trapping-math` (no effect), multiversioning the log-likelihood (its
AVX2 version is slower).

Checks: the full test suite passes on both code paths of the combination
(AVX2: 3586 expectations, with 19 clones in the test build; baseline, with
`-DRADMIX_NO_MULTIVERSION`: 3586, and the radmix tests again: 2481);
R CMD check (no vignettes): 1 NOTE, a hidden `.git` file, which is an
artefact of building from a git worktree (where `.git` is a file); "checking
compiled code" is OK.

**Decisions for you:**
1. *Multiversioning and reproducibility.* With it, the same seed gives the
   same results on a given machine, but results can differ by rounding
   between machines with and without AVX2 (today they are bit-identical on
   any x86-64 machine). Compiling with `-DRADMIX_NO_MULTIVERSION` restores
   that. Is the 14-29% worth this? (I think yes, with a sentence in the
   documentation.) It also only covers x86-64 Linux with GCC >= 12: clang
   (macOS) and Windows (no ifunc) get the baseline version.
2. *Default `rtol` 1e-8.* With the combination, `rtol = 1e-8` runs in
   0.65-0.74 of today's time at 1e-7, and gains 0.02-2.6 log-likelihood units
   (EUR K = 3: Q within 0.009 of a fully converged fit instead of 0.05). The
   speed-ups pay for the stricter threshold with room to spare.
3. *How to bring it in*: the worktree (and `kernel_speedups_combined.patch`)
   also has a NEWS entry and a paragraph in `?gt_radmix` (end of
   "Algorithm") on point 1 and the opt-out, so it can be committed as it
   is, on `svd-start` (or applied on dev once PR #31 and the SVD start are
   merged; the patch applies to a0097bf).

Leads not pursued: `-O3` is 5-10% faster on the Newton kernels with
bit-identical results (from inlining or unswitching, not vectorisation;
worth finding which loops); a vectorisable `log()` would let the
log-likelihood use AVX2 too; the CV refit policy (experiment 4 of the
speed-up plan) remains the biggest lever for CV.

## Plan

Base: branch `svd-start` at a0097bf (SVD start committed). Each experiment
gets its own git worktree and branch (`exp-divisions`, `exp-decode`,
`exp-flags`, `exp-hessian`), created from a0097bf, with its changes left
**uncommitted** (nothing is committed or pushed); each finished
implementation is also saved here as a patch (`<experiment>.patch`). The
adopted ones are combined at the end in a worktree `kernel-speedups`, again
uncommitted, with a combined patch.

1. **Fewer divisions** (`slice_coef`, 21% of a fit). Variants: (a) two
   reciprocals (`rx = 1/x; a = g*rx; wa = a*rx`, same for y); (b) one
   reciprocal (`r = 1/(x*y)`, `rx = y*r`, `ry = x*r`). Rounding-level
   changes.
2. **Faster decoding** (`slice_geno` + `GenoView::get`, 17%). A fast path
   for the common case (rows of the gen_tibble contiguous in the backing
   file, no cross-validation mask, uniform ploidy): read the bytes of the
   column directly, with two 256-entry tables giving g and h (0 for
   missing). Bit-identical results expected. If time allows, an arithmetic
   (vectorisable) decoding for files without imputed codes.
3. **Compiler flags** (measure only, then a portable route if worth it):
   `-O2` (baseline), `-O3`, `-O2 -fno-trapping-math`, `-O3 -march=native`,
   `-O3 -march=native -fno-trapping-math`, set through a private
   `R_MAKEVARS_USER` file (the user's `~/.R/Makevars` is not touched). If
   `-march=native` gains a lot, prototype function multiversioning
   (`target_clones("avx2","default")`, Linux/GCC only, guarded) on the hot
   kernels and check that it captures the gain.
4. **Hessian work** (23%): (a) ADAMIXTURE's expansion on the Q side
   ((1-p)(1-p)' = 11' - 1p' - p1' + pp': one multiply-add per pair plus K + 1
   per genotype, instead of two per pair); check cancellation against the R
   reference kernels; (b) if time allows, a fused P-side pass (all the
   gradient and Hessian sums of a locus in one pass over the individuals,
   instead of K + K(K+1)/2 short dot products, each with a horizontal
   reduction).
5. **Combine** the variants that pay off, run the full test suite and
   R CMD check, and benchmark end to end against the baseline (1 and 4
   cores, with and without cross-validation), also at `rtol = 1e-8` to
   inform the decision on the default threshold (not changed here).

## Protocol

- Each variant (and the baseline) is installed with
  `R CMD INSTALL --preclean -l <private library>` (never into the user's
  library; checks that the install log shows `-O2`), and benchmarked by
  loading `tidygenclust` from that library.
- Kernel benchmark (`bench_kernels.R`): `radmix_loglik_cpp`, `radmix_em_cpp`,
  `radmix_update_q_cpp`, `radmix_update_p_cpp` (and the training view of
  cross-validation), on EUR_subset (K = 4) and HGDP quick (K = 5 and 10, with
  missing data), fixed Q and P, 1 core, median of 5.
- Fit benchmark (`bench_fits.R`): `gt_radmix()` on EUR_subset (K = 3, 4) and
  HGDP quick (K = 5, 7), seed 1, 1 core, random start; seconds, iterations,
  log-likelihood, and the largest differences in Q and P from the baseline.
  Baseline and variant runs are interleaved (A B A B), with pauses, to limit
  the effect of thermal throttling on this laptop (i7-10510U).
- Correctness: the radmix tests (`devtools::test(filter = "radmix")`) in
  each worktree; the full suite and R CMD check on the combination.
- Adopt a change if it passes the tests and gains at least ~3% end to end.

## Results

### 1. Fewer divisions (done)

`results/e1/` (kernel medians of 3 runs x 5 repetitions; fits: medians of 2
runs, seed 1, 1 core; ratio = variant / baseline):

| | update_q | update_p | em, loglik | fits EUR 3 / EUR 4 / HGDP 5 / HGDP 7 |
|---|---|---|---|---|
| (a) two reciprocals | 0.91-0.95 | 0.95-1.00 | 0.98-1.02 | 0.97 / 0.96 / 0.95 / 0.99 |
| (b) one reciprocal | 0.94-0.97 | 0.95-0.98 | 0.99-1.04 | 0.98 / 1.04 / 0.99 / 0.99 |

Both change the results only by rounding (Q within 2e-12, the same
iterations and log-likelihoods to 4e-9). The gain is far below the 21% that
perf attributed to the lines of `slice_coef`: the division latency there was
mostly hidden behind other work (perf charges stalls to the waiting
instructions). (b) is no better than (a) and would turn the (theoretical)
case qp = 0 into NaN, so **(a) is kept**: `divisions.patch` (two
reciprocals in `slice_coef<true>`, the EM path unchanged; a and b stay
bit-identical, only wa and wb change by rounding). Gain ~3-4% of a fit.

### 2. Faster decoding (done)

`results/e24/` (same protocol; the run also re-measured the divisions:
fits 0.94-0.98, confirming 1a). Fast path in `slice_geno()`: when the
individuals of a block are adjacent in the backing file and decoded alike
(no pseudohaploid mixture), the bytes of the block are read directly and
decoded with two 256-entry tables (g, and h = n - g; 0 for missing).
Contiguity is precomputed per individual (`_contig`), so filtered datasets
use it wherever a block happens to be contiguous.

| | em | loglik | update_q | update_p | CV views | fits EUR 3 / EUR 4 / HGDP 5 / HGDP 7 |
|---|---|---|---|---|---|---|
| v1 (no CV mask) | 0.82-0.92 | 0.87-0.92 | 0.91-0.97 | 0.92-1.00 | 1.00 (not used) | **0.86 / 0.89 / 0.92 / 0.95** |

All outputs bit-identical (kernels and fits: Q, P and log-likelihood
differences exactly 0). Radmix tests pass. v2 extends the fast path to
cross-validation views (the bytes are read directly, and the fold mask
multiplies the table values, as in the general loop): see below. Patches:
`decode_v1.patch`, then `decode.patch` (v2).

### 3. Compiler flags (done) and multiversioning (in progress)

`results/e3/` (kernels only; medians of 3 runs x 5 repetitions; ratio to
the `-O2` baseline; flags set through a private `R_MAKEVARS_USER`):

| | update_q | update_p | em | loglik | CV views (update_q/p) |
|---|---|---|---|---|---|
| `-O3` | 0.90-0.96 | 0.93-0.94 | 1.01-1.11 | 0.90-0.94 | 0.94-0.99 |
| `-O2 -fno-trapping-math` | 0.96-0.99 | 0.98-1.00 | 0.97-0.99 | 0.98-1.00 | 0.97-1.02 |
| `-O2 -mavx2 -mfma` | **0.66-0.72** | **0.67-0.73** | 0.76-0.84 | 0.76-0.79 | 0.71-0.84 |
| `-O2 -march=native` | 0.63-0.68 | 0.67-0.70 | 0.71-0.78 | 0.81-0.84 | 0.68-0.83 |
| `-O3 -march=native` | 0.61-0.64 | 0.63-0.69 | 0.63-0.67 | 0.76-0.79 | 0.66-0.79 |
| `-O3 -march=native -fno-trapping-math` | 0.62-0.64 | 0.64-0.69 | 0.67-0.69 | 0.76-0.80 | 0.68-0.84 |

- **AVX2 + FMA is the big lever** (Newton kernels 28-34% faster), and
  `-mavx2 -mfma` alone gets nearly all of `-march=native`. Outputs differ from
  the baseline by rounding (up to ~2e-15 in the kernels: FMA and 4-wide
  reductions).
- `-O3` gains 5-10% on the Newton kernels with bit-identical outputs, but
  slows EM. It is not from vectorisation: compiling with
  `-fopt-info-vec-optimized` shows that the extra loops vectorised at `-O3`
  are tiny K-length loops and one-off passes (`vecinfo/`), so the gain is
  from inlining/unswitching. A lead for later, not pursued (the
  `optimize` attribute/pragma is not meant for production code).
- `-fno-trapping-math` does nothing here (our loops are already written
  without branches).

Since only R's flags can be used to build the package, the AVX2 + FMA gain
can only be had from the source: function multiversioning.
`__attribute__((target_clones("arch=x86-64-v3", "default")))` (GCC >= 12)
compiles each marked function twice, and an ifunc resolver picks the
version for the processor at load time, from its feature bits (x86-64-v3 =
AVX2, FMA, BMI1/2, ...; not by CPU model, unlike `arch=haswell`). Checked on
a toy example with GCC 13.3: the OpenMP-outlined parallel regions are cloned
too (`._omp_fn.0.arch_x86_64_v3`, with ymm registers and vfmadd). Prototype
`multiversion.patch`: a macro `RADMIX_MULTIVERSION` on the 7 kernel entry
points and the 3 block workers of `radmix_kernels.cpp`, enabled only for
GCC >= 12 on x86-64 Linux (elsewhere it expands to nothing, so nothing
changes). The built library has 29 `arch_x86_64_v3` clones.

### 4. Hessian work (a, b rejected; c + d adopted)

(a) ADAMIXTURE's expansion of the Q-side sums, and (b) larger P-side tiles
(`results/e24/`, labels h_tp64 = (a), h_tp128 and h_tp256 = (a) + (b)):

| | update_q | update_p | fits EUR 3 / EUR 4 / HGDP 5 / HGDP 7 |
|---|---|---|---|
| (a) | 1.05-1.10 (0.97 at K = 10) | ~1.00 | 1.06 / 1.04 / 1.01 / 1.01 |
| (a) + P tile 128 | 1.05-1.09 | 1.01 | 1.07 / 1.04 / 1.02 / 1.01 |
| (a) + P tile 256 | 1.05-1.07 | 1.03-1.06 | 1.08 / 1.09 / 1.05 / 1.05 |

Both are slower. (a) does fewer multiply-adds but more passes over the
accumulator vectors (K + 2 more per locus), so this loop is limited by the
loads and stores of the accumulators, not by the arithmetic; it only gains at
K = 10. Longer P-side tiles are slower too (the TILE_I x K corner of Q no
longer stays in L1). Rejected (`hessian_expansion_REJECTED.patch`, kept for
reference). Radmix tests passed with (a), so the cancellation of the expanded
sums was not a problem.

(c) follows from this: accumulate two loci per pass on the Q side (the
coefficients of loci j and j + 1 first, then one load and store of each
accumulator for both), halving the accumulator traffic. (d) does the same on
the P side, where the sums are scalars: the sums of loci j and j + 1 are
computed in the same short loops over the individuals of a tile, sharing the
loads of q and halving the number of loops (each ends with a horizontal
reduction). Both loops with two reductions are vectorised at `-O2` (checked
with `-fopt-info-vec`).

`results/e4c/` (labels h_pairs0 = (c), h_pairs1 = (c) + (d); mv = the
multiversioning prototype of experiment 3, alone):

| | update_q | update_p | em | loglik | CV views | fits EUR 3 / EUR 4 / HGDP 5 / HGDP 7 |
|---|---|---|---|---|---|---|
| (c) | 0.87-0.91 | 0.98-1.01 | 0.97-0.98 | 1.00 | 0.88-1.01 | 0.98 / 1.01 / 0.91 / 0.94 |
| (c) + (d) | 0.88-0.92 | **0.83-0.84** | 0.97-0.98 | 0.98-1.00 | 0.85-0.96 | **0.92 / 0.94 / 0.86 / 0.86** |
| multiversioning | **0.68-0.74** | **0.72-0.75** | 0.85-0.87 | **1.07-1.09** | 0.72-0.81 | **0.83 / 0.86 / 0.75 / 0.71** |

All rounding-level (Q within 5e-12, same iterations and log-likelihoods to
1e-8); radmix tests pass for both. **(c) + (d) adopted** (`pairs.patch`):
6-14% faster fits.

Multiversioning: 14-29% faster fits, but the AVX2 log-likelihood is 7-9%
slower (it calls the C library's scalar `log()` inside the loop, which costs
AVX-SSE transitions), so in the combination `loglik()` and `deviance()` stay
baseline-only. Radmix tests pass on the AVX2 path (including the
OpenADMIXTURE.jl fixtures with identical iteration counts).

## Combination

Worktree `kernel-speedups` (from a0097bf, uncommitted):
divisions (1a) + decoding v2 (2) + multiversioning (3, without `loglik()` and
`deviance()`) + two loci per pass (4c + 4d), plus an opt-out:
compiling with `-DRADMIX_NO_MULTIVERSION` (e.g. in `~/.R/Makevars`) builds
only the baseline version, whose results are bit-identical on every x86-64
machine. Patch: `kernel_speedups_combined.patch` (the two C++ files, plus a
NEWS entry and a paragraph in `?gt_radmix` added at the end).

Final benchmark (`results/final/summary.txt`; labels base, comb, comb_nomv =
the combination built with `-DRADMIX_NO_MULTIVERSION`): see the summary at the
top for the fits. Kernels, 1 core, ratio to the baseline:

| | update_q | update_p | em | loglik | CV views |
|---|---|---|---|---|---|
| comb | 0.57-0.58 | 0.58-0.65 | 0.74-0.84 | 0.92-1.02 | 0.60-0.67 |
| comb_nomv | 0.75-0.82 | 0.69-0.74 | 0.79-0.89 | 0.85-0.90 | 0.76-0.86 |

(the log-likelihood is not multiversioned; its ratio in `comb` is probably
lowered clock speed after the AVX2 kernels that precede it in the benchmark.)

Stricter threshold, `rtol = 1e-8` (seed 1, 1 core):

| | baseline 1e-7 | baseline 1e-8 | comb 1e-8 | log-lik gain from 1e-8 |
|---|---|---|---|---|
| EUR K = 3 | 19.8 s, 19 it. | 24.4 s, 23 it. (1.23) | 14.6 s (0.74) | +2.62 |
| EUR K = 4 | 30.6 s, 22 it. | 33.1 s, 24 it. (1.08) | 20.5 s (0.67) | +0.47 |
| HGDP K = 5 | 6.4 s, 16 it. | 7.4 s, 18 it. (1.15) | 4.5 s (0.70) | +0.07 |
| HGDP K = 7 | 12.9 s, 21 it. | 13.6 s, 22 it. (1.05) | 8.4 s (0.65) | +0.02 |

Tests: the full suite passes with the AVX2 versions (3586) and with the
baseline versions (3586; radmix tests again with `-DRADMIX_NO_MULTIVERSION`
in a fresh build: 2481). R CMD check: 1 NOTE (`.git`, from the worktree).

## Where things are

- Patches and results: this folder (`*.patch`, `results/`).
- Worktrees (uncommitted changes; under the session scratchpad in /tmp, so
  they do not survive a reboot, but the patches here do):
  `/tmp/claude-1000/-media-andrea-annex-git-tidygenclust/01fb417f-4d4c-4e63-92ee-6ee7e1f0a402/scratchpad/wt/`
  (`divisions`, `decode`, `flags` = multiversioning prototype, `hessian` =
  two loci per pass, `combined` = branch `kernel-speedups`). The branches
  `exp-divisions`, `exp-decode`, `exp-flags`, `exp-hessian` and
  `kernel-speedups` all point at a0097bf, with no commits; remove them with
  `git worktree remove <path>` and `git branch -D <branch>` (or
  `git worktree prune` after a reboot).
- Scripts: `bench_kernels.R`, `bench_fits.R`, `bench_fits2.R`, `compare.sh`,
  `final.sh`, `summarise.R`, `final_summary.R` in the scratchpad's `exp/`
  folder (copied to `scripts/` here).

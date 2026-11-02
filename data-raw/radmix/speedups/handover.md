> Copied from `data-raw/radmix/perf_verification/HANDOVER.md` on the branch `radmix-perf-verify` (commit f423896),
> which also holds the scripts, patches and results (CSV) it refers to.

# gt_radmix() speed-ups: handover for verification on a bigger machine

Written on 2026-11-01 at the end of the laptop work (Intel i7-10510U, 4
cores, 15 W, AVX2, no AVX-512; R 4.6.1, GCC 13.3, OpenBLAS 0.3.26), for a
Claude instance that will verify it on a desktop with many cores.
Read sections 1-3 for context, then follow sections 4-7. Section 8 lists the
files of this kit.

## 1. Branches

| branch | commit | content | to merge? |
|---|---|---|---|
| `svd-start` | 7baa46c, on `dev` (ddcea53) | `gt_radmix(init = "svd")`: SVD start (R/radmix_svd.R, src/radmix_als.cpp, tests, docs, NEWS) | yes, first (PR to `dev`) |
| `kernel-speedups` | 017e176, on `svd-start` | the combined kernel speed-ups (src/radmix_kernels.cpp, src/radmix_geno_access.h, NEWS, `?gt_radmix`) | yes, after `svd-start`, once the decisions in 2.4 are taken |
| `radmix-perf-verify` | on `kernel-speedups` | this kit (`data-raw/radmix/perf_verification/`), and the results of the verification | no: the coordinating branch, for the record |

`data-raw/` is in `.Rbuildignore`, so the kit never reaches the package.
The individual experiments are not branches: they are patches
(`patches/`), which `scripts/build_variants.sh` applies to `svd-start`.

## 2. What was done on the laptop

### 2.1 SVD start (`init = "svd"`, branch `svd-start`)

Why: ADAMIXTURE (Saurina-i-Ricos et al. 2026) needs 6 iterations where
`gt_radmix()` needs 16-21 from a random start, at the same cost per
iteration; the difference is its start (randomised SVD + ALS, as in
fastmixture). Design (decided with the user):
- An option only. Mixing SVD and random starts in one analysis would bias the
  exploration of modes by `gt_clumppling()` towards the modes the ALS
  favours, so random starts stay the default.
- The SVD (`tidypopgen::gt_pca_partialSVD()`, ~5x faster than the random
  SVD with the same subspace) is computed once per call, with max(k) + 5
  components, and centring 2f, scaling sqrt(2f(1-f)) (1 for monomorphic
  loci); it needs imputed genotypes (error otherwise; the fit itself still
  ignores the imputed values unless `use_imputed = TRUE`).
- Each run: 3 ALS solutions (ADAMIXTURE's refined ALS, ported to C++:
  Fisher-weighted Q step, momentum 0.8, bounded least squares by the direct
  solution where it is inside the bounds, coordinate descent otherwise),
  each from random allele frequencies; the one with the best log-likelihood
  is the start; no EM steps.
- Error if min(individuals, loci) <= max(k) + 5.

Laptop results (1 core, medians of 5 seeds, before the kernel speed-ups;
"first run" includes the SVD, "further runs" share it; scripts and CSVs in
`laptop_results/svd/` and `notes/speedup_plan.md`, experiment 5):

| data, K | random: time (iterations) | SVD: SVD + start + fit (iterations) | first run | each further run |
|---|---|---|---|---|
| HGDP quick 5 | 6.5 s (18) | 0.9 + 0.7 + 2.5 s (6) | 4.1 s, 1.6x faster | 3.2 s, 2.0x |
| HGDP quick 7 | 11.9 s (21) | 0.9 + 1.8 + 5.2 s (8) | 7.9 s, 1.5x | 7.0 s, 1.7x |
| EUR 3 | 22.7 s (19) | 2.2 + 1.4 + 7.7 s (7) | 11.3 s, 2.0x | 9.1 s, 2.5x |
| EUR 4 | 34.3 s (21) | 2.2 + 2.0 + 14.1 s (10) | 18.3 s, 1.9x | 16.1 s, 2.1x |
| Anolis 3 | 0.12 s (11) | 0.45 + 0.03 + 0.08 s (7) | 0.56 s, 4.6x slower | 0.11 s, same |
| Anolis 4 | 0.16 s (12) | 0.45 + 0.05 + 0.12 s (9) | 0.62 s, 3.9x slower | 0.17 s, same |

Quality, against the best log-likelihood of any run at that K:

| data, K | random starts | SVD starts |
|---|---|---|
| HGDP quick 5 | 5/5 at the best | 5/5 at the best |
| HGDP quick 7 | 1/5 at the best, worst 3,000 below | 0/5: all in a mode 7,800-8,100 below (ADAMIXTURE misses it too) |
| EUR 3 | 3/5 within 5, worst 63 below | 5/5 within 0.55 |
| EUR 4 | 5/5 within 1.5 | 5/5 within 1.7 |
| Anolis 3 | 0/5 within 5, worst 56 below | 5/5 within 2.8 |
| Anolis 4 | 0/5 within 5, worst 2,000 below | 1/5 within 5, worst 7.8 below |

Conclusions so far: 1.5-2.5x faster single fits on medium data, usually as
good or more consistent, but SVD starts can systematically miss the best
mode (HGDP K = 7) and explore few modes, hence option only. The partial SVD
is BLAS-bound (EUR, 9 components: 2.4 s with OpenBLAS, 10.1 s with the
reference BLAS). On tiny data the SVD costs more than a fit; the user
suggested a minimum size for the option (only the max(k) + 5 check exists).
Until this kit was prepared, it had not been tested beyond 209 x 25,000 and
379 x 54,051; see 2.5 for the first run on the full HGDP (and the bug it
found).

### 2.2 Convergence threshold (`rtol`, not changed)

`gt_radmix()` stops when an iteration improves the log-likelihood by less
than rtol x |loglik| (default 1e-7: 1.6 units per iteration on EUR).
ADMIXTURE uses an absolute 1e-4, ADAMIXTURE 0.1. Seed 1, 1 core, laptop,
before the kernel speed-ups:

| data, K | rtol 1e-7 (now) | 1e-8 | ADMIXTURE's 1e-4 absolute |
|---|---|---|---|
| EUR 3 | 21.9 s, 19 it., 2.7 below, Q max diff 0.050 | 26.2 s (+20%), 23 it., 0.08 below, 0.009 | 32.8 s (+50%), 29 it. |
| EUR 4 | 30.9 s, 22 it., 0.64 below, 0.026 | 33.7 s (+9%), 24 it., 0.17 below, 0.015 | 45.8 s (+48%), 31 it. |
| HGDP quick 5 | 7.1 s, 16 it., 0.07 below, 0.0001 | 8.0 s (+13%), 18 it., 0.00 below | 8.4 s (+18%), 19 it. |

("below" = log-likelihood short of the fit at ADMIXTURE's threshold; Q
difference from that fit.) Recommendation: 1e-8 (pending the user's
decision). On the Anolis, even 1e-10 changes nothing: runs there end in
genuinely distinct local optima, not in unconverged ones.

### 2.3 Kernel speed-ups (branch `kernel-speedups`)

perf profile of a fit (EUR K = 4, 1 core): Hessian accumulation 23%,
divisions in `slice_coef` 21%, genotype decoding 17%, gradients 11%, qp 9%,
log-likelihood 7%. Four experiments (details: `notes/kernel_experiments_laptop.md`;
results: `laptop_results/`):

| experiment (patch) | change | fits, ratio to baseline (EUR 3 / EUR 4 / HGDP quick 5 / 7) | outputs | decision |
|---|---|---|---|---|
| 1a `divisions.patch` | two reciprocals instead of four divisions | 0.94-0.98 | rounding | adopted |
| 1b (no patch kept) | one reciprocal | 0.98-1.04 | rounding, NaN risk | rejected |
| 2 `decode.patch` | decoding fast path: bytes read directly, two 256-entry tables (also for CV views) | 0.86 / 0.89 / 0.92 / 0.95 | bit-identical | adopted |
| 3 flags (`-O3`, `-mavx2 -mfma`, `-march=native`) | measurement only | Newton kernels 0.66-0.73 with AVX2+FMA; -O3 0.90-0.96 | rounding (-O3: identical) | diagnostic |
| 3 `multiversion.patch` | `target_clones("arch=x86-64-v3", "default")` on the kernels | 0.83 / 0.86 / 0.75 / 0.71 | rounding | adopted (without loglik/deviance: their AVX2 versions were 7-9% slower) |
| 4a `hessian_expansion_REJECTED.patch` | ADAMIXTURE's Hessian expansion | 1.01-1.06 | rounding | rejected (the loop is bound by the loads/stores of the accumulators) |
| 4b (no patch) | longer P-side tiles | 1.02-1.09 | | rejected |
| 4c+4d `pairs.patch` | two loci per pass, Q and P sides | 0.92 / 0.94 / 0.86 / 0.86 | rounding | adopted |

The combination (1a + 2 + 3 + 4c/d; `-DRADMIX_NO_MULTIVERSION` builds only
the baseline version), laptop, seed 1, rtol 1e-7:

| fits (seconds) | baseline | combination | without multiversioning |
|---|---|---|---|
| EUR K = 3, 1 core | 19.8 | 12.2 (0.62) | 14.0 (0.70) |
| EUR K = 4, 1 core | 30.6 | 18.7 (0.61) | 22.4 (0.73) |
| HGDP quick K = 5, 1 core | 6.4 | 4.0 (0.63) | 4.9 (0.76) |
| HGDP quick K = 7, 1 core | 12.9 | 7.9 (0.61) | 9.8 (0.76) |
| EUR K = 3, 4 cores | 7.4 | 3.5 (0.47) | 4.4 (0.60) |
| EUR K = 4, 4 cores | 13.2 | 8.5 (0.65) | 9.5 (0.72) |
| HGDP quick K = 5, 4 cores | 2.3 | 1.2 (0.52) | 1.5 (0.63) |
| HGDP quick K = 7, 4 cores | 5.3 | 3.1 (0.59) | 4.1 (0.78) |
| HGDP quick K = 5, CV (5 folds), 1 core | 31.2 | 20.5 (0.66) | 25.3 (0.81) |

Same iterations, log-likelihoods within ~1e-8, Q within ~1e-11, identical
CV errors. rtol 1e-8 with the combination: 0.65-0.74 of the baseline's time
at 1e-7 (log-likelihood +0.02 to +2.6). Full test suite passed on both code
paths (3586 expectations); R CMD check: 1 NOTE (hidden `.git`, an artefact
of building from a worktree). The 4-core numbers are unreliable: the
laptop throttles (~60% slower after 10-20 s with 4 busy cores).

Leads not pursued: `-O3` gains 5-10% on the Newton kernels with identical
results (inlining, not vectorisation; which loops?); a vectorisable `log()`
would let the log-likelihood use AVX2; admixer-style ideas (mini-batch
warm-up, log-likelihood folded into the P step); the CV refit policy
(experiment 4 of `notes/speedup_plan.md`), the biggest lever for CV.

### 2.4 Decisions pending with the user

1. Multiversioning: results reproducible on a machine, but different by
   rounding between machines with and without AVX2 (opt-out
   `-DRADMIX_NO_MULTIVERSION`); only on x86-64 Linux with GCC >= 12.
2. Default `rtol` 1e-8 instead of 1e-7.
3. Merging `svd-start` and `kernel-speedups`.
The verification should give the evidence for these: report, do not decide.

### 2.5 Found while preparing this kit: the SVD start on the full HGDP

The first run of `svd_vs_random.R` on the full HGDP (1043 x 130,943, K = 5,
seeds 1-2, 1 core, `comb`) found a bug, now fixed in `svd-start`:
- **Bug.** bigstatsr computes the SVD's scaled products in blocks of
  `bigstatsr::block_size(n_indiv)` loci (128,684 for 1043 individuals; the
  whole dataset fits in one block below ~1.3e8 genotypes, which is why the
  smaller datasets were unaffected), calling the scaling function once per
  block with that block's `ind.col`. Our scaling closure in
  `radmix_svd_factors()` returned the centres and scales of all the loci
  whatever `ind.col`, so the loci of the second block were centred and
  scaled with those of the first (R warned "number of items to replace is
  not a multiple of replacement length"). Effect: a distorted SVD; the ALS
  never converged (1000 iterations, ~81 s each on 1 core), 66% of its
  allele frequencies sat at the bounds, its start was far worse than a
  random one (log-likelihood -2.4e8 against -1.02e8 at the optimum), and the
  fit needed 23 iterations, like random starts.
- **Fix** (folded into the `svd-start` commit, 7baa46c): the closure
  returns the values of the loci in `ind.col` (matched to the FBM columns
  of the loci of `x`), and a regression test makes bigstatsr use blocks of
  50 loci (`options(bigstatsr.block.sizeGB = ...)`) and checks that the
  factors are unchanged (it fails without the fix).
- **After the fix** (`laptop_results/svd_hgdp/`; K = 5, 1 core, `comb`, 4
  OpenBLAS threads available to the SVD):

  | start | SVD | start | fit | iterations | total | log-likelihood below the best |
  |---|---|---|---|---|---|---|
  | random, seeds 1-2 | - | - | 147 / 190 s | 20 / 26 | 147 / 190 s | 0.09 / 0.00 |
  | SVD, seeds 1-2 | 11.8 s (once) | 8.0 / 9.5 s | 42 s | 6 / 6 | 62 / 64 s | 0.10 / 0.10 |

  So at this size the SVD start is ~2.7x faster for one run (3.3x for
  further runs, which share the SVD), with the same optimum.
- **Also seen (before the fix):** the ALS was *slower* with 4 OpenBLAS
  threads than with 1 (126 s against 81 s for 1000 iterations, with 277 s
  of system time): its Armadillo products call R's BLAS, whose threads
  compete with each other (and, with `n_cores` > 1, with OpenMP's). With
  the fix the ALS converges in tens of iterations, so this matters less,
  but check it on a many-core machine (OpenBLAS uses all cores by default):
  time `radmix_svd_start()` with `RhpcBLASctl::blas_set_num_threads(1)` and
  with all threads.

## 3. What to verify, and what is new

- **V1. Correctness** on this machine: the full test suite against the
  installed builds `base`, `comb` and `comb_nomv` (and `comb_v4` if the CPU
  has AVX-512), and R CMD check of `kernel-speedups`.
- **V2. Each experiment alone** (div, dec, pairs, mv, hexp) against `base`,
  1 core: kernels and fits, as on the laptop. Do the gains (and the
  rejections) hold on a different CPU?
- **V3. Flags** (f_o3, f_avx2fma, f_native) against `base`, and whether
  multiversioning (`comb`) captures what `-march=native` gives (`comb_native`).
  If the CPU has AVX-512: is an x86-64-v4 clone (`comb_v4`) worth adding?
- **V4. The combination** against the baseline (`final_compare.sh`): 1 core,
  several cores, CV, rtol 1e-8. Also confirms 2.2 on this machine.
- **V5. Multi-core scaling (new)**: kernels and fits of `base` and `comb` on
  the HGDP (1043 x 130,943) and the unpruned HGDP (1043 x 642,794), from 1
  core to all physical cores (and SMT threads): speed-up, efficiency,
  thread invariance (identical results), OpenMP binding, NUMA. Where does
  scaling stop, and why?
- **V6. The SVD start at scale (new)**: random against SVD starts on the HGDP
  (1043 x 130,943), K = 3, 5, 7, 10, seeds 1-5: time (SVD, start, fit),
  iterations, quality (as in 2.1). Does the ~2.7x of the laptop's first
  full-HGDP run (2.5, K = 5 only) hold at other K and on more cores, and
  does the quality problem of HGDP quick K = 7 appear at full size? How do
  the SVD and the ALS scale with cores and BLAS threads (2.5)? Also at
  rtol 1e-8.
- **V7 (optional)**: ADAMIXTURE (`gt_adamixture()`) on the same HGDP and K
  as a reference for time and log-likelihood, if its environment can be
  installed (`tgc_tools_install()`, Linux only).

## 4. Setup

```bash
# the repository, on the coordinating branch
git clone git@github.com:EvolEcolGroup/tidygenclust.git   # or git fetch
cd tidygenclust
git switch radmix-perf-verify
git log --oneline -4    # the kit commit(s), 017e176, 7baa46c, ddcea53
```

R packages: tidygenclust's dependencies, plus **tidypopgen from its `dev`
branch** (the laptop had 0.4.4.9008; CRAN's 0.4.4 may lack what the SVD
start needs), testthat, withr, RhpcBLASctl (BLAS thread counts), and for the
data script bigsnpr:

```r
install.packages(c("pak", "RhpcBLASctl"))
pak::pak("EvolEcolGroup/tidypopgen@dev")
pak::local_install_deps(dependencies = TRUE)   # from the repository root
```

Then, from the repository root (set `TGC_WORK` to a disk with a few GB;
`TGC_TMP` too if /tmp is a small tmpfs; export them in every shell):

```bash
export TGC_WORK=$(dirname $PWD)/tgc_perf_work TGC_TMP=$(dirname $PWD)/tgc_perf_work/tmp
K=data-raw/radmix/perf_verification; S=$K/scripts; R=$K/desktop_results
mkdir -p $R
bash $S/machine_info.sh > $R/machine.txt 2>&1
Rscript $S/prepare_data.R          # downloads EUR_subset and the HGDP (~190 MB), ~10-20 min
bash $S/build_variants.sh          # the 11 default builds, ~1-2 min each
bash $S/build_variants.sh comb_native comb_o3   # and comb_v4 if avx512f is listed in machine.txt
```

Check: `prepare_data.R` should say "same file as on the laptop" for all four
datasets; `build_variants.sh` should show `-O2` and no `-march` for `base`,
`v3 0` clones for base/div/dec/pairs/hexp/f_*/comb_nomv/comb_native, 23
for comb and comb_o3, 29 for mv, and 23 + 23 (v3 + v4) for comb_v4 (as on
the laptop, where all 14 builds were made with these scripts).

## 5. Protocol

- Nothing else heavy running (check `uptime`). Do not change system settings
  (CPU governor, turbo, perf permissions) without asking the user; record
  them (`machine.txt`).
- Never time a `devtools::load_all()` build (it compiles at -O0). Only the
  builds of `build_variants.sh` (fresh sources, `--preclean`, private
  `R_MAKEVARS_USER`, so `~/.R/Makevars` is ignored).
- Alternate the builds (A B C A B C), with pauses (`PAUSE`), and use
  medians. On a desktop, pin the 1-core runs to one core (`PIN=2`) to cut
  migrations; check with a repeat that the noise is below the effects.
- Write everything to `data-raw/radmix/perf_verification/desktop_results/`
  (the `.rds` files are git-ignored), one subdirectory per phase, and
  commit the CSV, the summaries and the notes as you go.
- Compare with `laptop_results/` (there, "hgdp" is `hgdp_quick`).
- Physical cores vs SMT threads: report both, but judge scaling on physical
  cores.

## 6. Steps

Times are guesses for a fast desktop; run long phases in the background
(`nohup ... &` or the harness's background mode) and summarise as they end.

**A. Correctness (V1), ~1 h.** Against the installed builds:

```bash
bash $S/run_tests.sh base; bash $S/run_tests.sh comb; bash $S/run_tests.sh comb_nomv
```

Expect 0 failures (tests that need conda tools are skipped). Then R CMD
check of `kernel-speedups` from its exported sources:

```bash
cd $TGC_WORK && R CMD build --no-build-vignettes src/comb && \
  _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --ignore-vignettes tidygenclust_*.tar.gz; cd -
```

**B. Experiments alone (V2), ~1.5-2 h.**

```bash
bash $S/compare.sh $R/B_experiments 3 2 base div dec pairs mv hexp
```

`summary_<label>.txt`: kernel and fit ratios to `base`, and output
differences (dec must be bit-identical: 0).

**C. Flags (V3), ~1 h.**

```bash
bash $S/compare.sh $R/C_flags 3 2 base f_o3 f_avx2fma f_native comb comb_native comb_o3 [comb_v4]
```

**D. Combination (V4), ~1-1.5 h.** `NCORES` = e.g. "4 <physical cores>":

```bash
NCORES="4 16" bash $S/final_compare.sh $R/D_final
```

**E. Multi-core scaling (V5), ~2-3 h.** Kernels (one pass each, fixed Q
and P; `diff_vs_1core` must be 0), for `base` and `comb`, on both HGDPs, K = 5
and 10; cores 1, 2, 4, 8, ..., the number of physical cores, and the number
of logical CPUs:

```bash
C=1,2,4,8,16,32    # adapt to the machine
for v in base comb; do for d in hgdp hgdp_full; do
  Rscript $S/scaling_kernels.R $TGC_WORK/libs/$v $v $R/E_scaling $d 5,10 $C 3
done; done
# binding: the same for comb with threads pinned (and OMP_PROC_BIND=spread on NUMA)
OMP_PROC_BIND=close OMP_PLACES=cores Rscript $S/scaling_kernels.R $TGC_WORK/libs/comb comb $R/E_scaling hgdp 5,10 $C 3
```

Then whole fits (seed 1; thread invariance with `check_thread_invariance.R`):

```bash
for n in 1 4 16 32; do for v in base comb; do   # adapt the core counts
  Rscript $S/bench_fits.R $TGC_WORK/libs/$v $v $R/E_fits $n nocv 1e-7 hgdp:5,hgdp:10
done; done
for n in 8 32; do Rscript $S/bench_fits.R $TGC_WORK/libs/comb comb $R/E_fits $n nocv 1e-7 hgdp_full:5; done
Rscript $S/check_thread_invariance.R $R/E_fits comb
```

If scaling flattens, find out why before proposing anything (perf is
allowed if `perf_event_paranoid` <= 1; ask the user otherwise). Things to
know: Q-side work is parallelised over blocks of 64 individuals (1043
individuals = 17 blocks), split into segments of loci when there are fewer
than 32 blocks (`q_segments()`, `accum_q_segmented()` in
src/radmix_kernels.cpp), so the Q side has ~34 tasks here; the P side has
one task per 64 loci; the QP solves, the quasi-Newton step and R code are
serial or lightly parallel (Amdahl); the genotypes are 1 byte each (137 MB
per pass for the HGDP). A plot of speed-up against cores (kernels and
fits, base and comb) is welcome in the report.

**F. SVD start at scale (V6), ~2-3 h**, with the combination, on e.g. 8
cores (and once on 1 core for K = 5, to compare with the laptop's 1-core
figures):

```bash
Rscript $S/svd_vs_random.R $TGC_WORK/libs/comb $R/F_svd hgdp 3,5,7,10 1:5 8
Rscript $S/svd_vs_random.R $TGC_WORK/libs/comb $R/F_svd_1e8 hgdp 5,7 1:3 8 1e-8 random,svd
Rscript $S/svd_vs_random.R $TGC_WORK/libs/comb $R/F_svd_1core hgdp 5 1:2 1 1e-7 random,svd
```

The summary (`svd_vs_random_summary_hgdp.txt`) gives, per K and start:
time of the SVD, start and fit; iterations; runs at the best (< 0.1) and
near it (< 5); distinct modes; the largest Q difference from the best run.
Also time the SVD alone with 1 and with all BLAS threads
(`RhpcBLASctl::blas_set_num_threads()`; `radmix_svd_factors()` does not
pass `n_cores` to `gt_pca_partialSVD()`, so its threads come from the BLAS
only), and the start (`radmix_svd_start()`, whose ALS products also use the
BLAS) likewise (see 2.5); note the BLAS library. Compare
`F_svd_1core` with `laptop_results/svd_hgdp/` (same command). Optional (V7): `gt_adamixture()` on the
same data, K and seeds (`threads = 8`), for time and log-likelihood.

**G. Report.** `desktop_results/REPORT.md`: the machine (from
`machine.txt`); for each of V1-V7, a table with the laptop's figures next to
this machine's, and a short conclusion; the scaling plots; anything
unexpected (with numbers); and, for the user's decisions (2.4), the
evidence for and against. Keep the user's style: plain, numbers first, no
hype. Commit `desktop_results/` (not the `.rds`) to `radmix-perf-verify`
and push that branch.

## 7. Ground rules

- Do not change the commits of `svd-start` and `kernel-speedups`. If you
  find a bug, describe it in the report with a reproducible example, and
  put a proposed fix on a new branch from `kernel-speedups` (not pushed
  without asking).
- Commit and push only to `radmix-perf-verify`. No PRs, merges, or pushes to
  other branches without the user's go-ahead.
- New experiments (e.g. an AVX-512 clone, scaling fixes) are welcome as
  patches in `patches/` with their results, clearly marked as new.
- Keep the laptop's protocol, or say where and why you deviate.

## 8. Files of this kit

- `HANDOVER.md`: this file.
- `PROMPT.md`: the prompt the user gives to the Claude instance on the
  desktop.
- `notes/kernel_experiments_laptop.md`: the full log of the four kernel
  experiments (plan, protocol, per-experiment tables).
- `notes/speedup_plan.md`: the history: comparison with ADAMIXTURE and
  admixer, the perf profile, the SVD start exploration (experiment 5), the
  BLAS/Armadillo/Eigen measurements, the protocol and the correctness gates
  of the test suite.
- `patches/`: one patch per experiment, each applying to `svd-start`
  (7baa46c): `divisions.patch` (1a), `decode_v1.patch` (2, first version,
  no CV views), `decode.patch` (2), `multiversion.patch` (3, prototype,
  also on loglik and deviance), `pairs.patch` (4c + 4d),
  `hessian_expansion_REJECTED.patch` (4a).
- `laptop_results/`: CSV and summaries of the laptop runs. `e1` divisions
  (`div1` = one reciprocal, 1b; `div2` = two reciprocals, 1a); `e24`
  decoding (`dec`) and Hessian (`h_tp64` = 4a, `h_tp128`/`h_tp256` = 4a + 4b),
  divisions again (`div` = 1a); `e3` flags (`f_*`) and decoding v2 (`dec2`);
  `e4c` two loci per pass (`h_pairs0` = 4c, `h_pairs1` = 4c + 4d) and
  multiversioning (`mv`); `final` the combination (`comb`, `comb_nomv`) and
  rtol 1e-8; `vecinfo` GCC vectorisation reports (-O2 vs -O3); `svd` the SVD
  start exploration (`svd_cpp_results*.csv`: the C++ version, 5 seeds;
  others: the R prototypes); `svd_hgdp` the first run on the full HGDP,
  after the fix of 2.5 (`svd_vs_random.R ... hgdp 5 1:2 1 1e-7 random,svd`).
- `scripts/`:
  - `common.R` dataset paths and helpers;
  - `machine_info.sh` the machine;
  - `prepare_data.R` the datasets, with md5 checks;
  - `build_variants.sh` the builds;
  - `run_tests.sh` the test suite against an installed build;
  - `bench_kernels.R` and `bench_fits.R` the benchmarks;
  - `compare.sh` and `summarise.R` A/B comparisons;
  - `final_compare.sh` and `final_summary.R` the final comparison;
  - `scaling_kernels.R` and `check_thread_invariance.R` the scaling;
  - `svd_vs_random.R` the SVD start against random starts.
- `desktop_results/`: to be filled on the desktop.

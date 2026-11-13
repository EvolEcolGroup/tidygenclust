# Per-iteration speed-ups of gt_radmix(): laptop report (branch `newspeed`)

Night of 2026-11-13, about 03:18 to 04:57, on the laptop (i7-10510U). The
work was done by an agent, and this report was saved by the session that
coordinated it. The numbers also appear in `dev_notes/radmix_performance.md`
("The changes of 2026-11-13") and in `laptop_results/`. The commands for
the desktop are in `README.md`.

## Commits (on b5999cc)

- `378d830` The fused pass: the log-likelihood of a point is computed in
  the same pass as the next update of Q (`update_q_loglik()`).
- `e2e3661` Cross-validation: each fold fit builds a bit mask of its hidden
  genotypes once, instead of hashing every genotype in every pass.
- `d358c02` `slice_geno()` and `slice_coef()` are forced inline into the
  AVX2 clones. Before, the clones called the SSE2 baseline code.
- `373ae6c` A cheaper log-likelihood of a slice, inlined into the AVX2
  clone of the fused pass by default. `-DRADMIX_FUSED_LL_NOINLINE` turns
  the inlining off.
- `13f89be` `-DRADMIX_TILE_J`, `-DRADMIX_MIN_Q_TASKS` and
  `-DRADMIX_OMP_SCHEDULE` can be set at compile time, for benchmarks only.
  The defaults and the results are unchanged.
- `5a2b4d1` The desktop kit (this directory) and the update of
  `dev_notes/radmix_performance.md`.

## How the changes were checked

- **Tests:** the following all pass on the final code, and also passed on
  the build with only the fused pass:
  - `test_radmix_kernels` (108);
  - `test_radmix_cv` (47);
  - `test_gt_radmix` (106);
  - `test_radmix_vs_julia` (11), with identical iteration counts;
  - `test_radmix_vs_admixture` (42, with NOT_CRAN).
- **Builds:** all radmix objects are compiled at -O2. `objdump` shows that
  the AVX2 clones now call only `memset()` and `log()`.
- **Same iterations as before:** whole fits take the same iterations as on
  b5999cc.
  - HGDP pruned, K = 5, seeds 1-3: 21, 27 and 21 iterations, with no
    rejected steps. Along the whole history the log-likelihoods differ by
    at most 1.3e-14 (relative), and Q by at most 4e-13.
  - EUR, K = 3 and 4: 23 and 24 iterations, with 3 rejected steps each.
    The log-likelihoods differ by about 1e-16 (relative).
- **Thread invariance:** the histories and Q are identical on 1 and 2
  threads, and the `expect_identical()` checks of the tests pass.

## Measurements

Laptop, core 0 only, 1 thread. Other agents kept cores 1-3 busy (load
3-11), so all times are indicative only. The builds were alternated A B A B.

Passes over the genotypes (independent of the machine):

| | before | after |
|---|---|---|
| accepted iteration | 5 | 4 |
| rejected iteration | 6 | 5 |
| HGDP K = 5 fit, 21 iterations | 106 | 85 |
| HGDP K = 5 fit, 27 iterations | 136 | 109 |
| EUR K = 3 fit | 119 | 96 |
| fit of 8 iterations | 41 | 33 |

Single passes, HGDP K = 5, seconds (median of 2 runs, each a median of 3):

| build | loglik | update_q | update_p | update_q_loglik | accepted iteration |
|---|---|---|---|---|---|
| base | 2.13 | 1.84 | 2.14 | - | 10.09 |
| fused pass only | 2.20 | 2.08 | 2.30 | 3.31 | 9.99 (0.99) |
| final, loglik in the baseline code | 1.83 | 1.65 | 1.89 | 2.77 | 8.20 (0.81) |
| final (the default) | 1.75 | 1.68 | 1.94 | 2.27 | 7.83 (0.78) |

The fused pass on its own saved about 0.9 s of log-likelihood per
iteration. In that build, however, `update_q` and `update_p` measured 5-10%
slower in every run, although their code had not changed. This is probably
a code-layout effect, and it disappears in the final build.

Whole fits, final against base (including the 5 EM steps):

| data, K, seed | base | final | ratio |
|---|---|---|---|
| HGDP 5, seed 1 | 225.9 s | 175.6 s | 0.78 |
| HGDP 5, seed 2 | 277.8 s | 191.4 s | 0.69 |
| HGDP 5, seed 3 | 219.4 s | 160.1 s | 0.73 |
| EUR 3 (median of 2) | 26.3 s | 21.6 s | 0.82 |
| EUR 4 (median of 2) | 30.4 s | 26.8 s | 0.88 |

Cross-validation, HGDP K = 5, fold fits of 2 iterations from the same warm
start:

| view | time | relative to plain |
|---|---|---|
| plain | 25.0-27.0 s | 1 |
| fold hidden, hashed (as before) | 39.0 s | 1.44-1.56 |
| fold hidden, bit mask (building the mask included) | 25.6 s | 0.95 |

The mask takes 17 MB for the pruned HGDP and 87 MB for the unpruned HGDP
(computed, not measured).

The log-likelihood kernel ("candidate 2"): most of its cost was in code
outside the clones, not in `log()`.

- Its scalar fix-up loop ran on every slice, although only dosages, or
  likelihoods below 1e-30, need it.
- `slice_geno()` ran as baseline (SSE2) code in every pass.

With the log-likelihood inlined into the clone, the fused pass costs
`update_q` + 0.6 s (it was `update_q` + 1.25 s). Only 8 `log()` calls per
64 genotypes remain, so a vectorisable `log()` was not tried. One lesson is
recorded in the notes: in baseline code, `(double)(x > y)`,
`std::min`/`std::max` and `omp simd reduction(max:)` stop a loop from
vectorising.

Not measured on the laptop: 4 and 8 threads, HGDP K = 10, and the unpruned
HGDP. They are in the desktop kit.

## Decisions for the user

1. **The log-likelihood is inlined into the AVX2 clone by default**
   (373ae6c): 0.78 of base per iteration on the laptop, against 0.81 without
   the inlining. Confirm on the desktop by comparing `final` with
   `final_llnoinl`.
2. **The CV bit mask is used automatically** when the mask takes less than
   1 GiB (8.6e9 genotypes); above that, folds are hashed as before. The
   recommendation is to keep it automatic: the mask costs 1/8 of the size of
   the genotypes, and the results are bit-identical.
   `radmix_set_cv(..., mask = FALSE)` exists internally, for the tests.
3. **The results differ from b5999cc by rounding**, because qp is now
   accumulated with FMA in the clones. This is the same kind of change as the
   original multiversioning. The iterations are the same in every case
   tested, and the log-likelihoods agree to about 1e-14 (relative).
4. **Tile sizes, `MIN_Q_TASKS` and the OpenMP schedule** for many cores
   can only be judged on the desktop. Changing `MIN_Q_TASKS` changes the
   rounding.

## Open problems

- The 5-10% slowdown of `update_q` and `update_p` in the build with only
  the fused pass is unexplained. Check `fused` against `base` on the
  desktop.
- On the laptop, the inlining commit (d358c02) and the log-likelihood
  commit (373ae6c) were timed together, not separately.
- `bench_fits.R` calls the internal `radmix_fit()`. The commits of the
  algorithmic work (`data-raw/radmix/newspeed_algo/`) were later added to
  this branch; they leave `radmix_fit()` unchanged, and the kit builds the
  tip of the branch as `head`.

## Desktop tests

See `README.md` for the full commands. In short, with
`TGC_WORK=$HOME/tgc_newspeed_work` and `K=data-raw/radmix/newspeed`:

1. **Builds:** `MAKEFLAGS=-j8 PIN=0-7 bash $K/build_variants.sh
   base=b5999cc fused=378d830 cvmask=e2e3661 inline=d358c02 final=13f89be`,
   plus `head=newspeed` (the tip of the branch). Add the flag variants of `13f89be`
   with `EXTRA_FLAGS=...`: `final_llnoinl`, `final_q64`, `final_q128`,
   `final_guided` and `final_tj32`.
2. **Single passes:** `MODE=kernels CASES=hgdp:5,hgdp:10,eur:4 THREADS=1
   ROUNDS=3`, for base, fused, inline, final and final_llnoinl.
3. **Fixed 8 iterations:** `MAX_ITER=8 CASES=hgdp:5,hgdp:10,hgdp_full:5
   THREADS="1 4 8" ROUNDS=3`, for base, fused and final.
4. **Whole fits:** `CASES=hgdp:5,hgdp:10,hgdp_full:5,eur:3,eur:4
   THREADS="1 4 8" ROUNDS=3`, for base and final.
5. **Seeds:** `SEEDS=1,2,3 CASES=hgdp:5,hgdp:10,eur:3,eur:4 THREADS=8
   ROUNDS=1`, for base and final.
6. **Cross-validation:** `MODE=cv CASES=hgdp:5,hgdp_full:5 THREADS="1 8"
   ROUNDS=3`, for base, cvmask and final.
7. **Tiles and schedule:** `MAX_ITER=8 CASES=hgdp:5,hgdp:10,hgdp_full:5
   THREADS="1 8" ROUNDS=3`, for final, final_q64, final_q128, final_guided
   and final_tj32.

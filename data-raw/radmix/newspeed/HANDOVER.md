# Handover: the work of the night of 2026-11-13, to verify on the desktop

The night's work was done on the laptop, which has 4 cores, throttles,
and was shared by three agents. Its timings are therefore indicative only.
Everything that turns on a few percent, on scaling, or on more seeds is left
for the desktop (i9-13900, 8 P-cores). This file is the entry point. The
details are in the reports and READMEs listed below.

## Branches (local on the laptop, based on origin/dev b5999cc)

- **`newspeed`** has three parts:
  - the per-iteration speed-ups (`data-raw/radmix/newspeed/`: `REPORT.md`
    and `README.md`);
  - the algorithmic experiments and two opt-in options
    (`data-raw/radmix/newspeed_algo/`: `REPORT.md` and `README.md`);
  - this file.

  The commit `3cf356f` (`init = "subset"`) can be dropped on its own
  (`git rebase -i b5999cc`): it moves runs to lower modes.
- **`supervised`**: supervised analyses (`gt_radmix(supervised = )`, as
  ADMIXTURE's `--supervised`). See `data-raw/radmix/supervised/REPORT.md`.

The two branches change the same C++ files (`src/radmix_fit.cpp`,
`src/radmix_kernels.cpp`), so whichever is merged second needs its
conflicts resolved, and then the tests of both.

## What the laptop showed

- **Per iteration (newspeed):** an accepted iteration makes 4 passes over
  the genotypes instead of 5. Passes with a CV fold hidden cost the same as
  plain ones (they were 1.5x).
  - Whole fits on HGDP pruned K = 5 (1 core, seeds 1-3) took 0.69-0.78 of
    the time on `dev`; EUR K = 3 and 4 took 0.82-0.88.
  - The iterations were the same in every case, and the log-likelihoods
    agreed to about 1e-14 (relative).
- **CV refit policy (newspeed, option `cv_max_iter`/`cv_rtol`; the default
  is unchanged):**
  - On HGDP every policy tried is harmless.
  - On EUR the bias of short refits grows with K. rtol 1e-7 is the cautious
    middle: 30% of the CV passes saved, bias at most 0.66e-4.
  - This is a decision for you; it needs more seeds on the desktop.
- **Subset start (newspeed, `init = "subset"`):** about half the passes,
  but it is no better on the unpruned than on the pruned HGDP, and it moved
  runs to lower modes (HGDP K = 10, EUR K = 3). Probably drop it.
- **Supervised (supervised):**
  - It agrees with ADMIXTURE `--supervised` on EUR and on HGDP: Q of the
    individuals to estimate within 6e-5 on HGDP.
  - The unsupervised path is bit-identical to `dev`.
  - Starting P from the reference frequencies saved about 45% of the passes
    on HGDP, but nothing on EUR. It is documented only, not built in.

## Runs on the desktop, in order of priority

Set up once:

- Clone or fetch the repository with both branches (they must be pushed
  from the laptop first: `git push -u origin newspeed supervised`).
- Copy `vignettes/articles/data/hgdp/` (`hgdp_bench.*`, `hgdp_full.*`,
  `hgdp_bench_quick.*`, `hgdp650_id_pop.txt`; git-ignored) and the
  EUR_subset files (`~/.cache/R/tidygenclust/EUR_subset.*`).
- Run everything on an otherwise idle machine. Use installed `-O2` builds
  only, with threads pinned to the P-cores. The scripts do both.

Then:

1. **Per-iteration kit** (`data-raw/radmix/newspeed/README.md`,
   "Commands"):
   - the builds (base, fused, cvmask, inline, final, head, and the flag
     variants of final);
   - then `kernels`, `iter8`, `whole`, `seeds` and `cv`.

   This confirms the 0.7-0.8 on 1, 4 and 8 threads, on HGDP K = 5 and 10
   and on the unpruned HGDP. It also decides whether the log-likelihood is
   inlined into the clone (`final` against `final_llnoinl`), and explains
   or dismisses the 5-10% slowdown of `fused` alone.
2. **Tiles and schedule** (same kit, step `tiles`): only if step 1 leaves
   time. Note that `MIN_Q_TASKS` changes the rounding.
3. **CV policy** (`data-raw/radmix/newspeed_algo/README.md`): `cv_timing`,
   then `cv` (HGDP K = 2-10 with seeds 1-3, EUR K = 1-6 with seeds 1-5).
   This step is long; it gives the evidence for the CV default.
4. **Warm-ups** (same README): `warmup hgdp`, then `warmup eur`. These
   settle the question of the modes of the subset start, with seeds 1-5 at
   K = 10. Run `warmup hgdp_full` only if the subset start is still a
   candidate.
5. **Supervised** (`data-raw/radmix/supervised/REPORT.md`, "Desktop
   tests"): the HGDP timings at 1 and 8 threads, thread invariance at 8
   threads, ADMIXTURE `--supervised -j8` against `gt_radmix(n_cores = 8)`,
   and the full test suite with `NOT_CRAN=true`.

## Decisions waiting for the user

- newspeed:
  - the log-likelihood inlined into the AVX2 clone (the default now);
  - the CV bit mask, used automatically when it takes less than 1 GiB;
  - results that differ from `dev` by rounding (FMA in the clones);
  - the CV refit default;
  - keeping or dropping `init = "subset"`.
- supervised:
  - the API: a vector with `NA` for the individuals to estimate; the column
    order (first appearance, or the factor levels); and the `ref_pops`
    element;
  - K equal to the number of populations (as in ADMIXTURE; fastmixture
    allows extra clusters);
  - `init = "svd"` is rejected;
  - whether to build in the start from the reference frequencies;
  - CV that scores the genotypes of the reference individuals.

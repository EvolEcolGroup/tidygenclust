# Supervised analyses in `gt_radmix()`: report

Branch `supervised` (from `origin/dev`, b5999cc), night session of
2026-11-13, about 2 hours, unattended. Laptop (i7-10510U, 4 cores, no SMT),
all four cores busy with three agents (load 5-11); this work used core 1
only, so every time below is indicative, and the comparisons rest on
iteration counts, log-likelihoods and per-kernel ratios.

## Commits

- `9691b83` Supervised analyses in gt_radmix() (ADMIXTURE's --supervised):
  API, validation, the mask in the fit and the compact Q side, CV, docs,
  NEWS, dev notes, tests.
- `16df1ec` Read scattered free individuals of a supervised fit through
  their rows (fast path for non-consecutive free individuals; test of the
  general path with pseudohaploids).
- `bf384a7` Test supervised gt_radmix() against ADMIXTURE --supervised
  (fixture `admixture_eur_supervised.rds`, 4 KB; fixture script; README).
- `5b1c8fd` Document the start of supervised fits from reference
  frequencies (roxygen).
- the last commit: these scripts, this report, `Notes_to_do`.

## What was built

- **API** (decision for the user, see below): `gt_radmix(x, k,
  supervised = <vector>)`, one element per individual of `x`, character or
  factor: the population of each reference individual, `NA` for the
  individuals to estimate. `k` must be a single value equal to the number
  of populations (ADMIXTURE's rule). The columns of Q and P follow the order
  of first appearance of the labels (ADMIXTURE's order), or the levels of a
  factor. The result gets a new element `ref_pops` (the population of each
  column); the columns keep the names `.Q1`, `.Q2`, ..., which tidypopgen's
  `tidy()`, `augment()` and `autoplot()` select (tested).
- **Validation** (`radmix_check_supervised()` in `R/gt_radmix.R`): type,
  length, all `NA`, factor levels without a reference individual, number of
  populations different from `k` (the message lists them), several `k`;
  `fix_q`, `fix_p` and `init = "svd"` are rejected with `supervised`;
  `init_q` (rows of the references ignored) and `init_p` are accepted;
  `crossval = TRUE` works.
- **Fit** (`src/radmix_fit.cpp`): an I x K matrix of fixed values (NA for
  the free individuals) becomes a mask. The references are set to their
  fixed values right after the random draw (before EM; the values drawn for
  them are discarded, so the draws of the other individuals do not move),
  EM and the Q block updates copy their columns, they are restored exactly
  after each quasi-Newton step (as a fixed Q block already was), and P uses
  every individual. The log starts with "Supervised analysis: N reference
  individuals (fixed ancestry), M individuals to estimate". Fixed values:
  the indicator projected by `radmix_project_cpp()` onto
  {sum = 1, q >= 1e-5}, i.e. 1 - (K - 1) 1e-5 and 1e-5, which is exactly
  what ADMIXTURE reports.
- **Skipping the references** (`src/radmix_kernels.cpp`): the Q-side
  passes (`em()`, `update_q()`) build the list of the free individuals and
  work on a compact transpose of their columns of Q, in blocks of 64
  consecutive entries of the list, so neither the sums nor the quadratic
  programs of the references are computed, whatever their order. A slice of
  free individuals that are consecutive rows of the backing file keeps the
  fast path; scattered ones are read from the column of the locus through
  their rows (commit 2; before it they took the per-genotype general path,
  which made a mask with few references 9% slower than no mask on HGDP);
  pseudohaploids keep the general path. The sums of a free individual are
  the same terms in the same order as without the mask (tests compare at
  1e-12). A first version that only trimmed each block of 64 to the range
  of its free individuals was dropped before committing: HGDP is not sorted
  by population (609 runs of population labels in 1043 individuals), so it
  would have skipped almost nothing there.
- **Cross-validation**: ADMIXTURE accepts `--cv` with `--supervised`, so
  `crossval = TRUE` is allowed; the references stay fixed in the fit of
  each fold, and their genotypes belong to the folds and are predicted like
  the others (ADMIXTURE does not document this detail).
- **Unsupervised path unchanged**: `unsupervised_identical.R` run with the
  baseline build (`libs/base`) and this branch (commit 1): the whole
  gt_admix objects (Q, P, log, log-likelihood, CV) are `identical()` for
  EUR K = 3 (1 and 2 threads), EUR K = 2 with CV, pop_a K = 2:3 x 2 runs,
  pop_a with CV, pseudohaploids, EUR imputed with `init = "svd"`, `fix_p`
  and `init_q`. Commit 2 only touches the supervised branch of
  `slice_geno()`.
- **Tests**: `test_gt_radmix.R` (references kept exactly at their fixed
  values, column order and `ref_pops`, recovery of the simulated ancestry
  of the others without relabelling, factor levels, thread invariance,
  tidypopgen methods; CV with the references fixed in a fold fit;
  `init_q`; every error), `test_radmix_kernels.R` (masked `update_q()` and
  `em()` equal the unmasked ones for free individuals and copy the fixed
  ones, with runs and scattered free individuals, segments, a CV view,
  pseudohaploids, 1 vs 3 threads, all fixed), and
  `test_radmix_vs_admixture.R` (supervised fixture, see below). All pass
  (radmix test files, against an installed -O2 build; `NOT_CRAN=true` for
  the ADMIXTURE file). The new tests fail on the baseline (unknown argument).
- **Docs**: roxygen (`@param supervised`, "Supervised analyses" section,
  `ref_pops` in the value, an example, the "Comparison with ADMIXTURE"
  bullet), NEWS.md, `dev_notes/admixture_implementation.md` (new section 9),
  `data-raw/radmix/README.md` (the new fixture), `Notes_to_do`.

## ADMIXTURE's behaviour (established with the binary)

EUR_subset, K = 3, seed 1, `-j1`; references = individuals with q > 0.9 in a
cluster of the stored unsupervised ADMIXTURE K = 3 fit (15, 59, 28; 277 to
estimate), labelled mid/zeta/alpha for clusters 1-3 so that first
appearance (zeta, alpha, mid) differs from alphabetical order.

- **Q of the references**: exactly 0.999980 and 0.000010 for all 102
  references, i.e. the indicator projected onto {sum = 1, q >= 1e-5}; not
  1/0. Reproduced (identical to ADMIXTURE's 6 decimals).
- **Column order**: order of first appearance in the `.pop` file (columns
  zeta, alpha, mid). Reproduced for character labels; factors use their
  levels (a choice, see decisions).
- **K different from the number of populations** (K = 2 and 4 with 3
  populations): "Expected 4 populations, but found reference individuals
  for 3." then "Error reading .pop file!", and exit status 0. We error.
- **`--cv` works with `--supervised`**: "CV error (K=3): 0.43580", same
  estimates as without `--cv`.
- **Log and convergence**: "Supervised analysis mode. Examining .pop
  file...", five EM steps, then QN/Block iterations; 18 iterations to
  delta < 1e-4, log-likelihood -15600982.798954; ends with the Fst
  divergences between the estimated populations (as unsupervised).
- **Time**: 70.2 s (1 thread, laptop at load ~3).

## Validation against ADMIXTURE

**EUR_subset, K = 3** (`compare_eur.R`; the same comparison is a test,
with the fixture `admixture_eur_supervised.rds`):

| | ADMIXTURE | gt_radmix |
|---|---|---|
| log-likelihood | -15600982.7990 | -15600982.8121 |
| iterations | 18 | 14 |
| time, 1 thread (indicative) | 70.2 s | 14.6 s |

- Q of the 277 individuals to estimate: mean |difference| 3.4e-4, max
  2.3e-3 (columns compared without any matching: same order).
- Q of the 102 references: identical to ADMIXTURE's 6 decimals.
- The log-likelihood difference (0.013, 8e-10 relative) is the
  convergence criterion: ADMIXTURE stops at an absolute change of 1e-4,
  gt_radmix at a relative change of 1e-8 (0.16 units here).
- Cross-validation (5 folds, different from ADMIXTURE's): 0.43547 against
  ADMIXTURE's 0.43580 (0.08%); the test checks it with a relative
  tolerance of 1e-3. Q and P are identical with 1 and 2 threads.

**HGDP pruned, K = 5, "regions" references** (660 references, 383 to
estimate; `hgdp_supervised.R`, ADMIXTURE `--supervised -s 1 ... 5 -j1`,
`compare_hgdp_admixture.R`):

| | ADMIXTURE | gt_radmix |
|---|---|---|
| log-likelihood | -101972841.650 | -101972837.992 |
| iterations | 12 | 11 |
| time, 1 thread (indicative) | 1543 s (about 50 s per iteration once alone on the core) | 121 s (about 11 s per iteration with the EM steps) |

- Q of the 383 individuals to estimate: mean |difference| 7.5e-6, max
  5.8e-5; references identical to 6 decimals; P mean |difference| 1.0e-6,
  max 6.4e-4.
- gt_radmix's log-likelihood is 3.66 units (3.6e-8 relative) higher. Our
  log-likelihood of ADMIXTURE's written estimates is -101972841.643, so
  ADMIXTURE's point really is slightly lower: it stopped a little short
  (its last steps changed the log-likelihood by less than 1e-4).
- ADMIXTURE's time includes its first iterations on a core shared with my
  other jobs (105-240 s each); alone it took 45-55 s per iteration.

**HGDP pruned, K = 5, "pops" references** (139 references, 904 to
estimate; same scripts):

| | ADMIXTURE | gt_radmix |
|---|---|---|
| log-likelihood | -101817791.999 | -101817790.851 |
| iterations | 15 | 14 |
| time, 1 thread (indicative) | 1312 s (about 65 s per iteration) | 168-184 s |

- Q of the 904 individuals to estimate: mean |difference| 3.5e-6, max
  4.3e-5; references identical to 6 decimals; gt_radmix's log-likelihood
  1.15 units higher.

## Efficiency (HGDP pruned, 1043 x 130,943, K = 5, 1 core)

References: "regions" = every individual of five Pop7Groups regions
(Africa, Europe, Est_Asia, America, Oceania: 660 references, 383 to
estimate); "pops" = every individual of Yoruba, French, Han, Pima and Papuan
(139 references, 904 to estimate). HGDP's individuals are not sorted by
population, so the free individuals of a block are scattered in both cases.

**Whole fits** (`run_hgdp_timings.sh`, commit 1 build, i.e. before
scattered free individuals got their fast path, which mostly matters for
"pops"; alternated
U R P U R P, seed 1, core 1 while the three other cores ran other agents'
jobs; load in brackets):

| case | iterations | log-likelihood | time (s), rep 1 / rep 2 |
|---|---|---|---|
| unsupervised | 21 | -101809867.106 | 255.2 (8.3) / 262.5 (6.3) |
| supervised, regions | 11 | -101972837.992 | 121.5 (5.5) / 121.0 (5.9) |
| supervised, pops | 14 | -101817790.851 | 184.0 (5.2) / 168.0 (5.9) |

Repeats agree within 3% (regions, unsupervised) to 9% (pops). Most of the
gain of the supervised fits comes from needing fewer iterations (the
references pin the clusters down); the rest from skipping the references.

**Per pass** (`time_kernels.R`, median of 3 calls, 1 thread, the same
random Q and P; seconds):

| mask | free individuals | `update_q()` | EM step |
|---|---|---|---|
| none (unsupervised) | 1043 | 1.85 | 2.79 |
| regions | 383 (37%) | 0.72 (-61%) | 1.88 (-33%) |
| pops | 904 (87%) | 1.74 (-6%) | 2.66 (-5%) |

(load 4.1-4.5, commit 2 build). With the commit 1 build, before scattered
free individuals were read through their rows, the same calls took 1.92
(none), 0.84 (regions) and 2.08 s (pops: 9% slower than no mask), EM 2.78,
2.01 and 3.08 s.

`update_p()` (2.14 s) and the log-likelihood (2.03 s) are not affected by
the references. An iteration is two Q updates, two P updates and a log-likelihood:
about 9.9 s of passes unsupervised, 7.7 s for "regions" (-23%) and 9.7 s
for "pops" (-2%). The Q side scales with the number of free individuals
(0.72 s for 37% of them, i.e. about 1.05 times the cost per individual of
the unmasked pass); the P side and the log-likelihood, about 60% of an
unsupervised iteration, still read every genotype, which bounds the gain
per iteration.

## Starts

`compare_starts.R`: "random" = random start + 5 EM steps (ADMIXTURE's, the
default); "freq_em" = `init_p` from the allele frequencies of the
reference populations, Q of the others from the 5 EM steps; "freq" = the
same P, Q of the others random, no EM (fastmixture's start). Iterations and
final log-likelihood (times were taken on a shared core and are not
comparable):

| data | seed | random | freq_em | freq |
|---|---|---|---|---|
| EUR, K = 3 (102 refs of 379) | 1 | 14, -15600982.812 | 15, -15600982.827 | 14, -15600982.883 |
| EUR, K = 3 | 2 | 13, -15600982.835 | 15, -15600982.853 | (stopped) |
| HGDP regions, K = 5 (660 of 1043) | 1 | 11, -101972837.992 | 8, -101972837.990 | 7, -101972838.017 |

With many references, fastmixture's start saves the 5 EM steps (10 passes)
and about 4 iterations (20 passes): about 36 passes instead of 66, i.e.
about 45% of a fit, for the cost of computing the reference frequencies
(about one pass over the references). With fewer references it gains
nothing. All starts reach the same optimum (differences of 0.02-0.07
log-likelihood units, below the convergence tolerance of ~1 unit). It is
documented (commit 4) as `init_p` + `n_em_iter = 0`, not built in: see
decisions.

## Decisions for the user

1. **API**: `supervised =` a vector with one element per individual, `NA`
   for the individuals to estimate (my brief's suggestion, kept).
   Alternatives: also accept ADMIXTURE's `"-"` as "estimate" (so that
   `readLines("x.pop")` works directly; now `"-"` would be a population and
   the `k` error would list it), or a path to a `.pop` file.
2. **Column order**: first appearance for character labels (ADMIXTURE),
   the levels for a factor (R's way to choose an order; a level without
   references is an error suggesting `droplevels()`). Alternative: first
   appearance for both, or alphabetical.
3. **Recording the populations**: new element `ref_pops`; the Q and P
   columns keep `.Q1`/`.P1`... (tidypopgen's methods select `.Q*`).
   tidypopgen's `c()` for gt_admix drops unknown elements, so
   `c(sup_fit, other_fit)` loses `ref_pops` (as it loses `converged`).
   Alternative: name the P columns by population (no tidypopgen
   constraint there).
4. **K must equal the number of reference populations** (ADMIXTURE).
   fastmixture also allows K larger than the number of populations (extra
   clusters without references). The C++ would support it as is (a cluster
   without references is just estimated); only the R check forbids it.
5. **`init = "svd"` is rejected** with `supervised` (its clusters come in
   an arbitrary order). It could be supported by matching the ALS clusters
   to the populations from the references' Q, but the frequency start
   (item 6) is simpler and better suited.
6. **Start**: ADMIXTURE's random start + 5 EM steps is the default; the
   frequency start (fastmixture's) is documented as `init_p` +
   `n_em_iter = 0`. Building it in (e.g. `init = "freq"`, computing the
   reference frequencies in C++ in one pass) would save ~45% of a fit when
   the references are many (HGDP regions) and nothing with few (EUR).
7. **Cross-validation**: the references' genotypes are in the folds and
   scored like the others. Alternative: score only the individuals to
   estimate (ADMIXTURE's choice is unknown; its CV error on EUR, 0.43580,
   matches ours within 1e-3 with the references included).
8. `n_runs > 1` is allowed with `supervised` (runs differ by their random
   start of the free individuals and of P); nothing special was done.

## ADMIXTURE behaviour not reproduced (or not checkable)

- The Fst divergences between the estimated populations printed at the
  end of the log (not computed, as in the unsupervised mode).
- The `.pop` file and its `"-"` marker (a vector with `NA` instead).
- Exit status 0 with only a message when K differs from the number of
  populations: we stop with an error.
- ADMIXTURE's random folds and its convergence criterion (absolute
  delta < 1e-4; ours relative 1e-8), as in the unsupervised mode.
- Not checkable without the source: whether ADMIXTURE fixes the references
  before the EM steps (we do), and whether its CV folds include the
  references' genotypes (we include them).

## Open problems

- The P side and the log-likelihood still read the references' genotypes
  (they must, for P), so the gain per iteration is bounded at about 40%
  even with almost all individuals as references; the large gains come
  from the fewer iterations of supervised fits, and from the frequency
  start when the references are many.
- Timings are from a laptop with every core busy (throttling); only the
  ratios between alternated runs and the per-kernel ratios mean much.
- The multi-thread behaviour of the compact Q side (fewer blocks of 64
  when few individuals are free: 383 free individuals give 6 blocks,
  hence segments of loci) was only checked for identical results (1 vs 2
  or 3 threads), not timed.
- `tidypopgen`'s `c()` drops `ref_pops` (see decisions).

## Tests worth running on the desktop

From the root of a checkout of the branch, with the branch installed at -O2
in `<lib>` and a scratch dir `<s>` (the HGDP gen_tibble is made there by the
first call of `hgdp_supervised.R`):

1. Whole fits, 1 core, alternated, 3 repeats:
   `bash data-raw/radmix/supervised/run_hgdp_timings.sh <lib> <s> <core> 3`
   (cases unsupervised, regions, pops; results in `<s>/hgdp_supervised.csv`).
   The same on 8 cores: edit the script's `n_cores` argument (the fifth
   argument of `hgdp_supervised.R`) and `taskset` list.
2. Per-pass costs: `Rscript data-raw/radmix/supervised/time_kernels.R <lib> <s> 5`.
3. Starts on HGDP with more seeds:
   `Rscript data-raw/radmix/supervised/compare_starts.R <lib> <s> hgdp 1 1:3`
   (and `... eur 1 1:3`).
4. ADMIXTURE against gt_radmix on 8 cores, timed (the comparisons above
   are 1 thread on a shared laptop): in `<s>/admix` (the `.pop` files are
   written by `hgdp_supervised.R`; link `hgdp_<case>.bed/.bim/.fam` to
   `hgdp_bench.*`), `admixture --supervised -s 1 hgdp_<case>.bed 5 -j8`
   for `<case>` = regions and pops, then
   `Rscript data-raw/radmix/supervised/compare_hgdp_admixture.R <s> <case>
   <s>/hgdp_<case>_random_1_<lib basename>.rds`, against
   `hgdp_supervised.R` with `n_cores` 8.
5. Thread invariance on HGDP: `gt_radmix(hgdp, k = 5, supervised = labels,
   seed = 1, n_cores = 1)` and `n_cores = 8` give `identical()` Q and P.
6. The full test suite (`devtools::test()` with `NOT_CRAN=true`, or
   `R CMD check`), which I ran only file by file (`test_gt_radmix.R`,
   `test_radmix_kernels.R`, `test_radmix_cv.R`,
   `test_radmix_vs_admixture.R`).

## Scripts in this directory

- `compare_eur.R`: supervised `gt_radmix()` against ADMIXTURE
  `--supervised` on EUR_subset, K = 3 (runs ADMIXTURE if needed).
- `unsupervised_identical.R`: unsupervised fits of one build, saved for an
  `identical()` comparison between builds.
- `hgdp_supervised.R` and `run_hgdp_timings.sh`: HGDP pruned K = 5,
  supervised (two sets of references) and unsupervised; writes the `.pop`
  files for ADMIXTURE.
- `compare_hgdp_admixture.R`: a saved HGDP fit against ADMIXTURE's output.
- `time_kernels.R`: per-pass costs with and without the masks.
- `compare_starts.R`: random start + EM against starts from the reference
  frequencies.
- `results/`: the CSV outputs of the HGDP runs quoted above
  (`hgdp_supervised.csv`, `time_kernels.csv` of the commit 2 build,
  `compare_starts_hgdp.csv`); the load average of each run is recorded.

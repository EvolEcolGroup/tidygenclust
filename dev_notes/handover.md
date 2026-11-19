# tidygenclust: handover

Updated 2026-11-19 (manica01): version 0.1.5.9006 on `dev` (d889f63), the
first with #38's cross-validation and SVD starts for likelihoods, which
tidypopgl now requires (its PR #47); `radmix_extensions` deleted on origin
and manica01. `summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| laptopX1 local: `radmix_extensions` | squash-merged as #38 (df747aa); deleted on origin and manica01 (19 Nov) | | | yes: `git branch -D` |
| `articles_revised` | revision of the pseudohaploid article (2 commits) | `dev` (old) | yes | no: conflict in `pseudohaploids.Rmd` |
| `hgdp_plots`, `windows` (origin) | not part of the radmix work; commits not in `dev` | | yes | no: ask Andrea before touching |

Deleted 18 Nov (laptopX1): 17 local branches without origin and the
`coverage-speed` worktree; tips in
`/media/andrea/annex/tmp/claude/tidygenclust/deleted_branches_2026-11-18.txt`.
Deleted 17 Nov (origin and manica01): `newspeed`, `supervised` (tags
`archive/newspeed` ec96315, `archive/supervised` fc43f9b),
`newspeed_supervised`, `coverage-speed`, `newspeed_port` (#36), and the
local `adamixture`, `rewrite`, `openadmixture`, `rclump`.

## Next tasks

### 1. **P2** [review] tidypopgl's PR #47 (`radmix_source_speed` into `gl_workflow`)

- tidypopgl's `gl_workflow` works with tidygenclust `dev` except one test:
  its CI (PR #46) fails since #38, as `test_gt_radmix_data.R` still expects
  `gt_radmix(crossval = TRUE)` to stop on likelihoods.
- PR EvolEcolGroup/tidypopgl#47 fixes it: the faster source, `Imports:
  tidygenclust (>= 0.1.5.9006)`, a test that cross-validation through the
  source gives the genotypes' error, and the vignette choosing `k` by
  cross-validation. Full suite passed on manica01 (19 Nov); CI on the PR
  started on all three platforms (its AVX2 file had only been built on
  the desktop).
- **Done when:** CI passes and Andrea merges #47 into `gl_workflow`
  (tidypopgl's handover has the rest).

### 2. **P2** [claude] then [review] `articles_revised`

- **Branch:** `articles_revised`; rebase or merge `dev` into it and resolve
  the conflict in `vignettes/articles/pseudohaploids.Rmd`.
- **Check:** render the article (`TMPDIR` in a work folder).
- **Done when:** a PR into `dev` is open for Andrea's review.

### 3. **P3** [decide] The bounds of the estimates with genotype likelihoods

- Keep `[1e-5, 1 - 1e-5]` (ADMIXTURE's) or use NGSadmix's `1e-9`? With
  1e-5 the log-likelihood is about 1e-6 (relative) below NGSadmix's on data
  with many estimates at the bounds (`data-raw/radmix/gl/REPORT.md`).
  Recommendation: keep 1e-5 (one set of bounds for both kinds of data).

### 4. **P3** [cleanup] Branches and scratch files

- laptopX1: the local `radmix_extensions` (`git branch -D`; merged as
  #38, gone from origin and manica01), and the
  untracked `data-raw/radmix/perf_verification/` (a copy of what the tag
  `archive/radmix-perf-verify` holds).
- manica01: delete `~/tgc_newspeed_work` (2.2 GB of builds and cached
  data), `~/tgc_bench_article` (logs of the 16 Nov run; summaries are in
  `dev_notes/radmix_port_timing_2026-11/`), the untracked
  `vignettes/articles/data/hgdp/runs` (4 GB of raw outputs of the article
  run; the article's data are committed), and `~/tgc_supfused_work` (the
  builds and simulated data of 18 Nov; the scripts and results are in
  `dev_notes/radmix_extensions_2026-11/`). `tgc-desktop-checks.timer` is
  disabled (17 Nov).

## Desktop queue

Empty. The driver of the last unattended run (port timing and the full
article benchmark, 16 Nov) is `dev_notes/radmix_port_timing_2026-11/run_desktop.sh`.

## Backlog

- Under covr, `-fprofile-update=single` might make the multithreaded tests
  fast enough not to skip (`covr_timing_2026-11.md`, "Open questions").
- Make PR runs fail on Codecov upload errors:
  `fail_ci_if_error: ${{ github.event_name != 'pull_request' ||
  secrets.CODECOV_TOKEN != '' }}` (the current expression gives the token,
  not `true`, so PRs hid the deactivated repository).

- The vignette's comparison of P as code (tidypopgen PR #211 is merged).
- P as a matrix (`gt_radmix()`) or a data.frame (`gt_admixture()`):
  decide.
- Bootstrap standard errors, as ADMIXTURE's `-B`.
- Speed leads in `radmix_performance.md` ("Leads not pursued yet"):
  `-march=native` gap, `-O3`, admixer's damped retries, single-precision
  Hessians for K >= 8.
- Genotype likelihoods: an EM-based PCA (PCAngsd's iterations) for the
  SVD start; `loglik()` of
  likelihoods runs in baseline code (only used after a rejected
  quasi-Newton step, by verbose EM priming and by the SVD start).
- Move the notes still in `data-raw/radmix/` (`Notes_to_do`, `gl/REPORT.md`)
  into `dev_notes/`.

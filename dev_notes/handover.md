# tidygenclust: handover

Updated 2026-11-18 (manica01, overnight): branch `radmix_extensions` (fused
pass for supervised fits, cross-validation and SVD starts for genotype
likelihoods, inlined likelihood helpers), pushed for review.
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `radmix_extensions` | fused pass for supervised fits; cross-validation and SVD starts for genotype likelihoods; likelihood helpers inlined into the AVX2 clones (7 commits + notes) | `dev` 00deb4f | yes | no: Andrea's review |
| `articles_revised` | revision of the pseudohaploid article (2 commits) | `dev` (old) | yes | no: conflict in `pseudohaploids.Rmd` |
| `hgdp_plots`, `windows` (origin) | not part of the radmix work; commits not in `dev` | | yes | no: ask Andrea before touching |
| laptopX1 local: `coverage-speed` (worktree `/media/andrea/annex/tmp/claude/tidygenclust/wt-coverage`), `newspeed`, `supervised`, `newspeed_supervised`, `newspeed_port`, `gl_gl16`, `gl_private`, `radmix-speedups`, `radmix-perf-verify`, `kernel-speedups`, `svd-start` | merged or superseded; origin copies deleted | | | to delete |

Deleted 17 Nov (origin and manica01): `newspeed`, `supervised` (tags
`archive/newspeed` ec96315, `archive/supervised` fc43f9b),
`newspeed_supervised`, `coverage-speed`, `newspeed_port` (#36), and the
local `adamixture`, `rewrite`, `openadmixture`, `rclump`.

## Next tasks

### 1. **P2** [review] `radmix_extensions`, then tidypopgl's `radmix_source_speed`

- **Branch:** `radmix_extensions` (from `dev` 00deb4f, pushed; PR #38). Commits: ec7ff9a fused pass for supervised fits; ecc0841
  cross-validation for likelihoods; 4ed9ceb SVD start for likelihoods;
  7561833 NEWS; 9b05912 tests made obsolete by the two; 9659e51 inlined
  likelihood helpers; then the notes; 014eb30 `svd_power` argument
  (18 Nov, laptopX1, at Andrea's request). PR #38 into `dev`.
- **Read first:** `admixture_implementation.md` ("The fused pass of
  supervised fits", "Cross-validation", "The SVD start" under "Genotype
  likelihoods"), `radmix_performance.md` (the two sections of 18 Nov),
  `radmix_extensions_2026-11/` (scripts and results).
- **Checked:** full suite with `NOT_CRAN=true` on an `-O2` build (22 files,
  3874 expectations, after the fixes of 9b05912), `R CMD check --as-cran`
  (see `progress.md`); tidypopgl's suite passes with this branch.
- **Questions for the review:** the CV error of likelihoods is not
  comparable with that of genotypes (documented). The power iterations of
  the likelihood SVD (3 by default) are now set with `svd_power`.
- **Then:** merge tidypopgl's `radmix_source_speed` (independent; see its
  handover), and update tidypopgl's vignette sentence that cross-validation
  is not available for likelihoods.

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

- laptopX1: the local branches in "Branches" (`git branch -D`), the
  worktree (`git worktree remove
  /media/andrea/annex/tmp/claude/tidygenclust/wt-coverage`), and the
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

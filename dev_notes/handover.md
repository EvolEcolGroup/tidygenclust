# tidygenclust: handover

Updated 2026-11-17 (manica01): PR #36 merged into `dev` (47f5780) and the superseded branches deleted.
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `articles_revised` | revision of the pseudohaploid article (2 commits) | `dev` (old) | yes | no: conflict in `pseudohaploids.Rmd` |
| `hgdp_plots`, `windows` (origin) | not part of the radmix work; commits not in `dev` | | yes | no: ask Andrea before touching |
| laptopX1 local: `coverage-speed` (worktree `/media/andrea/annex/tmp/claude/tidygenclust/wt-coverage`), `newspeed`, `supervised`, `newspeed_supervised`, `newspeed_port`, `gl_gl16`, `gl_private`, `radmix-speedups`, `radmix-perf-verify`, `kernel-speedups`, `svd-start` | merged or superseded; origin copies deleted | | | to delete |

Deleted 17 Nov (origin and manica01): `newspeed`, `supervised` (tags
`archive/newspeed` ec96315, `archive/supervised` fc43f9b),
`newspeed_supervised`, `coverage-speed`, `newspeed_port` (#36), and the
local `adamixture`, `rewrite`, `openadmixture`, `rclump`.

## Next tasks

### 1. **P2** [claude] then [review] `articles_revised`

- **Branch:** `articles_revised`; rebase or merge `dev` into it and resolve
  the conflict in `vignettes/articles/pseudohaploids.Rmd`.
- **Check:** render the article (`TMPDIR` on the annex drive).
- **Done when:** a PR into `dev` is open for Andrea's review.

### 2. **P3** [claude] A fused pass for supervised fits

- **Where:** `update_q_impl()` in `src/radmix_kernels.cpp` and the
  `speculate` logic of `radmix_fit_cpp()` (`src/radmix_fit.cpp`).
- **Idea:** run the full Q side (every individual), which gives L(q, p),
  then restore the columns of the references; it reads the references'
  genotypes again, so it pays only when they are few.
- **Done when:** supervised fits take 4 passes per accepted iteration, with
  tests that their estimates are unchanged (to rounding).

### 3. **P3** [cleanup] Branches and scratch files

- laptopX1: the local branches in "Branches" (`git branch -D`), the
  worktree (`git worktree remove
  /media/andrea/annex/tmp/claude/tidygenclust/wt-coverage`), and the
  untracked `data-raw/radmix/perf_verification/` (a copy of what the tag
  `archive/radmix-perf-verify` holds).
- manica01: delete `~/tgc_newspeed_work` (2.2 GB of builds and cached
  data), `~/tgc_bench_article` (logs of the 16 Nov run; summaries are in
  `dev_notes/radmix_port_timing_2026-11/`) and the untracked
  `vignettes/articles/data/hgdp/runs` (4 GB of raw outputs of the article
  run; the article's data are committed). `tgc-desktop-checks.timer` is
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
- Genotype likelihoods (claude_jobs#3, merged from #12): cross-validation, SVD start, the
  1e-5 bounds, a faster source in tidypopgl.
- Move the notes still in `data-raw/radmix/` (`Notes_to_do`, `gl/REPORT.md`)
  into `dev_notes/`.

# tidygenclust: handover

Updated 2026-11-16 (manica01, night): the port is timed against `dev` and the benchmark article shows the full HGDP run; PR #36 waits only for review. For the next session: read this and
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `newspeed_port` | the port of `newspeed` + `supervised` onto `dev`, with the decisions of 15 Nov; also PR #37's test changes and the covr skip (16 Nov) | `dev` a20eb08 | yes | no (PR #36) |
| `coverage-speed` | covr at `-O1` (dropped), robust adamixture thread test, covr skip | `dev` a20eb08 | yes (abe2d2a local only) | superseded by `newspeed_port`; PR #37 closed 16 Nov; worktree `/media/andrea/annex/tmp/claude/tidygenclust/wt-coverage` |
| `newspeed` | per-iteration speed-ups, CV refit options, subset start; kits and desktop results | `dev` b5999cc | yes | superseded by `newspeed_port`; archived as `archive/newspeed` (ec96315) |
| `supervised` | supervised analyses; validation scripts and desktop results | `dev` b5999cc | yes | superseded by `newspeed_port`; archived as `archive/supervised` (fc43f9b) |
| `newspeed_supervised` | the merge of the two, before #35 | `newspeed` | yes | superseded by `newspeed_port` |
| `articles_revised` | revision of the pseudohaploid article (2 commits) | `dev` (old) | yes | no: conflict in `pseudohaploids.Rmd` |
| `gl`, `simplify` (origin), `gl_gl16`, `gl_private` (local) | merged as #35 and #34, or earlier versions of them | | | to delete |
| `radmix-speedups`, `radmix-perf-verify`, `kernel-speedups`, `svd-start` (local) | merged in #32/#33; origin branches gone | | | to delete |

## Next tasks

### 1. **P1** [review] Merge PR #36

- **Ready:** the port is 0.78-0.92 of `dev`'s time
  (`radmix_performance.md`), the benchmark article shows the full HGDP
  run (5eb3c73), and both R CMD checks pass on the desktop (16 Nov), see the comment on PR #36; only
  `test_gt_adamixture.R` was not run there (python `adamixture` not
  installed).
- **Read first:** the PR description; `dev_notes/radmix_performance.md`
  ("The changes of 13 November 2026", "Cross-validation refits");
  `dev_notes/covr_timing_2026-11.md`.
- **What Andrea checks:** the decisions as applied (`cv_rtol =
  max(rtol, 1e-7)`, `MIN_Q_TASKS` = 128); supervised analyses of genotype
  likelihoods, which the port added; the covr skip of
  `test_radmix_vs_admixture.R` and the 2-thread fold tests.
- **After the merge:** delete `newspeed`, `supervised`,
  `newspeed_supervised` and `coverage-speed` on origin and locally (the
  tags keep the records of the first two), and the worktree
  (`git worktree remove /media/andrea/annex/tmp/claude/tidygenclust/wt-coverage`);
  update `progress.md` and claude_jobs#3.
- **Gotcha:** to reproduce a bounds-check abort locally (GCC 13 here does
  not enable the checks), install a source copy with `R_MAKEVARS_USER`
  pointing to `CXXFLAGS = -g -O1 -D_GLIBCXX_ASSERTIONS` (and `CXX17FLAGS`,
  `CXX20FLAGS`), and run the test file under
  `R -d gdb --debugger-args="-batch -x cmds.gdb"` with `run ...` and `bt`.
- **Gotcha:** heavy multithreaded tests are ~40x slower under covr (gcov
  counters): skip them there (`skip_on_covr()`). Tests that run on CRAN use
  at most 2 threads.

### 2. **P2** [claude] then [review] `articles_revised`

- **Branch:** `articles_revised`; rebase or merge `dev` into it and resolve
  the conflict in `vignettes/articles/pseudohaploids.Rmd`.
- **Check:** render the article (`TMPDIR` on the annex drive).
- **Done when:** a PR into `dev` is open for Andrea's review.

### 3. **P3** [claude] A fused pass for supervised fits

- **Where:** `update_q_impl()` in `src/radmix_kernels.cpp` and the
  `speculate` logic of `radmix_fit_cpp()` (`src/radmix_fit.cpp`).
- **Idea:** run the full Q side (every individual), which gives L(q, p),
  then restore the columns of the references; it reads the references'
  genotypes again, so it pays only when they are few.
- **Done when:** supervised fits take 4 passes per accepted iteration, with
  tests that their estimates are unchanged (to rounding).

### 4. **P3** [cleanup] Branches and scratch files

- Origin: `gl`, `simplify`; after task 1, `newspeed`, `supervised`,
  `newspeed_supervised`, `coverage-speed` and its worktree.
- Local: `gl_gl16`, `gl_private`, `radmix-speedups`, `radmix-perf-verify`,
  `kernel-speedups`, `svd-start`; the untracked
  `data-raw/radmix/perf_verification/` (a copy of what the tag
  `archive/radmix-perf-verify` holds).
- manica01: `systemctl --user disable --now tgc-desktop-checks.timer`, and
  delete `~/tgc_newspeed_work` and `~/tgc_bench_article` (the builds and raw
  outputs of the 15 and 16 Nov runs; the summaries are in the repository).

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

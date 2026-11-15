# tidygenclust: handover

Updated 2026-11-15 (laptopX1), after PR #37. For the next session: read this and
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `newspeed_port` | the port of `newspeed` + `supervised` onto `dev`, with the decisions of 15 Nov | `dev` a20eb08 | yes | no (PR #36) |
| `coverage-speed` | covr at `-O1` with libstdc++ bounds checks; robust adamixture thread test | `dev` a20eb08 | yes | no (PR #37); worktree `/media/andrea/annex/tmp/claude/tidygenclust/wt-coverage` |
| `newspeed` | per-iteration speed-ups, CV refit options, subset start; kits and desktop results | `dev` b5999cc | yes | superseded by `newspeed_port`; archived as `archive/newspeed` (ec96315) |
| `supervised` | supervised analyses; validation scripts and desktop results | `dev` b5999cc | yes | superseded by `newspeed_port`; archived as `archive/supervised` (fc43f9b) |
| `newspeed_supervised` | the merge of the two, before #35 | `newspeed` | yes | superseded by `newspeed_port` |
| `articles_revised` | revision of the pseudohaploid article (2 commits) | `dev` (old) | yes | no: conflict in `pseudohaploids.Rmd` |
| `gl`, `simplify` (origin), `gl_gl16`, `gl_private` (local) | merged as #35 and #34, or earlier versions of them | | | to delete |
| `radmix-speedups`, `radmix-perf-verify`, `kernel-speedups`, `svd-start` (local) | merged in #32/#33; origin branches gone | | | to delete |

## Next tasks

### 1. **P1** [review] Review and merge PR #36 (`newspeed_port`)

- **Branch:** `newspeed_port`, PR #36 into `dev`.
- **Read first:** the PR description; `dev_notes/radmix_performance.md`
  ("The changes of 13 November 2026", "Cross-validation refits").
- **What Andrea checks:** the decisions as applied (`cv_rtol =
  max(rtol, 1e-7)`, `MIN_Q_TASKS` = 128); supervised analyses of genotype
  likelihoods, which the port added (tests: the fixed-individual kernel
  test and a supervised fit of hard-call likelihoods against genotypes).
- **CI:** every check passes, `test-coverage` included since 3192e45
  (the out-of-bounds fix in `set_fixed()`, `src/radmix_fit.cpp`).
- **After the merge:** delete `newspeed`, `supervised` and
  `newspeed_supervised` on origin and locally (the tags keep their
  records); update `progress.md` and claude_jobs#3.

### 2. **P2** [review] Review and merge PR #37 (`coverage-speed`)

- **Branch:** `coverage-speed`, PR #37 into `dev`; all checks pass, and
  Codecov accepts uploads again (the repository was deactivated on Codecov
  from 5-12 Nov until Andrea reactivated it on 15 Nov).
- **What Andrea checks:** covr at `-O1` (coverage 96.77% against 96.76%
  at `-O0`) with `-D_GLIBCXX_ASSERTIONS` kept: the coverage job is the
  only CI job with bounds checks, which caught the bug in #36; the
  rewritten thread test of `test_gt_adamixture.R`.
- **After the merge:** delete `coverage-speed` and its worktree
  (`git worktree remove /media/andrea/annex/tmp/claude/tidygenclust/wt-coverage`).
- **Gotcha:** to reproduce a bounds-check abort locally (GCC 13 here does
  not enable the checks), install a source copy with `R_MAKEVARS_USER`
  pointing to `CXXFLAGS = -g -O1 -D_GLIBCXX_ASSERTIONS` (and `CXX17FLAGS`,
  `CXX20FLAGS`), and run the test file under
  `R -d gdb --debugger-args="-batch -x cmds.gdb"` with `run ...` and `bt`.

### 3. **P2** [claude] then [review] `articles_revised`

- **Branch:** `articles_revised`; rebase or merge `dev` into it and resolve
  the conflict in `vignettes/articles/pseudohaploids.Rmd`.
- **Check:** render the article (`TMPDIR` on the annex drive).
- **Done when:** a PR into `dev` is open for Andrea's review.

### 4. **P3** [desktop] Time the port against `dev` and `newspeed`

See "Desktop queue".

### 5. **P3** [claude] A fused pass for supervised fits

- **Where:** `update_q_impl()` in `src/radmix_kernels.cpp` and the
  `speculate` logic of `radmix_fit_cpp()` (`src/radmix_fit.cpp`).
- **Idea:** run the full Q side (every individual), which gives L(q, p),
  then restore the columns of the references; it reads the references'
  genotypes again, so it pays only when they are few.
- **Done when:** supervised fits take 4 passes per accepted iteration, with
  tests that their estimates are unchanged (to rounding).

### 6. **P3** [cleanup] Branches and scratch files

- Origin: `gl`, `simplify`; after task 1, `newspeed`, `supervised`,
  `newspeed_supervised`; after task 2, `coverage-speed` and its worktree.
- Local: `gl_gl16`, `gl_private`, `radmix-speedups`, `radmix-perf-verify`,
  `kernel-speedups`, `svd-start`; the untracked
  `data-raw/radmix/perf_verification/` (a copy of what the tag
  `archive/radmix-perf-verify` holds).
- manica01: `systemctl --user disable --now tgc-desktop-checks.timer`, and
  delete `~/tgc_newspeed_work` once the timing check (task 4) is done.

## Desktop queue

### Timing of the port (P3)

- **Branch and set-up:** on manica01, `~/tgc_newspeed_work` still has the
  builds and the cached data of the run of 15 Nov. The kit is in the tag:
  `git worktree add ../tgc_kit archive/newspeed`, then run from that
  worktree, which builds any commit given as `label=commit`.
- **Run:** `MAKEFLAGS=-j8 PIN=0-7 bash $K/build_variants.sh
  dev=a20eb08 port=newspeed_port`, then `MAX_ITER=8
  CASES=hgdp:5,hgdp:10 THREADS="1 8" ROUNDS=3 bash $K/run_bench.sh
  $TGC_WORK/results/port dev final port` (`K=data-raw/radmix/newspeed`,
  `TGC_WORK=$HOME/tgc_newspeed_work`; 8 P-cores, about 30 min,
  timing-sensitive: idle machine, warn the other sessions).
- **Results:** `$TGC_WORK/results/port/summary.md`; expected: `port` close
  to `final` (0.76-0.80 of `dev`'s time before #35) apart from
  `MIN_Q_TASKS` (3-6% at 8 threads).
- **Report:** a table in `radmix_performance.md`, a progress entry, and a
  log comment in claude_jobs#3.

## Backlog

- Where the covr overhead goes: the tests take 25-45 min under covr in CI
  (noisy runners, one run per setup) against 4.2 min in R CMD check. Time
  each test file under covr against a plain run on a quiet machine (the
  `benchmarking` skill), to find out whether the R instrumentation, the
  gcov counters or a few heavy tests dominate. **Running overnight
  15-16 Nov on laptopX1** (Andrea's request): systemd user unit
  `tgc-covr-timing` (`overnight.sh` in
  `/media/andrea/annex/tmp/claude/tidygenclust/covrtime/`), starts after
  00:00 once idle, three balanced rounds of covr at -O0 and -O1 (OMP 4 and
  1), and plain -O0, -O0 --coverage and -O2 builds; no run starts after
  07:30. Results: `results/summary.md` and `results/status.md` there.
  Then: choose the coverage flags, apply them to the workflow, and push
  with the unpushed dev_notes commit on `newspeed_port`.
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
- Genotype likelihoods (claude_jobs#12): cross-validation, SVD start, the
  1e-5 bounds, a faster source in tidypopgl.
- Move the notes still in `data-raw/radmix/` (`Notes_to_do`, `gl/REPORT.md`)
  into `dev_notes/`.

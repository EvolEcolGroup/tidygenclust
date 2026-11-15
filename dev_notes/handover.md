# tidygenclust: handover

Updated 2026-11-15 (laptopX1). For the next session: read this and
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `newspeed_port` | the port of `newspeed` + `supervised` onto `dev`, with the decisions of 15 Nov | `dev` a20eb08 | yes | no (PR #36) |
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
- **After the merge:** delete `newspeed`, `supervised` and
  `newspeed_supervised` on origin and locally (the tags keep their
  records); update `progress.md` and claude_jobs#3.

### 2. **P2** [claude] Find out why `test-coverage` has failed on `dev` since #34

- **Read first:** the GitHub Actions logs of `test-coverage` on `dev`
  (`gh run list --branch dev --workflow test-coverage`), and
  `.github/workflows/`.
- **Done when:** the cause is known and fixed on a branch with a PR, or
  reported to Andrea if it is outside the package (e.g. covr and OpenMP or
  the ADMIXTURE binary).

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
  `newspeed_supervised`.
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

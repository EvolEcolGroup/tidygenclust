# tidygenclust: handover

Updated 2026-11-28 (laptopX1): 0.2.0 is on `main` (#44) and `dev`, ready
for the first CRAN submission: `main`'s #40 merged into `dev` (#43),
tidypopgen (>= 1.1.0) from CRAN without `Remotes`, `cran-comments.md`.
`R CMD check --as-cran` gives only the new-submission note; CI green.
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `hgdp_plots` | Eirlys's HGDP plotting script (Jan 2025); `dev` has her later version (June 2025) | old `dev` | yes | no: kept until some runs are redone (Andrea, 28 Nov) |
| `windows` | Ondra's experimental Windows support through WSL2 (`R/wsl_backend.R`, `vignettes/windows.Rmd`, a CI workflow, tests; Oct 2026); its adamixture commits are in `dev` (#30) | `dev` (June) | yes | no PR |

`main` and `dev` have the same tree (0.2.0; `main` has one extra merge
commit).

Deleted 28 Nov (Andrea's go-ahead): `ccr-fd217561-izn297` and
`adamixture-1.9.6-benchmark` (origin; the first is in `dev` through #41,
the second, PR #39, was integrated in the article's rewrite), and the
local `ada196_article`, `radmix_review_fixes` and
`adamixture-1.9.6-benchmark` (all inside #41).
Deleted 20 Nov: the local `radmix_extensions` and
`data-raw/radmix/perf_verification/` (its one file not in the tag
`archive/radmix-perf-verify` moved to
`/media/andrea/annex/tmp/claude/tidygenclust/kept/`). Deleted 18 and 17
Nov: see `progress.md`.

## Next tasks

### 1. **P2** [claude] examples for the exported functions

- On a branch from `dev`: add `@examples` to the exported functions that
  have none (`gt_radmix_data()`, `autoplot_gt_clumppling()`,
  `subset_gt_clumppling()`, `tidy_gt_clumppling()`, `gt_fastmixture()`,
  `gt_adamixture()`, the `tgc_*` installers and `tgc_tools_version()`),
  with `\donttest{}` for those that need python; CRAN usually asks for them
  on a first submission. Bump to 0.2.0.9000 (or 0.2.1 if resubmitting).
- **Done when:** every exported function has examples and
  `R CMD check --as-cran` stays clean.

### 2. **P2** [andrea] submit 0.2.0 to CRAN

- `devtools::submit_cran()` from `main` (after task 1 if wanted);
  `cran-comments.md` is ready. After acceptance: tag the release, then
  bump `dev` to 0.2.0.9000 (re-add `Remotes: EvolEcolGroup/tidypopgen@dev`
  only if `dev` needs unreleased tidypopgen).

### 3. **P3** [cleanup] the article's data on manica01

- `~/project_temp/claude/tidygenclust/stream/` (about 75 GB: builds,
  benchmark data, the article's `mid/` and `large/`): delete when Andrea
  no longer needs them.

## Desktop queue

- Nothing queued. On manica01, `~/project_temp/claude/tidygenclust/stream/`
  holds the builds (`work/`), the 1.5 GB benchmark data (`data/`), the
  article's data (`mid/`, 6 GB; `large/`, 62 GB) and the driver: delete
  them once the article is final.

## Backlog

- The chunk size of the Q side in memory at large K (10, 20): 16 against
  64 MB on 10,000 x 150,000 (`radmix_streaming.md`, "Chunks in memory").
- A real hard disk or a cluster file system for streaming (manica01's
  disks are not mounted; the throttle caps throughput, not latency).
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

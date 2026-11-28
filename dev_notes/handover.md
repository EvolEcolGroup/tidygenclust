# tidygenclust: handover

Updated 2026-11-28 (laptopX1): `ada196_article` was squash-merged into
`dev` as #41 (27 Nov, 6e33825): streaming of data larger than the memory,
the Q side in chunks in memory, the review fixes (`gl_curvature`), the
register-blocked Hessians and the benchmark article with its section on
data larger than the memory. `dev` is at 0.1.5.9010 (#42 after it).
tidypopgl's prefetch function for streaming is in its `dev` (#91), which
requires tidygenclust >= 0.1.5.9010. `summary.md` first. History is in
`progress.md`; designs are in the topic documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `hgdp_plots` | Eirlys's HGDP plotting script (Jan 2025); `dev` has her later version (June 2025) | old `dev` | yes | no: kept until some runs are redone (Andrea, 28 Nov) |
| `windows` | Ondra's experimental Windows support through WSL2 (`R/wsl_backend.R`, `vignettes/windows.Rmd`, a CI workflow, tests; Oct 2026); its adamixture commits are in `dev` (#30) | `dev` (June) | yes | no PR |

`main` holds #40 (`articles_revised`, the revised articles, merged into
`main` on 27 Nov) and two "Update conda snapshot" commits that `dev` lacks.

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

### 1. **P2** [claude] bring `main`'s #40 into `dev`

- `git merge origin/main` into a branch from `dev` (#40 revised
  `vignettes/articles/pseudohaploids.Rmd` and the vignette; the conda
  snapshot commits), resolve the conflicts, render the articles (`TMPDIR`
  in a work folder), then a PR into `dev`.
- **Done when:** `main` is an ancestor of `dev` again.

### 2. **P3** [cleanup] the article's data on manica01

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

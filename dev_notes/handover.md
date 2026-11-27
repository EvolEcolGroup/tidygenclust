# tidygenclust: handover

Updated 2026-11-27 (laptopX1, with runs on manica01): streaming is merged
into the local `ada196_article` (version 0.1.5.9007), measured on the
desktop (`radmix_streaming.md`, "Desktop") and shown in the benchmark
article's section "Data larger than the memory": 5 GB under 4 and 1 GiB
caps, and 50 GB under 8 GiB (34% longer than in memory on the NVMe, 4.3
times at 500 MB/s), identical results. Andrea kept the streaming defaults
and asked for the Q side in chunks in memory too (26 Nov): done (c521e36,
79e21c4), 15-21% faster with 8,192 individuals or more; the article's
fits in memory were rerun with it (4503fdd). The article waits for
Andrea's review. `summary.md` first. History is in `progress.md`; designs
are in the topic documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `ada196_article` | the benchmark article rerun with ADAMIXTURE 1.9.6, the review fixes (`gl_curvature`), the register-blocked Hessians, and (merged 24 Nov) streaming, with its desktop benchmarks and the scripts of the article's section on data 5 times the memory | `dev` | yes (27 Nov, Andrea: a branch, no PR) | no |
| `ccr-fd217561-izn297` | prototype: `gt_radmix()` streams data larger than the memory (chunks, prefetch, alternating directions), with tidypopgl's branch of the same name (its prefetch function) | `dev` (`acff37a`) | yes | into `ada196_article` (349c929) |
| `radmix_review_fixes`, `adamixture-1.9.6-benchmark` | local branches merged into `ada196_article` | `dev` | no | into `ada196_article` |
| `articles_revised` | revision of the pseudohaploid article (2 commits) | `dev` (old) | yes | no: conflict in `pseudohaploids.Rmd` |
| `hgdp_plots`, `windows` (origin) | not part of the radmix work; commits not in `dev` | | yes | no: ask Andrea before touching |

Deleted 20 Nov: the local `radmix_extensions` and
`data-raw/radmix/perf_verification/` (laptopX1; its one file not in the tag
`archive/radmix-perf-verify`, `svd_vs_random_q_hgdp.rds`, moved to
`/media/andrea/annex/tmp/claude/tidygenclust/kept/`), and the scratch
folders of 16-18 Nov and the article's `hgdp/runs` (manica01).
Deleted 18 Nov (laptopX1): 17 local branches without origin and the
`coverage-speed` worktree; tips in
`/media/andrea/annex/tmp/claude/tidygenclust/deleted_branches_2026-11-18.txt`.
Deleted 17 Nov (origin and manica01): `newspeed`, `supervised` (tags
`archive/newspeed` ec96315, `archive/supervised` fc43f9b),
`newspeed_supervised`, `coverage-speed`, `newspeed_port` (#36), and the
local `adamixture`, `rewrite`, `openadmixture`, `rclump`.

## Next tasks

### 0. **P1** [review] the article's section on data larger than the memory

Read `dev_notes/radmix_streaming.md`, "The article's runs" and "Chunks in
memory", first. Scripts and raw results: `data-raw/radmix/stream/desktop/`
and `data-raw/radmix/stream/results/desktop_*`; the article's results are
in `vignettes/articles/data/admixture_benchmark/{mid,large}_*.csv` (the
`memory` rows from the rerun of 27 Nov, `phases/inmem_cq.sh` on manica01).

- **[review] (Andrea):** the section "Data larger than the memory" of
  `vignettes/articles/admixture_benchmark.Rmd` (results after "SVD and
  random starts"; code under "How the benchmark was run"; summary bullet;
  a line in the introduction and in the memory section). Rendered copy:
  `/media/andrea/annex/tmp/claude/tidygenclust/render_test/out/admixture_benchmark.html`
  (laptopX1). To render again: from `vignettes/articles`,
  `R_LIBS=<a library with tidypopgen 1.0.0.9015> TMPDIR=<work folder>
  Rscript -e 'rmarkdown::render("admixture_benchmark.Rmd")'`. The main
  benchmark (HGDP, 1,043 individuals: chunks in memory need 8,192) is not
  rerun.
- **Then [claude]:** a PR into `dev` only when Andrea says so
  (`ada196_article` is pushed, without a PR, by his instruction). Delete the article's
  data on manica01 when he no longer needs them (below).
- **Done when:** the article is reviewed.

### 1. **P2** [claude] then [review] `articles_revised`

- **Branch:** `articles_revised`; rebase or merge `dev` into it and resolve
  the conflict in `vignettes/articles/pseudohaploids.Rmd`.
- **Check:** render the article (`TMPDIR` in a work folder).
- **Done when:** a PR into `dev` is open for Andrea's review.

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

# tidygenclust: handover

Updated 2026-11-25 (laptopX1, with runs on manica01): the streaming of
`ccr-fd217561-izn297` is merged into the local `ada196_article` and
measured on the desktop (`radmix_streaming.md`, "Desktop"): identical
estimates everywhere, no slowdown in memory except 3-4% on small data on
one thread, 1.1-1.6 times the time in memory out of core on the NVMe
(`dev` up to 16 times). The prefetch read 1/16 of each chunk on the
desktop and is fixed (b07f087). Andrea decides next what to keep of the
streaming changes; the article's section on data 5 times the memory waits
for that. `summary.md` first. History is in `progress.md`; designs are in
the topic documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `ada196_article` | the benchmark article rerun with ADAMIXTURE 1.9.6, the review fixes (`gl_curvature`), the register-blocked Hessians, and (merged 24 Nov) streaming, with its desktop benchmarks and the scripts of the article's section on data 5 times the memory | `dev` | **no: local only** (laptopX1; on manica01 by git bundle). Andrea: keep it local | no |
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

### 0. **P1** [decide] then [claude] streaming of data larger than the memory

Read `dev_notes/radmix_streaming.md`, "Desktop", first. Measured on
manica01 (8 P-cores, NVMe, also throttled to 500 and 150 MB/s): see the
tables there. Scripts and raw results: `data-raw/radmix/stream/desktop/`
and `data-raw/radmix/stream/results/desktop_*`.

- **Decisions (Andrea):**
  1. Keep the streaming changes as they are (stream whenever the file
     exceeds 0.8 times the memory available, few individuals included:
     never worse than `dev` by more than 3% on the desktop).
  2. The prefetch on by default (fixed: 20-27% faster with many
     individuals on the NVMe; 7-9% slower with few individuals and a third
     of the file cached), and the alternation (needed: 22-51% faster on
     the slow disks).
  3. Chunked Q side in memory too (15% faster with 10,000 individuals, 3%
     slower with 1,000).
  4. Test on a real hard disk (one of manica01's, mounted inside
     `~/project_temp/claude/`, which is exempt from the antivirus; needs
     sudo) and on a cluster's network file system before the article.
- **[claude] the 3-4% on small data on one thread** (eur, `update_p` and
  `loglik`; not hgdp, not 8 threads): find the cause with `perf stat`
  (instructions and cycles of base 51b4e9f against the merge) on
  manica01, when no timing runs; the libraries are in
  `~/project_temp/claude/tidygenclust/stream/work/libs/{base,new,new2}`.
- **[desktop] then the article** (Andrea, 25 Nov: 250,000 individuals,
  three fits: in memory, capped at 8 GiB on the NVMe, capped and throttled
  to 500 MB/s; memory records shown): simulate 250,000 x 200,000 (50 GB
  backing file of the gen_tibble, read from a 12.5 GB bed) and fit it
  (`vignettes/articles/code/admixture_benchmark_large_data.R`, then
  `admixture_benchmark_large.sh`, from `vignettes/articles`); the driver's
  phases are `large_smoke`, `large_data`, `large_fits` in
  `data-raw/radmix/stream/desktop/phases/` (build `art`). Then the section
  "Data larger than the memory" of `admixture_benchmark.Rmd` (after "SVD
  and random starts"), its summary bullet, and a subsection of "How the
  benchmark was run" reading the chunks `large-*` of the three scripts;
  render the article.
- **Done when:** the decisions are taken, the article section is written
  and rendered, and `ada196_article` is committed locally (no push, by
  Andrea's instruction).

### 1. **P2** [claude] then [review] `articles_revised`

- **Branch:** `articles_revised`; rebase or merge `dev` into it and resolve
  the conflict in `vignettes/articles/pseudohaploids.Rmd`.
- **Check:** render the article (`TMPDIR` in a work folder).
- **Done when:** a PR into `dev` is open for Andrea's review.

## Desktop queue

- Task 0's article runs (after Andrea's decisions): about 3-4 h, the
  simulation on 16 cores and the fits on the P-cores 0-7; no sudo.
- If Andrea wants it: a hard disk test (mount one of manica01's disks
  inside `~/project_temp/claude/`, sudo), the phases `ooc_third_pf` and
  `ooc_pf` of `data-raw/radmix/stream/desktop/` with `BENCH_DATA` on it
  (about 30 min, plus the `dev` runs if wanted).

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

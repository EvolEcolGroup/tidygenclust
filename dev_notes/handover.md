# tidygenclust: handover

Updated 2026-11-24 (cloud session): the streaming prototype of
`gt_radmix()` for data larger than the memory is on `ccr-fd217561-izn297`
(and tidypopgl's branch of the same name), measured in the cloud
container; task 0 is its measurement on the desktop. Before that
(2026-11-20, laptopX1): the bounds of likelihood estimates stay
`[1e-5, 1 - 1e-5]` (Andrea); the scratch files and the local
`radmix_extensions` are deleted on both machines.
`summary.md` first. History is in `progress.md`; designs are in the topic
documents listed in `README.md`.

## Branches

| branch | purpose | base | pushed | merged |
|---|---|---|---|---|
| `articles_revised` | revision of the pseudohaploid article (2 commits) | `dev` (old) | yes | no: conflict in `pseudohaploids.Rmd` |
| `hgdp_plots`, `windows` (origin) | not part of the radmix work; commits not in `dev` | | yes | no: ask Andrea before touching |
| `ccr-fd217561-izn297` | prototype: `gt_radmix()` streams data larger than the memory (chunks, prefetch, alternating directions; `dev_notes/radmix_streaming.md`), with tidypopgl's branch of the same name (its prefetch function) | `dev` (`acff37a`) | yes | no: task 0 (desktop) first |

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

### 0. **P2** [claude] then [review] streaming of data larger than the memory (desktop)

Read `dev_notes/radmix_streaming.md` first: the design, step 1 (what each
pass of `dev` costs out of core), the full-fit matrices (67%, 92% and 33% of
the file cacheable) and the cap of the Q side's prefetch. In short, from
the cloud container (4 cores, a disk backed by the host's cache, noisy):
with 8,192 individuals or more, `dev`'s Q-side passes read the file
`n_blocks / n_cores` times (40 times for 10,000 individuals on 4 threads)
and a fit takes 21-25 times its time in memory, against 1.1-1.4 times when
streaming; with fewer individuals `dev` is already within 1.1-1.5 times,
and streaming only helps on slow disks. The container could not resolve
differences under about 5%, and has no real disks: the desktop decides.

- **Branches:** `ccr-fd217561-izn297` of tidygenclust (base `dev` at
  `acff37a`) and of tidypopgl (its source's prefetch function; base `dev`).
- **Setup (manica01):**
  - `git fetch origin && git switch ccr-fd217561-izn297` in
    `~/git/tidygenclust`; the timing kit is the tag `archive/newspeed`
    (`git worktree add ../tidygenclust_newspeed archive/newspeed` if that
    worktree is gone).
  - Builds, from the kit's worktree:
    `bash data-raw/radmix/newspeed/build_variants.sh dev=acff37a
    stream=origin/ccr-fd217561-izn297` (into `$TGC_WORK/libs/{dev,stream}`;
    `TGC_WORK` defaults to `../tgc_newspeed_work`).
  - The scripts of `data-raw/radmix/stream/` (README there) take
    `BENCH_DATA` (the data, on the disk under test), `LIB_BASE` and
    `LIB_NEW` (the two libraries), `N_CORES=8`. They cap memory and read
    throughput with `run_capped.sh`, which on the desktop uses
    `sudo systemd-run --scope -p MemoryMax=... -p IOReadBandwidthMax=...`
    with R run as the user, and drops the page cache with sudo. That path
    was written in the container (cgroup v1) and is **untested**: check it
    first with `sh run_capped.sh 300 0 <file> "" cat /proc/self/cgroup`
    (it should print a `run-*.scope`) and one small `bench_stream.sh` run
    (a 1.5 GB file with a 600 MB cap should read several GB from disk).
  - sudo in an unattended run: either start the matrix in a terminal after
    `sudo -v` with a keep-alive (`while sleep 120; do sudo -n true; done &`),
    or allow `systemd-run` and the cache drop without a password (Andrea's
    call: a password-less `systemd-run` is root access).
  - Data: `Rscript data-raw/radmix/stream/make_data.R 10000 150000 big_a`
    and `... 1000 1500000 big_b` in `BENCH_DATA` (1.5 GB each, about 2 min),
    then the peak memory of the scripts: `sh measure_peak.sh anon_big_a.txt
    Rscript bench_fit.R big_a TRUE 8 4 2 1 64 peak` (and `big_b`;
    `anonp_<prefix>.txt` with `bench_passes.R big_a TRUE 8 4 1 64 peak`).
- **Steps, in order:**
  1. Conventional fits (the "no slowdown" check, quiet machine): the kit's
     `run_bench.sh` as for #36, `CASES=hgdp:5,hgdp:10,eur:3,eur:4
     THREADS="1 8" ROUNDS=3 MAX_ITER=1000 bash
     data-raw/radmix/newspeed/run_bench.sh <out> dev stream` (data that fit
     in memory never stream, so this times the usual layout of the branch;
     the kit converts and caches HGDP and EUR in `$TGC_WORK/data`), and
     `MODE=kernels` for the passes. A slowdown above about 5% is asked
     about (the only change in the usual layout is three arguments of
     `accum_q_block()`).
  2. Step 1 on the real disks: `N_CORES=8 sh run_passes.sh passes.csv
     "mem:0 33:0 67:0 92:0"` with `BENCH_DATA` on the NVMe, then on a hard
     disk if there is one (e.g. `/media/andrea/data1`, if it is not an SSD);
     `Rscript summarise_passes.R passes.csv`. The question: are `dev`'s
     passes bound by the disk or by computation on 8 cores, and from what
     fraction cached.
  3. Full fits: `N_CORES=8 sh run_ooc_third.sh third.csv 32` (33%) and
     `N_CORES=8 sh run_ooc.sh ooc.csv 2 1 64 ""` (67% and 92%, `FRACS`) on
     each disk; the throttled rows (500 and 150 MB/s) emulate slower disks
     on the NVMe. The `dev` runs with 10,000 individuals take about 10 min
     each.
  4. A real file larger than the desktop's memory, no caps: `make_data.R`
     with, e.g., 20,000 individuals and enough loci for 1.5 times the RAM;
     one fit with each build (`bench_fit.R`, `TRUE`/`base`; or
     `gt_radmix()` with the default options, checking that
     `radmix_stream_mode()` chooses to stream and the chunk size).
  5. Decide, with Andrea: (i) stream only with one Q segment (8,192
     individuals or more, the case that needs it) or whenever the file
     exceeds 0.8 times the memory; (ii) prefetch and alternation on by
     default (alternation needs the prefetch; with few individuals, no
     prefetch was 7% faster on a 150 MB/s disk); (iii) the chunk size (a
     sixteenth of the memory, 16-256 MB); (iv) whether the chunked Q side
     of one segment becomes the default in memory too (it was 15-20%
     faster there in the container).
- **Then:** the follow-ups of "Left to do" in `radmix_streaming.md` that are
  wanted (the SVD start of likelihoods re-reads the file; Windows has no
  prefetch and never streams by itself), the version bump, and a PR into
  `dev`; tidypopgl's branch after it (it works with any tidygenclust: the
  prefetch elements are optional).
- **Done when:** the decisions of step 5 are made and a PR into `dev` is
  open for Andrea's review.

### 1. **P2** [claude] then [review] `articles_revised`

- **Branch:** `articles_revised`; rebase or merge `dev` into it and resolve
  the conflict in `vignettes/articles/pseudohaploids.Rmd`.
- **Check:** render the article (`TMPDIR` in a work folder).
- **Done when:** a PR into `dev` is open for Andrea's review.

## Desktop queue

Task 0, steps 1-4 (about 4-5 hours, most of it the `dev` runs with 10,000
individuals out of core), once `run_capped.sh` is checked on manica01 and
the sudo question is settled. A driver can follow the last unattended run's
(`dev_notes/radmix_port_timing_2026-11/run_desktop.sh`: waits for an idle
machine, one phase at a time, a status file).


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

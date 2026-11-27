# tidygenclust: progress

The log of major activity, newest first. Never rewrite old entries: add a
correction as a new entry. Earlier history (before 15 Nov 2026) is in the
comments of dramanica/claude_jobs#3 (radmix) and #12 (genotype
likelihoods), and in `radmix_performance.md`.

## 2026-11-26/27 (laptopX1, manica01): chunks for the Q side in memory; the article's fits in memory again

- Andrea, 26 Nov: keep the streaming defaults (stream above 0.8 x the
  memory available, prefetch in pieces of 128 KB, alternation on), and
  take the Q side in chunks for data in memory too.
- c521e36: `chunk_q`, for a Q side of one segment (8,192 individuals or
  more) when the data are not streamed; identical results (tests with
  chunks of one byte, K = 2 and 8, with and without fixed individuals).
- manica01, quiet, 8 P-cores: chunks of 8-64 MB are 15-19% faster than
  the usual layout on 10,000 x 150,000, 50,000 x 100,000 and
  250,000 x 200,000, and within about 2% of each other; 64 MB kept
  (79e21c4, best on the largest data). NEWS item.
- The article's fits in memory rerun with 79e21c4 (manica01, quiet,
  06:24-07:02; 4503fdd): 5 GB 214 s (was 255), 50 GB 1,711 s (was 2,173),
  same log-likelihoods; the streamed runs, always chunked, now take
  1.29-1.55 times the time in memory on the NVMe (50 GB under 8 GiB:
  1.34, was 1.05). The first comparison understated the cost of the disk.
- Changed: `radmix_streaming.md` ("Chunks in memory", "The article's
  runs").

## 2026-11-25/26 (laptopX1, manica01): the article's section on data larger than the memory

- 524b347: the log-likelihood pass inlines `loglik_slice()` again (the
  merge had made it two calls per slice): 7-11% faster than before the
  merge on one thread; whole fits unchanged. The 3-4% slowdown of eur on
  one thread in the gate did not reproduce (5 rounds: equal).
- Version 0.1.5.9007 (870cff4). HGDP check on manica01: 5.07 s on 8
  threads and 40.2 s on 1, as in the article (Andrea: no rerun).
- Article (Andrea, 25 Nov: 250,000 individuals; a bed file; in memory,
  SSD and throttled): 5 GB (50,000 x 100,000) under 4 and 1 GiB caps on
  the NVMe, 500 and 150 MB/s; 50 GB (250,000 x 200,000) under 8 GiB on
  the NVMe and 500 MB/s. 1 GiB of 5 GB: 30% longer than in memory on the
  NVMe; 8 GiB of 50 GB: 5% longer, 3.4 times at 500 MB/s; identical
  results; R at most 0.40 GB. Section written and rendered (03c7161).
- Changed: `radmix_streaming.md` ("The article's runs").

## 2026-11-24/25 (laptopX1, manica01): streaming merged and measured on the desktop

- `ccr-fd217561-izn297` merged into the local `ada196_article` (349c929;
  the blocked Hessians combined with the streaming accumulators, tests at
  K = 2, 3, 7, 8 bit-identical: 24908ae). Nothing pushed (Andrea).
- manica01, 8 P-cores: in memory, identical results, equal times except
  eur on one thread, 3-4% slower (open). Out of core, identical estimates
  in 72 runs; streaming 1.1-1.6 times its time in memory on the NVMe,
  `dev` up to 16 times (it reads a 1.5 GB file 100 times).
- Fixed: the prefetch read 1/16 of each chunk (madvise pieces of 2 MB,
  128 KB honoured per call: b07f087), and `run_capped.sh` gave systemd a
  relative path for the throttle, which it ignores (df19269).
- The article's section on data 5 times the memory is postponed until
  Andrea decides on the streaming changes (`handover.md`, task 0).
- Changed: `radmix_streaming.md` ("Desktop", "Left to do"),
  `data-raw/radmix/stream/` (desktop scripts and results).

## 2026-11-20 (laptopX1): bounds of likelihoods decided; clean-up

- Andrea decided to keep `[1e-5, 1 - 1e-5]` for estimates from genotype
  likelihoods, as for genotypes, not NGSadmix's 1e-9 (recorded in
  `admixture_implementation.md`, "Genotype likelihoods").
- Deleted on laptopX1: the local `radmix_extensions` (ae6bb5c, the head of
  #38) and the untracked `data-raw/radmix/perf_verification/`; its
  `svd_vs_random_q_hgdp.rds` (97 KB, not in the archive tag) moved to
  `/media/andrea/annex/tmp/claude/tidygenclust/kept/`.
- Deleted on manica01: `~/tgc_newspeed_work`, `~/tgc_bench_article`,
  `~/tgc_supfused_work`, `vignettes/articles/data/hgdp/runs` (about 6.4 GB).

## 2026-11-20 (manica01): the website shows main at the root, dev at /dev/

- Since 9989231 (31 Oct) the pkgdown workflow also built pushes to `dev`
  into the root of `gh-pages`, so the landing page showed `dev`.
- ac79eb4 (`dev`): the workflow sets `PKGDOWN_DEV_MODE` from the branch
  (`dev` -> `/dev/`; main and releases -> root); `_pkgdown.yml` uses
  `development: mode: auto` for local builds. f600f3f: README links to
  the `/dev/` site.
- `gh-pages` root restored to main's last build (efec401, tree of 56c5d2a
  from main 0fc9c8b, 17 Oct) (Andrea). Main gets the new settings when
  `dev` is merged; until then its builds still go to the root, as wanted.
- Checked live: root 0.1.5 (no dev-only articles), `/dev/` 0.1.5.9006
  with `noindex` (deploy 4981fef touched only `dev/`).

## 2026-11-19 (manica01): tidypopgl's PR #47 merged

- Andrea squash-merged tidypopgl's #47 into `gl_workflow` (201f29f):
  `gl_workflow` now requires tidygenclust 0.1.5.9006. Its squash message
  held the notes' skip marker, so GitHub ran no CI on it; tidypopgl's
  60d956d runs the CI of its #46 on the merged code.

## 2026-11-19 (manica01): tidypopgl's `gl_workflow` with `dev`; version 0.1.5.9006

- tidypopgl's `gl_workflow` works with `dev` (b14fdeb) on all three CI
  platforms except one test, which still expected `gt_radmix(crossval =
  TRUE)` to stop on likelihoods (#38 made it work).
- `dev` bumped to 0.1.5.9006 (d889f63), so that tidypopgl can require the
  version with cross-validation and SVD starts for likelihoods (Andrea).
- tidypopgl PR #47 (`radmix_source_speed` into `gl_workflow`): requires
  0.1.5.9006, tests cross-validation through its source, and its vignette
  chooses `k` by cross-validation; full suite passed on manica01, CI on
  the PR running.
- `radmix_extensions` deleted on origin and manica01 (Andrea).

## 2026-11-18 (laptopX1): PR #38 merged

- CI of PR #38 green on all platforms; squash-merged into `dev` (df747aa)
  at Andrea's request.
- 17 local branches without origin and the `coverage-speed` worktree
  deleted on laptopX1 (Andrea's request); tips saved in
  `/media/andrea/annex/tmp/claude/tidygenclust/deleted_branches_2026-11-18.txt`.

## 2026-11-18 (laptopX1): `svd_power`, PR #38

- At Andrea's request, new `gt_radmix()` argument `svd_power` (default 3)
  for the power iterations of the randomised SVD of likelihoods (014eb30,
  branch `radmix_extensions`); stale `@param init` sentence fixed.
- `test_radmix_svd.R`, `test_gt_radmix.R`, `test_radmix_cv.R` pass (-O2).
- PR #38 (`radmix_extensions` into `dev`) opened for Andrea's review.

## 2026-11-18 (manica01, overnight): extensions of `gt_radmix()` (branch `radmix_extensions`)

- Unattended evening/night session (claude_jobs#3, P3 "Extensions"),
  branch `radmix_extensions` from `dev` 00deb4f, pushed, not merged:
  fused pass for supervised fits (ec7ff9a; fits 0.88-0.89 of the time,
  HGDP K = 5, 139 or 660 references, 1 and 8 threads, Q and P
  bit-identical); cross-validation (ecc0841) and SVD starts (4ed9ceb) for
  genotype likelihoods; obsolete test expectations removed (9b05912);
  `slice_gl()` and `loglik_gl_slice()` inlined into the AVX2 clones (9659e51).
- tidypopgl branch `radmix_source_speed` (pushed): a faster source; with
  both, likelihood fits take 0.80 (1 thread) and 0.83 (8) of the time.
- CV of likelihoods chose the simulated K in 12/12 cases; SVD starts of
  likelihoods converged faster and higher at K = 5.
- Checks: full suite NOT_CRAN (22 files, 3871 expectations), `R CMD check
  --as-cran` (0 errors, 0 warnings, 3 environment notes); tidypopgl's suite
  passes with the branch (44 files, 7119).
- Notes: `admixture_implementation.md`, `radmix_performance.md`,
  `radmix_extensions_2026-11/`.

## 2026-11-18 (manica01): attribution merged into `dev`

- Andrea approved the new `cph` roles (Ko and Lange restored; Ioannidis
  Lab, Meisner, Lopatin, Kassambara) and asked to merge: `code_attribution`
  (4c8464d attribution, d41248d `radmix_cv.R` comment, then these notes)
  fast-forwarded into `dev` and pushed; the branch is deleted.

## 2026-11-18 (manica01): provenance of the clumppling code; text for the paper

- Same audit for `gt_clumppling()` (branch renamed `code_attribution`,
  still uncommitted): the Hungarian solver is adapted from cp-algorithms
  (Andrey Lopatin, CC BY-SA 4.0, compatible with GPLv3) and `.get_legend()`
  is copied from ggpubr (GPL >= 2); neither was credited. Both added to
  `COPYRIGHTS` and DESCRIPTION (`cph`); file headers now name the source
  functions and licences of the clumppling, python-louvain and
  markov_clustering ports (comments only; builds). Fixed "anative" in
  `?gt_clumppling`.
- `radmix_methods_text.md` became `paper_methods_text.md`, with a short
  paragraph on `gt_clumppling()` (Andrea asked for one short paragraph).

## 2026-11-18 (manica01): provenance of the radmix code; text for the paper

- Branch `radmix_attribution` (from `dev` 54ce6cd, not committed yet):
  audited the code of `gt_radmix()` against OpenADMIXTURE.jl (401ac5c),
  ADAMIXTURE 1.8.9 and fastmixture 1.3.0. Translated from OpenADMIXTURE.jl
  (MIT): the QP solver, the projections, the quasi-Newton step, the
  structure of the fitting loop; adapted from ADAMIXTURE (BSD-3): the ALS
  of the SVD start; derived from fastmixture (GPL-3.0-only):
  `py_rfastmixture.py`. The kernels are original.
- Gaps fixed: ADAMIXTURE and fastmixture were missing from `COPYRIGHTS`
  and DESCRIPTION; Ko and Lange (dropped as `cph` in #33) restored; each
  derived file now names its source, licence and changes (comments only;
  builds). Also a stale comment in `radmix_cv.R`.
- New `paper_methods_text.md`: two paragraphs for the paper.

## 2026-11-17 (manica01): PR #36 merged; branches deleted

- Andrea merged PR #36 into `dev` (47f5780).
- Deleted on origin: `newspeed`, `supervised` (records in the tags
  `archive/newspeed`, `archive/supervised`), `newspeed_supervised`,
  `coverage-speed`; locally on manica01 also `newspeed_port`, `adamixture`,
  `rewrite`, `openadmixture`, `rclump`, and the worktrees of the old
  branches. Kept: `articles_revised`, `hgdp_plots`, `windows` (not merged).
- `tgc-desktop-checks.timer` disabled.

## 2026-11-16 (manica01): full HGDP benchmark in the article

- `admixture_benchmark.Rmd` now reads `data/admixture_benchmark` (full
  HGDP, 1043 x 130,943, K = 5; desktop P-cores, `newspeed_port` 2e2f5b4 at
  -O2), and its provisional note is gone (5eb3c73). The runs took 2.4 h.
- `gt_radmix()` from the command line: 39.8 s on 1 core, 6.6 s on 8
  (ADMIXTURE 326 / 360 s, fastmixture 455 / 69 s, ADAMIXTURE 75 / 10.7 s,
  OpenADMIXTURE.jl 220 / 50 s). ADMIXTURE ends 9.6 below the others in
  log-likelihood; the article's sentence on local optima now says so.
- ADMIXTURE's `-j` does not speed it up here (as in the quick run).

## 2026-11-16 (manica01): the port is faster than `dev`

- Whole fits, `dev` (a20eb08) against `newspeed_port` (2e2f5b4), desktop,
  idle, 3 alternated rounds: 0.78-0.85 of `dev`'s time on EUR K = 3-4 and
  HGDP K = 5-10 at 1 and 8 threads, except HGDP K = 10 on 1 thread (0.92:
  `dev` gained there with #35). Same iterations, log-likelihoods within
  7.4e-16, thread-invariant. Table in `radmix_performance.md`; data in
  `radmix_port_timing_2026-11/`.
- The full HGDP benchmark of `admixture_benchmark.Rmd` runs next, unattended
  (systemd `tgc-port-bench`, `~/tgc_bench_article/status.md`).

## 2026-11-16 (laptopX1): PR #37 closed; plan before the merge of #36

- PR #37 closed as superseded by #36 (its test change is in #36 unchanged;
  covr keeps `-O0`, with `test_radmix_vs_admixture.R` skipped under covr).
- **Decided (Andrea):** before merging #36, time `newspeed_port` against
  `dev` on the desktop to confirm the speed gain, then run the full HGDP
  benchmark of `admixture_benchmark.Rmd` (which shows only the `quick` run
  so far) and update the article on `newspeed_port`. The old P3 timing
  item (port against `newspeed`) is replaced by the first.

## 2026-11-16 (manica01): newspeed_port verified on the desktop

- `R CMD check --as-cran` of 3f2a8e4: 0 ERROR, 0 WARNING, 5 known NOTEs;
  tests FAIL 0, PASS 3517, CPU/elapsed 1.07.
- `NOT_CRAN=true R CMD check`: 1 NOTE (`.git`, an artefact of building
  from a worktree); FAIL 0, SKIP 1, PASS 3832, with the ADMIXTURE
  comparison and fastmixture.
- Not covered: `test_gt_adamixture.R` (python `adamixture` not installed on
  manica01). Results: comment on PR #36. PR #36 now waits only for
  Andrea's review.

## 2026-11-16 (laptopX1): why covr was slow; PR #37 folded into #36

- **Timing** (overnight, one run per set-up, laptop): under covr,
  `test_radmix_vs_admixture.R` (4-thread fits of EUR_subset) took 45 of the
  50 min of the suite, against 1 min at -O2. A plain `-O0 --coverage` build
  is as slow, plain `-O0` takes 9 min: gcov's counters in the threads, not
  covr nor -O0. `-O1` gains nothing. `OMP_NUM_THREADS` does not limit the
  tests (they pass `n_cores`). Details: `covr_timing_2026-11.md`.
- **Decided (Andrea):** skip that file under covr (850185b): the suite takes
  5 min under covr, coverage unchanged (96.83%); keep covr's default flags
  (PR #37's `-O1` dropped); tests run on CRAN use at most 2 threads (the
  fold tests of `test_radmix_cv.R`, 99eef49); PR #37's test changes
  cherry-picked onto `newspeed_port` (7a3c206, 850185b), so that one branch
  is verified on the desktop; #37 to be closed.
- **CI** on 0a6848e: every check passes; the coverage job takes 6.5 min
  (was 25-45), coverage 96.83%.

## 2026-11-15 (laptopX1): test-coverage fixed; out-of-bounds index in #36

- **Codecov:** `test-coverage` on `dev` failed only at the upload: Codecov
  answered "This repository is deactivated" from 12 Nov (last upload 5 Nov),
  although its API said active. Not the secret (unchanged since June 2025,
  accepted) nor the org (tidypopgen uploaded on 14 Nov). Andrea deactivated
  and reactivated the repository on Codecov; uploads work again (PR #37,
  20:48). PR runs had hidden it: `fail_ci_if_error` is off on PRs.
- **Bug in #36:** the coverage job aborted in `test_radmix_cv.R`: libstdc++
  bounds checks (on by default at -O0 with GCC 15, as under covr) caught
  `&q_sup[K * (i + 1)]` in `set_fixed()` (`radmix_fit.cpp`), one past the
  end for a last individual that is a reference. Fixed in 3192e45
  (`std::copy_n`); full suite clean with `-D_GLIBCXX_ASSERTIONS` under gdb.
- **PR #37** (`coverage-speed`, open): covr at `-O1` plus
  `-D_GLIBCXX_ASSERTIONS`; and `test_gt_adamixture.R` compares the thread
  variables before and after the fit (python's copy of the environment
  restored an `OMP_NUM_THREADS` that R had unset). `OMP_NUM_THREADS: 1` was
  tried and removed (Andrea chose A+B).
- **Result:** CI test times under covr are too noisy to judge `-O1`: the same
  `-O0` setup took 45.6 and 33.8 min, `-O1` 25.4 and 44.3 min; R CMD check
  runs the same tests in 4.2 min. Coverage 96.77% at `-O1`, 96.76% at `-O0`.

## 2026-11-15 (laptopX1): newspeed and supervised ported onto dev

- **Done:** rebuilt the package commits of `newspeed` (fused pass, CV bit
  mask, AVX2 inlining, cheaper log-likelihood, tile settings,
  `cv_max_iter`/`cv_rtol`) and `supervised` (supervised analyses, ADMIXTURE
  fixture, docs) on top of #35, with the resolutions of the merge
  `newspeed_supervised`; the kernels keep the GL template parameter, and
  `slice_gl()` reads the free individuals of a supervised fit through their
  index, so both features also work on genotype likelihoods. Then the two
  decided defaults, and the notes.
- **Where:** PR #36 (open), branch `newspeed_port`, commits 32aeab2 to
  the commit of this entry, on `dev` a20eb08; tags `archive/newspeed` (ec96315) and `archive/supervised`
  (fc43f9b) keep the old branches' records.
- **Result:** the full test suite passes on an -O2 build (3918
  expectations, none skipped, with the ADMIXTURE, OpenADMIXTURE.jl and
  NGSadmix fixtures); the new likelihood tests fail without the indexed
  read. Not timed yet (desktop item in `handover.md`).
- **Decisions:** (Andrea) newspeed defaults as recommended: `cv_rtol` 1e-7
  (implemented as `max(rtol, 1e-7)`), no subset start, FMA rounding
  accepted, inlined log-likelihood, `MIN_Q_TASKS` 128, automatic bit mask;
  the supervised API as built. (Claude) supervised analyses allowed for
  genotype likelihoods, for Andrea's review in the PR.
- **Notes changed:** `radmix_performance.md`, `admixture_implementation.md`;
  new `README.md`, `summary.md`, `handover.md`, `progress.md`.

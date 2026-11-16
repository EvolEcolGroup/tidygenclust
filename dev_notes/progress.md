# tidygenclust: progress

The log of major activity, newest first. Never rewrite old entries: add a
correction as a new entry. Earlier history (before 15 Nov 2026) is in the
comments of dramanica/claude_jobs#3 (radmix) and #12 (genotype
likelihoods), and in `radmix_performance.md`.

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

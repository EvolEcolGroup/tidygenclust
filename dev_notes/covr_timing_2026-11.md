# Why covr was slow: the timing of November 2026

The `test-coverage` job took 25-45 min in CI against 4.2 min for the same
tests in R CMD check. On the night of 15-16 Nov 2026 each test file was
timed under covr and in plain builds on laptopX1, to find where the time
went. Answer: almost all of it went to one file,
`tests/testthat/test_radmix_vs_admixture.R`, whose fits of EUR_subset run on
4 threads. The slowdown comes from gcov's line counters, not from covr's R
instrumentation or from compiling at -O0. That file is now skipped under
covr (`skip_on_covr()`), and covr keeps its default flags. Raw data and
scripts: `dev_notes/covr_timing/`.

## Design decisions

| decision | alternatives considered | why | date | decided by |
|---|---|---|---|---|
| Skip `test_radmix_vs_admixture.R` under covr (`skip_on_covr()`, 850185b on `newspeed_port`) | run it on 1 thread under covr; `-fprofile-update=single` (untested) | 45 of the 50 min of the suite under covr; coverage unchanged without it (96.83%), as the other tests run the same code | 2026-11-16 | Andrea |
| Keep covr's default flags (`-O0 --coverage`) | `-O1` plus `-D_GLIBCXX_ASSERTIONS` (PR #37, `coverage-speed`); `OMP_NUM_THREADS: 1` | `-O1` is no faster (4.9 vs 5.0 min with the skip; 52 vs 45 min on the file without it); GCC 15 on the CI runners turns the libstdc++ bounds checks on at -O0, which found the bug of #36 | 2026-11-16 | Andrea, on Claude's proposal |
| Tests that run on CRAN use at most 2 threads (fold tests of `test_radmix_cv.R` moved from 3 to 2 threads, 99eef49) | skip them on CRAN | CRAN's policy; 2 threads still test that fold fits do not depend on the number of threads | 2026-11-16 | Andrea |
| PR #37's test changes moved onto `newspeed_port` (#36), and #37 to be closed | merge #37 separately | one branch to verify on the desktop | 2026-11-16 | Andrea |

## How it works

- `OMP_NUM_THREADS` does not limit these tests: they pass `n_cores`
  explicitly, which radmix's kernels use in `num_threads()`. A covr run at
  `OMP_NUM_THREADS=1` still used 3.6 cores.
- The plain `-O0 --coverage` build (gcov counters, no covr) is as slow as
  covr on the file (45 min), while a plain `-O0` build takes 9 min and
  `-O2` 1 min. Every other file takes 4-6 min in total in every set-up, so
  the counters only hurt where several threads run hot loops on large
  data. The likely mechanism, not verified: the threads update the same
  counters (possibly atomically, `-fprofile-update=prefer-atomic`), and the
  cache lines bounce between cores.
- Heavy multithreaded tests should be skipped under covr, or run on one
  thread there.

## Numbers

Minutes, laptopX1 (4 cores), one run each (indicative only; other sessions
ran niced single-core jobs during the `-O0 --coverage` run). Source: commit
347189c (`newspeed_port`) with the `test_gt_adamixture.R` of 882ad09; the
"skip" runs add `skip_on_covr()` and the 2-thread fold tests.

| set-up | whole suite | `test_radmix_vs_admixture.R` | other files | CPU | coverage |
|---|---|---|---|---|---|
| plain `-O2` (release) | 5.0 | 1.1 | 3.9 | 8.7 | |
| plain `-O0` | 13.7 | 9.2 | 4.6 | 37.7 | |
| plain `-O0 --coverage` | 51.2 | 45.4 | 5.8 | 170.9 | |
| covr `-O0` (default) | 49.3 | 45.4 | 4.0 | 184.4 | 96.83% |
| covr `-O1` + bounds checks | 56.1 | 52.1 | 4.0 | 211.3 | 96.82% |
| covr `-O1`, `OMP_NUM_THREADS=1` | 53.3 | 49.3 | 4.0 | 190.4 | |
| covr `-O0`, skip | 4.1 (5.0 with the build) | skipped | 4.1 | 5.0 | 96.83% |
| covr `-O1`, skip | 4.1 (4.9 with the build) | skipped | 4.1 | 5.0 | 96.82% |

In CI, before the change: covr `-O0` took 45.6 and 33.8 min on `dev`
(14-15 Nov), `-O1` 25.4 and 44.3 min (#37); R CMD check 4.2 min.

To reproduce (the scripts have the laptop's paths in `W`): install the
package into one library per set-up with `R_MAKEVARS_USER` set to
`Makevars_O0` or `Makevars_O0cov` (`-O2`: the default flags), then
`TGC_LIB=<lib> TGC_TESTS=<tests> TGC_OUT=<csv> Rscript run_files.R` for a
plain build, or `Rscript run_covr.R <label> O0|O1` under covr. The
unattended night is `overnight.sh` (a systemd user unit), the last two runs
`focus.sh`; `status.md` is its log.

## Open questions

- Would `-fprofile-update=single` in `covr.flags` make the multithreaded
  tests fast under covr, so that the file need not be skipped? One covr
  run of that file would tell; not needed while the skip loses no
  coverage.
- The CI time after the change: expected about 5-10 min (the first run of
  #36 after the push will tell).

#!/bin/bash
# Runs the test suite of a build against its installed package (the one in
# $TGC_WORK/libs/<label>, compiled as for the benchmarks), with the tests of
# the same commit ($TGC_WORK/src/<label>/tests/testthat), as R CMD check
# does; unlike devtools::test(), this never compiles at -O0 or reuses stale
# object files, so it tests exactly the code that was timed.
#
#   bash run_tests.sh <label> [filter]
#
# e.g. bash run_tests.sh comb; bash run_tests.sh comb_nomv radmix
# Writes $TGC_WORK/logs/tests_<label>.log and prints the totals.
set -uo pipefail
REPO=$(git rev-parse --show-toplevel)
WORK=${TGC_WORK:-$(dirname "$REPO")/tgc_perf_work}
LABEL=$1; FILTER=${2:-}
LOG=$WORK/logs/tests_$LABEL${FILTER:+_$FILTER}.log
cd "$WORK/src/$LABEL/tests/testthat" || exit 1
R_LIBS="$WORK/libs/$LABEL${R_LIBS:+:$R_LIBS}" NOT_CRAN=true Rscript -e "
  library(testthat)
  suppressMessages(library(tidygenclust))
  cat('testing', dirname(getNamespaceInfo(asNamespace('tidygenclust'), 'path')), '\n')
  filter <- '$FILTER'
  res <- as.data.frame(test_dir('.', package = 'tidygenclust',
    load_package = 'installed', filter = if (nzchar(filter)) filter,
    reporter = ProgressReporter\$new(show_praise = FALSE), stop_on_failure = FALSE))
  cat('\nfailed:', sum(res\$failed), ' errors:', sum(res\$error), ' skipped:', sum(res\$skipped),
    ' passed:', sum(res\$passed), '\n')
" > "$LOG" 2>&1
grep -E "^testing|^failed:" "$LOG"

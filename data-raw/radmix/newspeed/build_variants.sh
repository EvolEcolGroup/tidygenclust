#!/bin/bash
# Builds variants of tidygenclust from commits, each in its own library, with
# R CMD INSTALL --preclean from a fresh `git archive` of the commit (never
# from a working tree, whose leftover .o files of load_all() would be -O0),
# and checks that the kernels were compiled with -O2 and with the AVX2
# clones:
#   $TGC_WORK/src/<label>        sources (git archive of the commit)
#   $TGC_WORK/libs/<label>       the installed package
#   $TGC_WORK/logs/<label>.log   the install log
#
# Run from the root of the repository (any branch: the commits only need to
# be known to git, e.g. after `git fetch`):
#   bash data-raw/radmix/newspeed/build_variants.sh <label>=<commit> ...
# e.g.
#   bash data-raw/radmix/newspeed/build_variants.sh base=b5999cc \
#     fused=newspeed~1 algo=newspeed-algo
# Environment: TGC_WORK (default ../tgc_newspeed_work), MAKEFLAGS (e.g.
# -j8; the default is R's), PIN (CPUs for the compiler, e.g. PIN=0-7: taskset
# -c $PIN), EXTRA_FLAGS (compiler flags for every build, through a private
# R_MAKEVARS_USER file, so that ~/.R/Makevars is never used).
set -euo pipefail
REPO=$(git rev-parse --show-toplevel)
WORK=${TGC_WORK:-$(dirname "$REPO")/tgc_newspeed_work}
mkdir -p "$WORK"/src "$WORK"/libs "$WORK"/logs
[ $# -ge 1 ] || { echo "usage: build_variants.sh <label>=<commit> ..."; exit 1; }
PRE=""; [ -n "${PIN:-}" ] && PRE="taskset -c $PIN"

build() {
  local label=${1%%=*} commit=${1#*=} src lib log mk so producer
  [ "$label" != "$1" ] || { echo "argument '$1' is not <label>=<commit>"; exit 1; }
  src=$WORK/src/$label; lib=$WORK/libs/$label; log=$WORK/logs/$label.log
  mk=$WORK/logs/$label.mk
  echo "== $label: commit $(git -C "$REPO" rev-parse --short "$commit") ($(git -C "$REPO" log -1 --format=%s "$commit"))"
  rm -rf "$src" "$lib"; mkdir -p "$src" "$lib"
  git -C "$REPO" archive "$commit" | tar -x -C "$src"
  : > "$mk"
  if [ -n "${EXTRA_FLAGS:-}" ]; then
    for v in CXXFLAGS CXX11FLAGS CXX14FLAGS CXX17FLAGS CXX20FLAGS CXX23FLAGS; do
      echo "$v += $EXTRA_FLAGS" >> "$mk"
    done
  fi
  R_MAKEVARS_USER=$mk $PRE R CMD INSTALL --preclean --no-test-load -l "$lib" "$src" > "$log" 2>&1 \
    || { echo "   INSTALL FAILED, see $log"; return 1; }
  # the optimisation level of the kernels, from the debugging information of
  # the object file (the last -O of the command line wins)
  producer=$(readelf --debug-dump=info "$src/src/radmix_kernels.o" 2>/dev/null | grep -m1 DW_AT_producer || true)
  echo "   radmix_kernels.o: $(echo "$producer" | grep -oE ' -O[0-3s]' | tail -1)${EXTRA_FLAGS:+ (extra: $EXTRA_FLAGS)}"
  echo "$producer" | grep -oE ' -O[0-3s]' | tail -1 | grep -q -- '-O2\|-O3' \
    || { echo "   NOT OPTIMISED: do not time this build"; return 1; }
  so=$lib/tidygenclust/libs/tidygenclust.so
  echo "   AVX2 clones in the library: $(nm "$so" | grep -c 'arch_x86_64_v3' || true)"
  (cd "$REPO" && Rscript -e "suppressMessages(library(tidygenclust, lib.loc = '$lib')); cat('   loads from', dirname(getNamespaceInfo(asNamespace('tidygenclust'), 'path')), '\n')")
  # the commit of the build, for the summaries
  git -C "$REPO" rev-parse "$commit" > "$lib/COMMIT"
}

for v in "$@"; do build "$v"; done
echo "finished"

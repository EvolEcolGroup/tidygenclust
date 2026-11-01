#!/bin/bash
# Builds the variants of tidygenclust compared in the verification, each in
# its own library, from the commits of the SVD start and of the kernel
# speed-ups (patches from ../patches, compiler flags through a private
# R_MAKEVARS_USER file, so that ~/.R/Makevars is never used):
#   $TGC_WORK/src/<label>        sources (git archive of the commit, patched)
#   $TGC_WORK/libs/<label>       the installed package
#   $TGC_WORK/logs/<label>.log   the install log
#   $TGC_WORK/makevars/<label>.mk
#
# Run from the root of the repository (any branch that has the commits):
#   bash data-raw/radmix/perf_verification/scripts/build_variants.sh [label ...]
# Without labels, builds the default set (all but the optional ones). Each
# build is checked: the compile line of radmix_kernels.cpp is printed, and
# the number of x86-64-v3 / v4 clones in the shared library is counted.
set -euo pipefail
REPO=$(git rev-parse --show-toplevel)
KIT=$REPO/data-raw/radmix/perf_verification
WORK=${TGC_WORK:-$(dirname "$REPO")/tgc_perf_work}
# svd-start: "Add an SVD start to gt_radmix() (init = "svd")"
BASE=${TGC_BASE:-7baa46c}
# kernel-speedups: "Speed up the kernels of gt_radmix()"
COMB=${TGC_COMB:-017e176}
mkdir -p "$WORK"/src "$WORK"/libs "$WORK"/logs "$WORK"/makevars

# label | commit | patch | extra compiler flags | note
VARIANTS="
base|$BASE|||baseline: the SVD start, kernels as on dev
div|$BASE|divisions.patch||exp 1a: two reciprocals instead of four divisions
dec|$BASE|decode.patch||exp 2: decoding fast path (v2, also for CV views)
pairs|$BASE|pairs.patch||exp 4c+4d: two loci per pass, Q and P sides
mv|$BASE|multiversion.patch||exp 3: x86-64-v3 multiversioning prototype (incl. loglik)
hexp|$BASE|hessian_expansion_REJECTED.patch||exp 4a: ADAMIXTURE's Hessian expansion (rejected)
comb|$COMB|||the combination (kernel-speedups)
comb_nomv|$COMB||-DRADMIX_NO_MULTIVERSION|the combination without multiversioning
f_o3|$BASE||-O3|flags: -O3
f_avx2fma|$BASE||-mavx2 -mfma|flags: AVX2 + FMA
f_native|$BASE||-O3 -march=native|flags: -O3 -march=native
comb_native|$COMB||-O3 -march=native -DRADMIX_NO_MULTIVERSION|optional: the combination built for this CPU
comb_o3|$COMB||-O3|optional: the combination at -O3
comb_v4|$COMB|V4||optional: the combination with an extra x86-64-v4 (AVX-512) clone
"
DEFAULT="base div dec pairs mv hexp comb comb_nomv f_o3 f_avx2fma f_native"
LABELS=${*:-$DEFAULT}

build() {
  local label=$1 line commit patch flags note src lib mk log
  line=$(echo "$VARIANTS" | grep "^$label|") || { echo "unknown variant $label"; exit 1; }
  IFS='|' read -r _ commit patch flags note <<< "$line"
  src=$WORK/src/$label; lib=$WORK/libs/$label
  mk=$WORK/makevars/$label.mk; log=$WORK/logs/$label.log
  echo "== $label ($note): commit $commit${patch:+, $patch}${flags:+, flags $flags}"
  rm -rf "$src" "$lib"; mkdir -p "$src" "$lib"
  git -C "$REPO" archive "$commit" | tar -x -C "$src"
  if [ "$patch" = "V4" ]; then
    # one more clone, for AVX-512 (x86-64-v4), before the AVX2 one
    grep -q '"arch=x86-64-v3"' "$src/src/radmix_kernels.cpp"
    sed -i 's/"arch=x86-64-v3"/"arch=x86-64-v4", "arch=x86-64-v3"/' "$src/src/radmix_kernels.cpp"
  elif [ -n "$patch" ]; then
    patch -s -p1 -d "$src" --no-backup-if-mismatch < "$KIT/patches/$patch"
  fi
  : > "$mk"
  if [ -n "$flags" ]; then
    for v in CXXFLAGS CXX11FLAGS CXX14FLAGS CXX17FLAGS CXX20FLAGS CXX23FLAGS; do
      echo "$v += $flags" >> "$mk"
    done
  fi
  R_MAKEVARS_USER=$mk R CMD INSTALL --preclean --no-test-load -l "$lib" "$src" > "$log" 2>&1 \
    || { echo "   INSTALL FAILED, see $log"; return 1; }
  echo "   compile: $(grep -E '(g\+\+|clang\+\+).*radmix_kernels\.cpp' "$log" | head -1 | grep -oE ' -(O[0-3s]|march=[^ ]+|m(avx|fma)[^ ]*|D ?RADMIX[^ ]*)' | tr -d '\n')"
  local so=$lib/tidygenclust/libs/tidygenclust.so
  echo "   clones: v3 $(nm "$so" | grep -c 'arch_x86_64_v3' || true), v4 $(nm "$so" | grep -c 'arch_x86_64_v4' || true)"
  (cd "$REPO" && Rscript -e "suppressMessages(library(tidygenclust, lib.loc = '$lib')); cat('   loads from', dirname(getNamespaceInfo(asNamespace('tidygenclust'), 'path')), '\n')")
}

for label in $LABELS; do build "$label"; done
echo "finished"

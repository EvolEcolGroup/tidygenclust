# A check that the machine gives the article's times: HGDP at K = 5 on 1 and
# 8 threads, the build before streaming (base) against the article's build
# (art), 2 balanced rounds; the article (24 Nov) has radmix_cli 40.6 s on 1
# thread and 6.5 s on 8.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
set -e
cd "$KIT_REPO"
CASES=hgdp:5 THREADS="1 8" ROUNDS=2 \
  bash data-raw/radmix/newspeed/run_bench.sh "$RES/hgdp_check" base art
grep -E "^\| hgdp" "$RES/hgdp_check/summary.md"

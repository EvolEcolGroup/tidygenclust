# Phase 1, the gate: no slowdown on data in memory. The kit's balanced
# rounds (base, new alternated), whole fits and single passes, HGDP at K = 5
# and 10 (blocked Hessians), EUR at K = 3 and 4, on 1 and 8 threads.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
cd "$KIT_REPO"
CASES=hgdp:5,hgdp:10,eur:3,eur:4 THREADS="1 8" ROUNDS=3 \
  bash data-raw/radmix/newspeed/run_bench.sh "$RES/gate_fits" base new
MODE=kernels CASES=hgdp:5,hgdp:10,eur:3,eur:4 THREADS="1 8" ROUNDS=3 \
  bash data-raw/radmix/newspeed/run_bench.sh "$RES/gate_kernels" base new

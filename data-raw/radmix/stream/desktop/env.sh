# Environment of the desktop benchmarks of streaming (sourced by the phases
# and by hand): the work folder, the builds (base: ada196_article before the
# merge of streaming, 51b4e9f; new: ada196_article with it), tidypopgen
# 1.0.0.9015, the timing kit (tag archive/newspeed) and the scripts of
# data-raw/radmix/stream from the new build's sources.
W=$HOME/project_temp/claude/tidygenclust/stream
export TGC_WORK=$W/work
export TMPDIR=$W/rtmp
export R_LIBS=$W/lib_tpg
export TGC_HGDP_DIR=$HOME/git/tidygenclust/vignettes/articles/data/hgdp
KIT_REPO=$HOME/git/tidygenclust_newspeed
STREAM=$TGC_WORK/src/new/data-raw/radmix/stream
export BENCH_DATA=$W/data
export LIB_BASE=$TGC_WORK/libs/base:$W/lib_tpg
export LIB_NEW=$TGC_WORK/libs/new:$W/lib_tpg
export N_CORES=8 FILE_MB=1430
export OMP_PROC_BIND=close OMP_PLACES=cores
RES=$W/results

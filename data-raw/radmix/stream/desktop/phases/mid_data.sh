# Data about 5 GB (50,000 x 100,000), for the article's section on data
# larger than the memory, with the article's script.
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
set -e
export R_LIBS=$TGC_WORK/libs/art:$W/lib_tpg LARGE_DIR=$W/mid LARGE_CORES=16
export LARGE_N_INDIV=50000 LARGE_N_LOCI=100000
rm -rf "$LARGE_DIR"
cd "$TGC_WORK/src/art/vignettes/articles"
Rscript code/admixture_benchmark_large_data.R
ls -la "$LARGE_DIR"

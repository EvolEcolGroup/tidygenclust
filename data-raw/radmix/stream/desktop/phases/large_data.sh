# The article's data: 250,000 x 200,000 simulated genotypes (50 GB backing
# file, 12.5 GB bed), on 16 cores (the data do not depend on them).
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
set -e
export R_LIBS=$TGC_WORK/libs/art:$W/lib_tpg LARGE_DIR=$W/large LARGE_CORES=16
cd "$TGC_WORK/src/art/vignettes/articles"
df -h "$W"
Rscript code/admixture_benchmark_large_data.R
ls -la "$LARGE_DIR"

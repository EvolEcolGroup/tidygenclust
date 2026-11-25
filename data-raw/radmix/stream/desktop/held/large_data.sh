# The article's data: 220,000 x 200,000 simulated genotypes (44 GB backing
# file, 11 GB bed), on 16 cores (the data do not depend on them).
source $HOME/project_temp/claude/tidygenclust/stream/env.sh
export R_LIBS=$TGC_WORK/libs/new2:$W/lib_tpg LARGE_DIR=$W/large LARGE_CORES=16
cd "$TGC_WORK/src/new2/vignettes/articles"
df -h "$W"
Rscript code/admixture_benchmark_large_data.R
ls -la "$LARGE_DIR"

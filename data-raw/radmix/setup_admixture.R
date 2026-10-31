# Install ADMIXTURE, used by the vignette of gt_radmix() to compare
# gt_radmix() with the original program, and by
# data-raw/radmix/regenerate_admixture_fixtures.R to create the
# ADMIXTURE results stored for
# tests/testthat/test_radmix_vs_admixture.R (the tests themselves do
# not need ADMIXTURE).
#
# ADMIXTURE is installed by tgc_tools_install() in the conda environment
# "ctidygenclust", which tidypopgen::gt_admixture() finds automatically.
# Alternatively, put the admixture binary on the PATH.
#
# Run from the root of the package, with tidygenclust installed:
#   Rscript data-raw/radmix/setup_admixture.R

tidygenclust::tgc_tools_install()

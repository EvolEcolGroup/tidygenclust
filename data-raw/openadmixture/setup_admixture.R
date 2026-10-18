# Install ADMIXTURE, used by the OpenADMIXTURE vignette to compare
# gt_openadmixture() with the original program, and by
# data-raw/openadmixture/regenerate_admixture_fixtures.R to create the
# ADMIXTURE results stored for
# tests/testthat/test_openadmixture_vs_admixture.R (the tests themselves do
# not need ADMIXTURE).
#
# ADMIXTURE is installed by tgc_tools_install() in the conda environment
# "ctidygenclust", which tidypopgen::gt_admixture() finds automatically.
# Alternatively, put the admixture binary on the PATH.
#
# Run from the root of the package, with tidygenclust installed:
#   Rscript data-raw/openadmixture/setup_admixture.R

tidygenclust::tgc_tools_install()

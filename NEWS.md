# tidygenclust (development version)
* `gt_clumppling()` is now a native R and C++ reimplementation of clumppling
  2.0.10, previously developed as the separate package `rclump`: it gives the
  same results as the python implementation (which the test suite checks
  against stored results), without Python, and is much faster. Its input is
  always a `gt_admix` object: Q files can be read with
  `tidypopgen::read_q_files()`, so the arguments `input_path`, `input_format`
  and `extension` are replaced by `x`. New arguments `cd_res`, `test_comm`,
  `comm_min` and `comm_max` expose further options of clumppling, and
  `mode_stats` gives the statistics of each mode. The example data are now
  `gt_admix` objects (`inst/extdata/capeverde_admix.rds` and
  `chicken_admix.rds`).
* The conda environment is no longer needed by default, as the core
  functionality is native. `tgc_tools_install()` gains a `tools` argument to
  choose which optional tools to install (`fastmixture` and ADMIXTURE), and
  `tgc_fastmixture_install()` and `tgc_admixture_install()` add a tool to an
  existing environment. clumppling is no longer installed, and the package no
  longer warns at startup if the conda environment is missing.
* New `gt_openadmixture()`: a native C++ (OpenMP) reimplementation of
  OpenADMIXTURE.jl, previously developed as the separate package
  `openadmixturer`. It estimates ancestry proportions directly from the
  backing files of `gen_tibble` objects (no Python or Julia needed), supports
  diploid and pseudohaploid data, missing genotypes and ADMIXTURE-style
  cross-validation, and returns `gt_admix` objects. Its results are checked
  against OpenADMIXTURE.jl and ADMIXTURE in the test suite.
* Fixed the colours of the clusters in the plots of `gt_clumppling` objects
  with K >= 10: clusters were ordered as text ("1", "10", "2", ...), so the
  same aligned cluster was drawn in a different colour (and position) for
  different values of K.
* `gt_fastmixture()` now uses the individuals and loci selected in a
  `gen_tibble`: it read the whole backing file, so a filtered `gen_tibble`
  failed, and a reordered one gave the rows of Q in the order of the backing
  file while `id` followed the `gen_tibble`.

# tidygenclust 0.1.4
* Updated clumppling to v2.0.10, which is now installed in the same conda
  environment as fastmixture (`ctidygenclust`); the separate `cclumppling`
  environment is no longer needed.

# tidygenclust 0.1.3
* Fixed `gt_fastmixture()` wrapper for fastmixture 1.3.0: P matrices are now
  returned in the original locus order, and `gen_tibble` input is shuffled and
  normalised as in fastmixture, giving the same results as PLINK input.

# tidygenclust 0.1.2
* Updated version of fastmixture to 1.3.0
* Use pre-tested conda environments to avoid dependency hell.

# tidygenclust 0.1.1
* Updated versions of fastmixture and clumppling

# tidygenclust 0.1.0
* Initial public version

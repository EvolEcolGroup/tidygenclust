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
* `gt_fastmixture(threads = )` now limits the number of threads (with the
  python package threadpoolctl, which `tgc_fastmixture_install()` now
  installs): the environment variables used before had no effect once R had
  started OpenMP, and they stayed set in the R session.
* `gt_clumppling()` checks its options: a non-finite `cd_res` made community
  detection loop forever, and invalid logical options gave cryptic errors.
* `subset_gt_clumppling(k = )` selects the right values of K: it failed when
  `k` covered all of `K_range`, selected K = 20 to 29 together with K = 2 (and
  so on), could not select values of K with gaps (e.g. 17 and 19 in the
  chicken example), and silently returned an empty object for values of K that
  were not in `K_range`.
* `subset_gt_clumppling(indivs = )` checks that the individuals are between 1
  and N, rather than silently returning rows of NA.
* `gt_openadmixture()` with `fix_p = TRUE` can project fewer individuals than
  `k` (e.g. a single new sample) onto a reference panel.
* `gt_openadmixture()` matches the rows of `init_p` to the loci by their
  names, when it has row names (as its P matrices do): a P with the loci in a
  different order was silently accepted, giving wrong results, and a P with
  additional loci was rejected.
* `gt_openadmixture(crossval = TRUE)` no longer changes the estimates of the
  following runs when `seed` is not given.
* `tgc_tools_install(reset = TRUE)` with some of the tools no longer removes
  the others: the tools that were installed in `ctidygenclust` are installed
  again, and on macOS resetting ADMIXTURE alone leaves `ctidygenclust`
  untouched.
* Installing fastmixture no longer replaces the numpy pinned by the conda
  environment with the latest version from pypi.
* On macOS, `tgc_admixture_install()` uses the conda binary of reticulate
  (rather than any `conda` on the PATH, which may be missing, e.g. in RStudio)
  and quotes the path of the yaml file.
* On macOS, a single path given to `tgc_tools_install(conda_yaml = )` is only
  used for `ctidygenclust`; it was also used to create the environment of
  ADMIXTURE.
* Plots of `gt_clumppling` objects subset by individuals draw the group
  boundaries and labels (and the number of runs of each mode) in the right
  place.
* The `group` argument of the plots of `gt_clumppling` objects can be a
  factor.
* `autoplot(type = "modes_within_k")` checks that `k` is a single value of
  `K_range`: other values gave an empty plot that failed when printed.
* `autoplot(type = "modes")` no longer warns for a `gt_clumppling` object with
  a single value of K.
* The vignette and the documentation of `autoplot()` no longer say that the
  edges of the 'modes' plot are labelled with their cost: they are coloured by
  it.
* The Newton update of Q in `gt_openadmixture()` stops cycling at vertex
  solutions, which made it up to about 20% slower near convergence (results
  are unchanged).
* `gt_openadmixture()` checks all its arguments before fitting: `fix_q = NA`
  was taken as TRUE, `rtol = NA` silently prevented convergence, non-integer
  values of `k`, `n_runs` and `n_folds` were truncated, invalid `max_iter`,
  `n_qn` and `n_cores` gave a generic error, and invalid seeds were only
  detected after the first fits.
* `gt_openadmixture()` on a `gen_tibble` grouped by more than one variable
  fails with a clear error before fitting, rather than after all the fits.
* `gt_fastmixture()` checks its arguments before starting python, with error
  messages naming the right arguments, and a single `seed` (such as the
  default) can be used with `n_runs > 1`: it is incremented for each run.

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

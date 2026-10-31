# tidygenclust (development version)
* New `gt_radmix()`: ancestry estimation with the ADMIXTURE method, using a
  modified version of the algorithm of OpenADMIXTURE (Ko et al. 2023).
  Its results are checked
  against OpenADMIXTURE.jl and ADMIXTURE in the test suite.
* `gt_clumppling()` is now a native R and C++ port, with modifications, of
  clumppling 2.0.10 (Liu, Kopelman & Rosenberg 2024, by Xiran Liu).
* The conda environment is no longer needed by default, as the core
  functionality is native. `tgc_tools_install()` gains a `tools` argument to
  choose which optional tools to install (`fastmixture` and ADMIXTURE), and
  `tgc_fastmixture_install()` and `tgc_admixture_install()` add a tool to an
  existing environment.
* Fixed the colours of the clusters in the plots of `gt_clumppling` objects
  with K >= 10: clusters were ordered as text ("1", "10", "2", ...), so the
  same aligned cluster was drawn in a different colour (and position) for
  different values of K.
* `gt_fastmixture()` now starts from the SVD-based initialisation of
  fastmixture by default (`random_init = FALSE`).
* `gt_fastmixture()` now uses the individuals and loci selected in a
  `gen_tibble`.
* `gt_fastmixture(threads = )` now limits the number of threads (with the
  python package threadpoolctl, which `tgc_fastmixture_install()` now
  installs): the environment variables used before had no effect once R had
  started OpenMP, and they stayed set in the R session.
* `subset_gt_clumppling(k = )` selects the right values of K: it failed when
  `k` covered all of `K_range`, selected K = 20 to 29 together with K = 2 (and
  so on), could not select values of K with gaps (e.g. 17 and 19 in the
  chicken example), and silently returned an empty object for values of K that
  were not in `K_range`.
* `subset_gt_clumppling(indivs = )` checks that the individuals are between 1
  and N, rather than silently returning rows of NA.
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
* `gt_fastmixture()` checks its arguments before starting python, with error
  messages naming the right arguments, and a single `seed` (such as the
  default) can be used with `n_runs > 1`: it is incremented for each run.
* On Intel Macs, fastmixture is compiled with the llvm installed by Homebrew
  (its location was hard-coded for Apple silicon).
* The `safety` argument of `gt_fastmixture()` is deprecated: it was never
  passed on to fastmixture, which turns on its safety updates by itself when
  an update is unstable.

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

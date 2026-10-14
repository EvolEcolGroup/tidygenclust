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

# fastmixture algorithm for population genetics clustering

This function implements the fastmixture algorithm for population
genetics clustering by calling the python module. If you use this
function, make sure that you cite the relevant paper by Santander,
Refoyo-Martínez, and Meisner (2024).

## Usage

``` r
gt_fastmixture(
  x,
  k,
  n_runs = 1,
  threads = 1,
  seed = 42,
  iter = 1000,
  tole = 1e-09,
  batches = 32,
  supervised = NULL,
  check = 5,
  power = 11,
  chunk = 8192,
  subsample = 0.7,
  min_subsample = 50000,
  max_subsample = 5e+05,
  als_iter = 1000,
  als_tole = 1e-04,
  no_freqs = TRUE,
  random_init = FALSE,
  safety = NULL,
  cv = NULL,
  cv_tole = 1e-07
)
```

## Arguments

- x:

  either a
  [`tidypopgen::gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.html),
  or the name of the binary plink file (without the .bed extension)

- k:

  the number of ancestral components (clusters), either a single value
  or a vector

- n_runs:

  the number of repeats for each k value

- threads:

  the number of threads to use (1)

- seed:

  the random seed (defaults to 42): a vector of length `n_runs`, with
  one seed per run (the same seeds are used for each value of `k`), or a
  single seed, used for the first run and incremented for the following
  ones

- iter:

  the maximum number of iterations (1000)

- tole:

  the tolerance in log-likelihood units between iterations (1e-9)

- batches:

  the number of maximum mini-batches (32)

- supervised:

  the name of the file with the supervised labels (NULL)

- check:

  the number of iterations to check for convergence (5)

- power:

  number of power iterations in randomised SVD (11)

- chunk:

  the number of SNPs in chunk operations (8192)

- subsample:

  Fraction of SNPs to subsample in SVD/ALS (0.7)

- min_subsample:

  Minimum number of SNPs to subsample in SVD/ALS (50000)

- max_subsample:

  Maximum number of SNPs to subsample in SVD/ALS (500000)

- als_iter:

  the maximum number of iterations in the ALS algorithm (1000)

- als_tole:

  the tolerance for the RMSE of P between iterations (1e-4)

- no_freqs:

  do not save P-matrix (TRUE)

- random_init:

  random initialisation of parameters (FALSE). By default, as in the
  fastmixture program, the parameters are initialised from a randomised
  SVD of the genotypes (followed by alternating least squares), which
  needs all the loci to be polymorphic; random starting values can take
  much longer to converge

- safety:

  deprecated and ignored: fastmixture turns on its safety updates by
  itself when an update is unstable

- cv:

  the number of cross-validation folds, greater than 1 (NULL, no
  cross-validation)

- cv_tole:

  the tolerance for the cross-validation error in scaled log-likelihood
  units (1e-7)

## Value

an object of class `gt_admix`. See
[`tidypopgen::gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_admixture.html)
for details.

## Details

This function returns a `gt_admix` object that can be plotted with
`autoplot`, and tidied with `tidy` methods from the `tidypopgen`
package. Cross-validation is not performed by default (`cv = NULL`); if
you want to include cross-validation you can set `cv` to the number of
folds, a value greater than 1 (ADMIXTURE performs 5-fold cv as default).

`fastmixture` is an optional tool, which is not needed by the rest of
`tidygenclust` (see
[`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix.md)
for a native implementation of the ADMIXTURE model). Install it with
`tgc_tools_install(tools = "fastmixture")`, or add it to an existing
`ctidygenclust` conda environment with
[`tgc_fastmixture_install()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_fastmixture_install.md)
(only on Linux at the moment, see
[`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_tools_install.md)).

## References

C. G. Santander, A. Refoyo Martinez, J. Meisner (2024) Faster
model-based estimation of ancestry proportions. bioRxiv
2024.07.08.602454; doi: https://doi.org/10.1101/2024.07.08.602454

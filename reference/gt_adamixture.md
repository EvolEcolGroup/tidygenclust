# ADAMIXTURE algorithm for population genetics clustering

This function runs ADAMIXTURE, a python package for biobank-scale
genetic clustering with the model of ADMIXTURE. It combines a spectral
initialisation (SVD and ALS) with accelerated block relaxation (BR-QN)
to reach high likelihoods with fast convergence. If you use this
function, make sure that you cite the paper by Saurina-i-Ricos, Mas
Montserrat and Ioannidis (2026).

## Usage

``` r
gt_adamixture(
  x,
  k,
  n_runs = 1,
  threads = 1,
  seed = 42,
  tol = 0.1,
  max_iter = 10000,
  no_freqs = TRUE,
  cv = NULL,
  verbose = TRUE
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

  the number of repeats for each k value (1)

- threads:

  the number of threads to use (1)

- seed:

  the random seed (defaults to 42): a vector of length `n_runs`, with
  one seed per run (the same seeds are used for each value of `k`), or a
  single seed, used for the first run and incremented for the following
  ones

- tol:

  the tolerance: the fit stops when the log-likelihood changes by less
  than `tol` between iterations (0.1)

- max_iter:

  the maximum number of iterations (10000)

- no_freqs:

  do not save P-matrix (TRUE)

- cv:

  the number of cross-validation folds, greater than 1 (NULL, no
  cross-validation)

- verbose:

  boolean; should the progress of ADAMIXTURE be printed (TRUE)?

## Value

an object of class `gt_admix`. See
[`tidypopgen::gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_admixture.html)
for details. As for
[`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/reference/gt_radmix.md),
the P matrices give the frequency of the alternate allele of each locus
(for PLINK files, of the first allele in the bim file).

## Details

This function returns a `gt_admix` object that can be plotted with
`autoplot`, and tidied with `tidy` methods from the `tidypopgen`
package. Cross-validation is not performed by default (`cv = NULL`); if
you want to include cross-validation you can set `cv` to the number of
folds, a value greater than 1 (e.g. 5 for 5-fold cross-validation).

`ADAMIXTURE` is an optional tool, which is not needed by the rest of
`tidygenclust` (see
[`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/reference/gt_radmix.md)
for a native implementation of the ADMIXTURE model). It is installed by
[`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_tools_install.md),
or it can be added to an existing `ctidygenclust` conda environment with
[`tgc_adamixture_install()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_adamixture_install.md)
(only on Linux at the moment, see
[`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_tools_install.md)).
Its CPU version, used here, does not need PyTorch.

## References

Saurina-i-Ricos, J., Mas Montserrat, D., & Ioannidis, A. G. (2026).
ADAMIXTURE: adaptive first-order optimization for biobank-scale genetic
clustering. Bioinformatics, 42(Supplement_1). doi:
https://doi.org/10.1093/bioinformatics/btag236

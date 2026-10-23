# tidygenclust <img src="./man/figures/logo.png" align="right" width="150"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/EvolEcolGroup/tidygenclust/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EvolEcolGroup/tidygenclust/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/EvolEcolGroup/tidygenclust/graph/badge.svg)](https://app.codecov.io/gh/EvolEcolGroup/tidygenclust)

<!-- badges: end -->

`tidygenclust` provides functions and methods to run genetic clustering
in R, and to align and compare multiple runs of the same or different K.
`tidygenclust` builds on `tidypopgen`, enhancing the grammar of population
genetics with a focus on genetic clustering.

Its core functionality is native, and needs neither Python nor Julia:

- `gt_openadmixture()`, a C++ reimplementation of
  [OpenADMIXTURE.jl](https://github.com/OpenMendel/OpenADMIXTURE.jl), to
  estimate ancestry proportions under the ADMIXTURE model;
- `gt_clumppling()`, an R and C++ reimplementation of the python package
  [clumppling](https://github.com/PopGenClustering/Clumppling), to align
  multiple clustering runs within and across K.

Optionally, `tidygenclust` can also run the python package
[fastmixture](https://github.com/Rosemeis/fastmixture) (`gt_fastmixture()`)
and the original ADMIXTURE program (`tidypopgen::gt_admixture()`), which it
installs in a conda environment on request.

## Installation

Currently installation is limited to Linux and OSX, but NOT Windows (you
can use the WSL). You can install the development version of
`tidygenclust` from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("EvolEcolGroup/tidygenclust")
```

`gt_openadmixture()` and `gt_clumppling()` work straight away. The
optional tools (`fastmixture` and ADMIXTURE) are installed with
`reticulate` in a conda environment, `ctidygenclust`, when you ask for them:

``` r
# both tools
tgc_tools_install()
# or only some of them
tgc_tools_install(tools = "fastmixture")
# add a tool to an existing environment
tgc_admixture_install()
# check what is installed
tgc_tools_version()
```

**Note for MacOS user**

To install `fastmixture` on MacOS, you need `llvm` installed via `brew` in
the first place. In `bash`, you will need:\

``` bash
brew install llvm libomp
```

## Overview of functionality

On overview of the functionality of `tidygenclust` and its integration
with `tidypopgen` is provided in the "Get Started"" vignette. This also
includes a detailed description of how to align multiple runs of
clustering and compare them.

## Clumppling

`gt_clumppling()` aligns the runs stored in a `gt_admix` object (created by
any of the clustering functions, or read from Q files with
`tidypopgen::read_q_files()`), detects modes among runs with the same K, and
aligns the modes across K. It reproduces clumppling 2.0.10, against which the
test suite checks it, and it is much faster (about 0.1 seconds, rather than
about 9 seconds with python, for the 200 runs of the capeverde example):

```r
capeverde <- readRDS(
  system.file("extdata", "capeverde_admix.rds", package = "tidygenclust")
)
clump_res <- gt_clumppling(capeverde, seed = 1)
autoplot(clump_res, type = "modes", group = capeverde$group)
```

See `data-raw/clumppling/README.md` for details of the port.

## OpenADMIXTURE

`gt_openadmixture()` is a native C++ (OpenMP) reimplementation of
[OpenADMIXTURE.jl](https://github.com/OpenMendel/OpenADMIXTURE.jl) (Ko et al.
2023), a fast algorithm for estimating individual ancestry proportions (Q)
and ancestral allele frequencies (P) under the ADMIXTURE model. It needs
neither Python nor Julia: genotypes are read directly from the backing file
of a `gen_tibble`, and results are returned as `gt_admix` objects, which work
with all the methods of `tidypopgen` and with `gt_clumppling()`.
Cross-validation, as in ADMIXTURE, is available with `crossval = TRUE`.

```r
adm <- gt_openadmixture(anole_gt, k = 2:4, n_runs = 3, seed = c(1, 2, 3),
                        crossval = TRUE, n_cores = 4)
autoplot(adm, type = "cv")
```

The implementation reproduces OpenADMIXTURE.jl to numerical precision and
agrees with ADMIXTURE 1.3.0, which the test suite checks against stored
results; on the 379 individuals and 54,051 SNPs of the `EUR_subset` data
used to test OpenADMIXTURE.jl, it is 2-3 times faster than ADMIXTURE on a
single thread, and about 10 times faster on 4 threads. See the
"OpenADMIXTURE" vignette, `?gt_openadmixture`, and
`data-raw/openadmixture/README.md` for details.

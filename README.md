
# tidygenclust

<!-- badges: start -->
[![R-CMD-check](https://github.com/EvolEcolGroup/tidygenclust/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EvolEcolGroup/tidygenclust/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`tidygenclust` provides functions and methods to run genetic clustering in R,
using the commonly used ADMIXTURE program as we as the python package `fastmixture`.
It also help to align and compare multiple runs of the same or differen K using
the functionalities of the python package `clumppling`. This package builds on
`tidypopgen`, enhancing the grammar of population genetics with a focus on 
genetic clustering.

## Installation

You can install the development version of tidygenclust from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("EvolEcolGroup/tidygenclust")
```

## Overview of functionality

On overview of the functionality of `tidygenclust` and its integration with `tidypopgen` is provided
in the overview vignette. A mode detailed description of how to align multiple runs of clustering
and compare them is provided in the vignette `clumppling overview`.


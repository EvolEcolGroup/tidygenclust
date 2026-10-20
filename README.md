# tidygenclust <img src="./man/figures/logo.png" align="right" width="150"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/EvolEcolGroup/tidygenclust/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EvolEcolGroup/tidygenclust/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/EvolEcolGroup/tidygenclust/graph/badge.svg)](https://app.codecov.io/gh/EvolEcolGroup/tidygenclust)

<!-- badges: end -->

`tidygenclust` provides functions and methods to run genetic clustering
in R, using the commonly used ADMIXTURE program as we as the python
package `fastmixture`. It also help to align and compare multiple runs
of the same or different K using the functionalities of the python
package `clumppling`. `tidygenclust` builds on `tidypopgen`, enhancing
the grammar of population genetics with a focus on genetic clustering.

## Installation

Fastmixture, ADAMIXTURE, and Clumppling run natively on Linux and
macOS. On Windows, `tidygenclust` runs them in WSL2 while R and RStudio remain on
Windows. Windows support is experimental and currently requires an x86-64
computer. Start with the [Windows setup guide](vignettes/windows.Rmd), which
covers WSL2, Ubuntu, and installation from scratch.
You can install the development version of
`tidygenclust` from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("EvolEcolGroup/tidygenclust")
```

The first time you use the package, install its analysis tools:

``` r
library(tidygenclust)
tgc_tools_install()
```

On Windows, install WSL2 and Ubuntu first, following the guide above. If your
distribution is not named `Ubuntu`, set `TIDYGENCLUST_WSL_DISTRO` to its name
in your Windows `.Renviron` file and restart R. The tool environments are
stored under `~/.local/share/tidygenclust` inside WSL. To include the optional
ADAMIXTURE algorithm, use `tgc_tools_install(adamixture = TRUE)`.
Classic ADMIXTURE remains a `tidypopgen` function and is not supported by this
Windows backend.

**Note for MacOS user**

We use `reticulate` to seamlessly integrate the required python packages
and dependencies
into R without additional installation steps for the user. 
However, for this to work correctly for **MacOS users**, you need 
`llvm` installed via `brew` in the first place. In `bash`, you will need:\

``` bash
brew install llvm libomp
```

## Overview of functionality

On overview of the functionality of `tidygenclust` and its integration
with `tidypopgen` is provided in the "Get Started"" vignette. This also
includes a detailed description of how to align multiple runs of
clustering and compare them.

# Install optional tools for `tidygenclust`

The core functionality of `tidygenclust` is native:
[`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix.md)
and
[`gt_clumppling()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_clumppling.md)
need neither Python nor conda. A few additional features rely on
external tools, which are optional and are only installed on request:

- `fastmixture`, a python package used by
  [`gt_fastmixture()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_fastmixture.md);

- `ADMIXTURE`, used by
  [`tidypopgen::gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_admixture.html);

- `ADAMIXTURE`, a python package used by
  [`gt_adamixture()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_adamixture.md).

## Usage

``` r
tgc_tools_install(
  tools = c("fastmixture", "admixture", "adamixture"),
  reset = FALSE,
  fastmixture_hash = "1f03a8980be4872341e5faca9d6f2dfc7ccbfe35",
  adamixture_version = "1.9.6",
  conda_yaml = "auto",
  ci_install = FALSE
)
```

## Arguments

- tools:

  a vector with the tools to install, any of "fastmixture", "admixture"
  and "adamixture". Defaults to all of them.

- reset:

  a boolean used to reset the virtual environment. Only set to TRUE if
  you have a broken virtual environment that you want to reset.

- fastmixture_hash:

  a string with the commit hash of the `fastmixture` version to install.
  Default is the latest tested version.

- adamixture_version:

  a string with the version of `ADAMIXTURE` to install from PyPI.
  Default is the latest tested version.

- conda_yaml:

  A vector of string with the path for the `yaml` files used for
  installation (one for `ctidygenclust` and, for osx only, one for
  admixture in x86 compatibility mode). The value of each string can
  either be "auto" (the default, select the latest available yaml that
  has been tested to work on a given operating system), NULL (let conda
  attempt to install the latest available packages from the internet),
  or a path to a specific conda yaml found in "inst/env_snapshots". For
  "auto" or NULL, it is also possible to provide a vector of length one,
  and the value will apply to all the conda environments (e.g.
  `conda_yaml="auto"` is equivalent to `conda_yaml = c("auto", "auto")`
  in osx); a single path is used for `ctidygenclust` (and "auto" for the
  environment of ADMIXTURE on osx). Tools that are not in `tools` are
  left out of the yaml files, as is `fastmixture`, which is always
  installed from github.

- ci_install:

  a boolean indicating if the installation is being run on continuous
  integration (CI) services. Default is FALSE. If TRUE, the function
  will look for the conda yaml file in the `inst/env_snapshots` folder
  of the package source directory, rather than in the installed package
  directory. This is useful when testing the package on CI services.

## Value

NULL (invisibly)

## Details

We use `reticulate` to install them in a conda environment,
`ctidygenclust`.

The optional tools are only available on Linux at the moment. On
Windows, they can be used through the Windows Subsystem for Linux (WSL).
On macOS, the OpenMP library of the python tools conflicts with the one
used by R packages such as `bigstatsr`, which crashes the R session.

`tgc_tools_install()` creates the conda environment with the chosen
`tools`. If the environment already exists, it only adds the tools that
are not installed yet (alternatively, use
[`tgc_fastmixture_install()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_fastmixture_install.md),
[`tgc_admixture_install()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_admixture_install.md)
and
[`tgc_adamixture_install()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_adamixture_install.md)
to add a tool to an existing environment).
[`tgc_tools_version()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_tools_version.md)
reports the tools that are installed.

For each tool, default to the latest version of these packages that has
been tested to work with `tidygenclust`. It is possible to provide a
more recent github commit for a specific tool, but this might lead to
incompatibilities and errors.

The conda environment is created from a conda yaml file with the
versions of all the packages that have been tested on each operating
system (see `conda_yaml`); `fastmixture` is then installed from its
github repository, and `ADAMIXTURE` from PyPI, with `pip`, run in the
environment with
[`reticulate::conda_run2()`](https://rstudio.github.io/reticulate/reference/conda_run2.html).
If the installation fails with the tested yaml, you can try
`conda_yaml = NULL` to let conda install the latest available packages.

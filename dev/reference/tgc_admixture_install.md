# Add `ADMIXTURE` to the conda environment of `tidygenclust`

`ADMIXTURE` (used by
[`tidypopgen::gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_admixture.html))
is an optional tool. This function adds it to an existing
`ctidygenclust` environment (created with
[`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_tools_install.md)).
Like the other optional tools, it is only available on Linux at the
moment (see
[`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_tools_install.md)).

## Usage

``` r
tgc_admixture_install(conda_yaml = "auto", ci_install = FALSE)
```

## Arguments

- conda_yaml:

  only used on OSX (not supported at the moment): the path of the `yaml`
  file used to create the `cadmixture86` environment, "auto" (the
  default, the latest tested yaml) or NULL (the latest version of
  `ADMIXTURE` in bioconda). See
  [`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_tools_install.md).

- ci_install:

  a boolean indicating if the installation is being run on continuous
  integration (CI) services. Default is FALSE. If TRUE, the function
  will look for the conda yaml file in the `inst/env_snapshots` folder
  of the package source directory, rather than in the installed package
  directory. This is useful when testing the package on CI services.

## Value

NULL (invisibly)

## Examples

``` r
if (FALSE) { # \dontrun{
# installs software in a conda environment (only on Linux at the moment)
tgc_admixture_install()
} # }
```

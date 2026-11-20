# Add `ADAMIXTURE` to the conda environment of `tidygenclust`

`ADAMIXTURE` (used by
[`gt_adamixture()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_adamixture.md))
is an optional tool. This function adds it to an existing
`ctidygenclust` environment (created with
[`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_tools_install.md)),
or installs another version, from PyPI. Its CPU version, used by
[`gt_adamixture()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_adamixture.md),
does not need PyTorch. Like the other optional tools, it is only
available on Linux at the moment (see
[`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/tgc_tools_install.md)).

## Usage

``` r
tgc_adamixture_install(adamixture_version = "1.8.9")
```

## Arguments

- adamixture_version:

  a string with the version of `ADAMIXTURE` to install from PyPI.
  Default is the latest tested version.

## Value

NULL (invisibly)

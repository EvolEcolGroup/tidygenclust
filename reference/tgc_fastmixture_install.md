# Add `fastmixture` to the conda environment of `tidygenclust`

`fastmixture` (used by
[`gt_fastmixture()`](https://evolecolgroup.github.io/tidygenclust/reference/gt_fastmixture.md))
is an optional tool. This function adds it to an existing
`ctidygenclust` environment (created with
[`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_tools_install.md)),
or reinstalls it. Like the other optional tools, it is only available on
Linux at the moment (see
[`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_tools_install.md)).

## Usage

``` r
tgc_fastmixture_install(
  fastmixture_hash = "1f03a8980be4872341e5faca9d6f2dfc7ccbfe35"
)
```

## Arguments

- fastmixture_hash:

  a string with the commit hash of the `fastmixture` version to install.
  Default is the latest tested version.

## Value

NULL (invisibly)

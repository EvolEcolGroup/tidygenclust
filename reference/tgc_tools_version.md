# Print version of the optional tools installed by `tidygenclust`

This function returns the version of the optional tools installed by
[`tgc_tools_install()`](https://evolecolgroup.github.io/tidygenclust/reference/tgc_tools_install.md):
`fastmixture`, `ADMIXTURE` and `ADAMIXTURE`. Tools that are not
installed have version NA.

## Usage

``` r
tgc_tools_version()
```

## Value

A list with the version of the optional tools installed by
`tidygenclust`

## Examples

``` r
# \donttest{
tgc_tools_version()
#> $fastmixture
#> [1] "v1.3.0"
#> 
#> $admixture
#> [1] "1.3.0"
#> 
#> $adamixture
#> [1] "1.9.6"
#> 
# }
```

# Tidy a `gt_clumppling` object

A `tidy` method to extract information from a
[`gt_clumppling`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_clumppling.md)
object, and return it as a tibble. It can extract:

- 'modes': all the modes detected by
  [`gt_clumppling()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_clumppling.md).
  The modes have label 'KxMy', where 'x' and 'y' represent the K value
  and the mode rank.

- 'major_modes': modes of rank 1 for each K.

- 'Q_modes', 'q_modes': a list of q matrices, one per mode, each tidied
  into a tibble

- 'Q_major_modes', 'q_major_modes': the same output as 'Q_modes' but
  subsetted to only the major modes.

## Usage

``` r
# S3 method for class 'gt_clumppling'
tidy(
  x,
  matrix = c("modes", "major_modes", "Q_modes", "q_modes", "Q_major_modes",
    "q_major_modes"),
  ...
)
```

## Arguments

- x:

  the
  [`gt_clumppling`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_clumppling.md)
  object

- matrix:

  a string defining the information to be extracted, one of: "modes",
  "major_modes", "Q_modes" (or "q_modes"), "Q_major_modes" (or
  "q_major_modes").

- ...:

  Additional arguments. Not used. Needed to match generic signature
  only.

## Value

a [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
of the information of interest (for "Q_modes" and "Q_major_modes", a
named list of tibbles, one per mode)

## Examples

``` r
# the runs of ADMIXTURE for K = 2 to 5 of the Cape Verde example of
# clumppling, aligned into modes
capeverde <- readRDS(
  system.file("extdata", "capeverde_admix.rds", package = "tidygenclust")
)
clump_res <- gt_clumppling(capeverde, seed = 1)

# the modes, with their K, their rank within K and their label
tidy(clump_res, matrix = "modes")
#> # A tibble: 7 × 3
#>       k     m label
#>   <dbl> <dbl> <chr>
#> 1     2     1 K2M1 
#> 2     3     1 K3M1 
#> 3     4     1 K4M1 
#> 4     4     2 K4M2 
#> 5     5     1 K5M1 
#> 6     5     2 K5M2 
#> 7     5     3 K5M3 

# the ancestry proportions of the major modes, one tibble per mode
q_major <- tidy(clump_res, matrix = "Q_major_modes")
names(q_major)
#> [1] "K2M1" "K3M1" "K4M1" "K5M1"
q_major[["K3M1"]]
#> # A tibble: 1,197 × 3
#>       id q     percentage
#>    <int> <chr>      <dbl>
#>  1     1 1      1.000    
#>  2     1 2      0.0000100
#>  3     1 3      0.0000100
#>  4     2 1      1.000    
#>  5     2 2      0.0000100
#>  6     2 3      0.0000100
#>  7     3 1      1.000    
#>  8     3 2      0.0000100
#>  9     3 3      0.0000100
#> 10     4 1      1.000    
#> # ℹ 1,187 more rows
```

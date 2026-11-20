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

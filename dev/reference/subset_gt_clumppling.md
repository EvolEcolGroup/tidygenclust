# Subset a `gt_clumppling` object

This function subsets `gt_clumppling` objects to a set of individuals or
a set of values of K. This is intended to create plot insets, or to
visualise a subset of individuals during data analysis. To understand
the modes within a subset of individuals in your data, you should subset
your `gt_admix` object and re-run `gt_clumppling`.

## Usage

``` r
subset_gt_clumppling(x, k = NULL, indivs = NULL)
```

## Arguments

- x:

  a gt_clumppling object

- k:

  a vector of values of K to subset to: consecutive values of `K_range`
  (e.g. `c(17, 19)` if `K_range` is 17, 19 and 21), in increasing order

- indivs:

  a vector of individual indices to keep

## Value

a gt_clumppling object subsetted to the individuals specified

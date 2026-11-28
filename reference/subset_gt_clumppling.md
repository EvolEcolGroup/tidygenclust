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

## Examples

``` r
# the runs of ADMIXTURE for K = 2 to 5 of the Cape Verde example of
# clumppling, aligned into modes
capeverde <- readRDS(
  system.file("extdata", "capeverde_admix.rds", package = "tidygenclust")
)
clump_res <- gt_clumppling(capeverde, seed = 1)

# the modes with K = 3 and 4
clump_k34 <- subset_gt_clumppling(clump_res, k = c(3, 4))
autoplot(clump_k34, type = "modes")


# the Cape Verdean individuals only, e.g. for an inset (their modes are
# still those of all the individuals, and autoplot() warns about it)
cape_verdeans <- which(capeverde_pops == "Cape Verdean")
clump_cv <- subset_gt_clumppling(clump_res, indivs = cape_verdeans)
clump_cv$N
#> [1] 44
```

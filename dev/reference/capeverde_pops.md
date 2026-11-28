# Population membership for the Capeverde dataset

Membership to populations for a dataset of 399 individuals, including 44
Cape Verdean. These are the populations of the Cape Verde example of
clumppling (Liu, Kopelman and Rosenberg 2024), whose runs of ADMIXTURE
are in `extdata/capeverde_admix.rds`.

## Usage

``` r
capeverde_pops
```

## Format

A vector of length 399

## Source

The example data of clumppling
(<https://github.com/PopGenClustering/Clumppling>), from the Cape Verde
data of Verdu P, Jewett EM, Pemberton TJ, Rosenberg NA, Baptista M
(2017) Parallel trajectories of genetic and linguistic admixture in a
genetically admixed creole population. Current Biology, 27(16),
2529-2535.

## Examples

``` r
table(capeverde_pops)
#> capeverde_pops
#>      British Cape Verdean       French      Gambian      Iberian     Mandenka 
#>           89           44           28          109          107           22 
```

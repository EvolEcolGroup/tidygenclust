# Run clumppling

This function aligns multiple runs of a clustering algorithm (e.g.
ADMIXTURE or STRUCTURE) with the method of clumppling (Liu, Kopelman and
Rosenberg 2024). `gt_clumppling()` is a native R and C++ port, with
modifications to enhance performance, of the python package clumppling
(version 2.0.10, by Xiran Liu;
<https://github.com/PopGenClustering/Clumppling>), and gives the same
results without the need for Python. The main modifications are that the
alignment of each pair of runs, an integer linear program in clumppling,
is solved exactly as an equivalent assignment problem (with an
implementation of the Hungarian algorithm adapted from cp-algorithms),
and that community detection uses ports of the implementations that
clumppling uses through cdlib (python-louvain, in C++, and
markov_clustering, in R) or the `igraph` package. The licences of the
ported code are in the file `COPYRIGHTS` of the package
(`system.file("COPYRIGHTS", package = "tidygenclust")`). If you use this
function, make sure that you cite the paper on clumppling by Liu,
Kopelman, and Rosenberg (2024).

## Usage

``` r
gt_clumppling(
  x,
  use_rep = TRUE,
  merge = TRUE,
  cd_method = "louvain",
  use_best_pair = TRUE,
  seed = NULL,
  output_path = NULL,
  cd_res = 1,
  test_comm = FALSE,
  comm_min = 1e-06,
  comm_max = 0.01
)
```

## Arguments

- x:

  a `gt_admix` object with at least two Q matrices

- use_rep:

  boolean, whether to use representative modes (alternative: average),
  defaults to TRUE

- merge:

  boolean, whether to merge two clusters when aligning K+1 to K,
  defaults to TRUE

- cd_method:

  the community detection method to use, one of 'louvain' (default),
  'leiden', 'infomap', 'markov_clustering', 'label_propagation',
  'walktrap', 'custom'. 'leiden', 'infomap', 'label_propagation' and
  'walktrap' require the 'igraph' package; as in clumppling (which calls
  them through 'cdlib'), only 'louvain', 'infomap' and
  'markov_clustering' use the similarity between replicates as edge
  weights, while 'leiden', 'label_propagation' and 'walktrap' only use
  the presence of edges (all pairs of replicates are connected, except
  for the most dissimilar pair). 'custom' reproduces the example custom
  method of clumppling (Markov clustering with a pruning threshold of
  0.2); alternatively, `cd_method` can be a function that takes the
  adjacency matrix between replicates (similarities rescaled to \[0,
  1\], with ones on the diagonal) and returns a vector of community
  labels.

- use_best_pair:

  boolean, whether to use best pair as anchor for across-K alignment
  (alternative: major), defaults to TRUE

- seed:

  (optional) an integer used to set the random seed, so that community
  detection (e.g. 'louvain') is reproducible. Defaults to NULL (no seed
  is set).

- output_path:

  deprecated and ignored. Results are returned as R objects rather than
  written to files.

- cd_res:

  the resolution of the 'louvain' community detection method, defaults
  to 1

- test_comm:

  boolean, whether to test for community structure (with a Tracy-Widom
  test) before running community detection; if there is no structure,
  all replicates form a single mode. Defaults to FALSE.

- comm_min, comm_max:

  cost thresholds for community detection: if all pairwise costs between
  replicates are below `comm_min`, they form a single mode; if they are
  all above `comm_max`, each replicate is a mode. Defaults to 1e-6 and
  1e-2.

## Value

a list of class `gt_clumppling` containing:

- aligned_modes: a list of data.frames, each data.frame is a Q-matrix

- cost_acrossK: a named list of costs for each pairwise K alignment

- mode_replicates: a list of replicate indices for each mode

- K_range: vector of K values analyzed

- N: number of individuals

- mode_stats: a data.frame with the size, cost and performance (average
  similarity among replicates) of each mode

## Details

The example data (`capeverde_admix.rds` and `chicken_admix.rds` in
`extdata`) are the example datasets of clumppling: runs of ADMIXTURE on
the Cape Verde data of Verdu et al. (2017), and runs of STRUCTURE on the
chicken data of Rosenberg et al. (2001).

The algorithm has three steps:

1.  all pairs of replicates with the same K are aligned, finding the
    matching of clusters that minimises the squared distance between
    their memberships;

2.  replicates with the same K are grouped into modes by community
    detection on the network of alignment costs; each mode is summarised
    by a representative replicate (the one with the smallest cost to the
    others) or by the average of its replicates;

3.  modes are aligned across K (merging two clusters of the larger K
    when two values of K differ by one), and their clusters are
    reordered so that they are consistent across K.

The replicates are the Q matrices of a `gt_admix` object, as returned by
the clustering functions of 'tidypopgen' (e.g.
[`tidypopgen::gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_admixture.html))
and 'tidygenclust' (e.g.
[`gt_radmix()`](https://evolecolgroup.github.io/tidygenclust/reference/gt_radmix.md));
to align replicates stored as Q files, read them into a `gt_admix`
object with
[`tidypopgen::read_q_files()`](https://evolecolgroup.github.io/tidypopgen/reference/read_q_files.html).
Replicates with the same K are numbered in the order in which they
appear in the `gt_admix` object. The rows of each Q matrix must sum to 1
(to within 1e-6), so memberships that were rounded (e.g. by STRUCTURE)
should be rescaled first.

If you would like to generate an annotated autoplot from your
`gt_clumppling` object, ensure that all individuals from the same
population are adjacent to one another in the `gt_admix` object supplied
to `gt_clumppling` (see
[`tidypopgen::gt_admix_reorder_q()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_admix_reorder_q.html)).
Autoplot 'group' argument requires that all individuals from the same
group are adjacent.

The only random step is community detection. For the default 'louvain'
method, the C++ implementation follows python-louvain (used by
clumppling), but random numbers are drawn from R's generator: the modes
are the same as in python unless the community structure is ambiguous,
but the same seed does not give the same random choices in R and python.
If `seed` is given, the state of R's random number generator is restored
on exit. When a mode includes replicates with identical Q matrices, the
representative of the mode can also differ from the python
implementation (it is then an identical replicate), as their costs only
differ by rounding errors.

## References

X. Liu, N. M. Kopelman, N. A. Rosenberg (2024) Clumppling: cluster
matching and permutation program with integer linear programming.
Bioinformatics 40(1): btad751. doi:
https://doi.org/10.1093/bioinformatics/btad751

N. A. Rosenberg, T. Burke, K. Elo, M. W. Feldman, P. J. Freidlin, M. A.
M. Groenen, J. Hillel, A. Mäki-Tanila, M. Tixier-Boichard, A. Vignal, K.
Wimmers, S. Weigend (2001) Empirical evaluation of genetic clustering
methods using multilocus genotypes from 20 chicken breeds. Genetics
159(2): 699-713.

P. Verdu, E. M. Jewett, T. J. Pemberton, N. A. Rosenberg, M. Baptista
(2017) Parallel trajectories of genetic and linguistic admixture in a
genetically admixed creole population. Current Biology 27(16):
2529-2535.

## Examples

``` r
capeverde <- readRDS(
  system.file("extdata", "capeverde_admix.rds", package = "tidygenclust")
)
clump_res <- gt_clumppling(capeverde, seed = 1)
clump_res$mode_stats
#>   Mode Representative Size         Cost Performance
#> 1 K2M1       39_K2R39   50 8.043783e-14   0.9999998
#> 2 K3M1       65_K3R15   50 3.067899e-09   0.9999592
#> 3 K4M1      112_K4R12   39 7.341008e-04   0.9788892
#> 4 K4M2      138_K4R38   11 4.944349e-02   0.7918295
#> 5 K5M1      192_K5R42   18 5.074185e-02   0.7939366
#> 6 K5M2      167_K5R17   18 3.704965e-02   0.8188414
#> 7 K5M3      160_K5R10   14 2.550706e-02   0.8478811
```

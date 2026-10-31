# data-raw/clumppling

`gt_clumppling()` implements the method of clumppling (Liu, Kopelman and
Rosenberg 2024). Its R and C++ code is a port, with the modifications listed
below, of the python package
[clumppling](https://github.com/PopGenClustering/Clumppling) by Xiran Liu
(version 2.0.10, MIT licence), which previous versions of `tidygenclust` ran
through `reticulate`. The Louvain method and Markov clustering are ported
from python-louvain and markov_clustering, which clumppling uses through
cdlib; the licences of all the ported code are in `inst/COPYRIGHTS`. The scripts in this directory create
the example data used by `gt_clumppling()`, and the reference results of the
python implementation used by the tests. Run every script from the root of
the package. None of this is needed to use the package or to run its tests.

## How the port reproduces clumppling

The R code is a close port of clumppling 2.0.10 (`clumppling/core.py` and
friends) and of `py_rclumppling.py` (the python wrapper used by previous
versions of `tidygenclust`, kept here); function names follow the python ones
(e.g. `align_across_k()` becomes `clumppling_align_across_k()`).

- **Alignment** (`src/clumppling_assignment.cpp`). clumppling solves an
  integer linear program with GLPK for every pair of matrices. The same
  problem (each cluster of Q is matched to one cluster of P, each cluster of P
  receives at least one) is solved exactly as a linear assignment problem,
  with extra "free" columns when Q has more clusters than P, using the
  Hungarian algorithm.
- **Mode detection** (`src/clumppling_louvain.cpp`,
  `R/clumppling_community.R`). The default `louvain` method is a C++ port of
  python-louvain, which clumppling uses through cdlib. Markov clustering (and
  the `custom` example of clumppling) is a port of cdlib's implementation.
  `leiden`, `infomap`, `label_propagation` and `walktrap` use the `igraph`
  package, with the same choice of weighted or unweighted graph as cdlib. The
  Tracy-Widom test for community structure (`test_comm`) uses the critical
  values computed by clumppling.
- **Everything else** (validation of the Q matrices, mode statistics,
  alignment across K and reordering of clusters) is ported line by line.

On the 200 capeverde replicates (K = 2 to 5) the results are identical to
those of the python implementation (same modes, representatives and cluster
order, with numerical differences of about 1e-16), and the run takes about
0.1 seconds instead of about 9 seconds. The article
`vignettes/articles/clumppling_benchmark.Rmd` compares the times and the
results of the two implementations on the example data.

### Differences from the python implementation

- Q files are not read directly: they are read into a `gt_admix` object with
  `tidypopgen::read_q_files()`. Formats that `tidypopgen` does not read (e.g.
  '.indivq' or STRUCTURE files) have to be turned into a `gt_admix` object by
  the user (see `create_example_gt_admix.R`); unlike clumppling, individuals
  are not reordered by their labels, and the rows of every Q matrix must sum
  to 1 to within 1e-6 (clumppling allows 1e-2 for STRUCTURE files, whose
  memberships are rounded).
- Community detection draws random numbers from R's generator, so the same
  seed does not give the same random choices as in python. Modes are the same
  unless the community structure is ambiguous, in which case python's own
  results change with the seed (e.g. K = 17 and K = 19 in the `chicken`
  example). When `seed` is given, R's random number generator is restored on
  exit.
- When a mode includes replicates with identical Q matrices, their costs to
  the other replicates differ only by rounding errors, so the representative
  of the mode can be a different (but identical) replicate.
- Additional arguments expose clumppling options that the python wrapper
  fixed: `cd_res`, `test_comm`, `comm_min` and `comm_max`. `cd_method` can also
  be a function, and `leiden` works (it needs `leidenalg` in python).

## Example data: `inst/extdata/*_admix.rds`

`create_example_gt_admix.R` creates two `gt_admix` objects from the example
datasets of clumppling, also analysed in the paper that describes it, and
stored here as the original Q files (which previous versions of
`tidygenclust` distributed in `inst/extdata`):

- `capeverde_admix.rds`, from `capeverde.zip`: 200 ADMIXTURE runs (50 for
  each K from 2 to 5) for 399 individuals, with their ids and populations
  (the names of the populations are those of the `capeverde_pops` dataset),
  from the Cape Verde data of Verdu et al. (2017);
- `chicken_admix.rds`, from `chicken_gapk.zip`: 60 STRUCTURE runs (K = 17, 19
  and 21) for 600 chickens from 20 breeds, from the data of Rosenberg et
  al. (2001). STRUCTURE rounds memberships to 3 decimals, so rows are
  rescaled to sum to 1.

Files are added in alphabetical order, as the python reference sorts them, so
that replicates are numbered in the same way.

```bash
Rscript data-raw/clumppling/create_example_gt_admix.R
```

## Reference results: `tests/testthat/fixtures/clumppling_python_reference.rds`

Used by `tests/testthat/test_*clumppling*.R`. It holds the results of the
python implementation for: optimal alignments of pairs of Q matrices, community
detection on simulated graphs and on special cases, alignment within K and
community detection with each method on the capeverde data, and complete runs
on the capeverde data (with different options, and with shuffled individuals
and runs with K = 1), on the chicken data (read from the original STRUCTURE
files) and on the anole data (`tests/data/anole_q`, which has a single run for
each K). Large Q matrices are stored as fingerprints (column sums, sums of
squares and weighted sums). It also records the versions of the main python
packages.

To regenerate it:

1. Create a conda environment (`clumppling_reference`) with clumppling 2.0.10
   and the pinned versions of its dependencies (`requirements.txt`):

   ```bash
   bash data-raw/clumppling/setup_clumppling_env.sh
   ```

   Set `CLUMPPLING_REFERENCE_ENV` to use a differently named environment
   (for this script and the ones below). The `infomap` python package is not
   installed, so that cdlib uses the implementation of infomap in igraph, as
   the native port does.

2. Run the python implementation and save the results (about a minute):

   ```bash
   Rscript data-raw/clumppling/regenerate_python_fixtures.R
   ```

   This uses `py_clumppling.R`, which runs `py_rclumppling.py` (complete runs)
   and `clumppling_reference.py` (unit-level cases) through `reticulate`, and
   the helpers of the tests (`tests/testthat/helper-clumppling.R`), so that
   both implementations get exactly the same inputs.

## Validation on random scenarios

`validate_vs_python.R` compares `gt_clumppling()` with the python
implementation on 60 random scenarios: random subsets of the capeverde runs
(2 to 4 values of K, 1 to 50 runs for each, with shuffled individuals in half
of the scenarios), and simulated runs with two planted modes for each K (50
to 200 individuals, 2 to 6 clusters), with random options (`use_rep`, `merge`,
`use_best_pair`, `test_comm`, and all the community detection methods).

```bash
Rscript data-raw/clumppling/validate_vs_python.R
```

The results are in `validate_vs_python.csv`: in all 60 scenarios, the modes,
their statistics, the costs across K and the aligned modes are the same; in
two of them, a mode had several identical replicates, and the two
implementations chose different (identical) representatives.

## References

Liu X, Kopelman NM, Rosenberg NA (2024) Clumppling: cluster matching and
permutation program with integer linear programming. Bioinformatics, 40(1),
btad751.

Rosenberg NA, Burke T, Elo K, Feldman MW, Freidlin PJ, Groenen MAM, Hillel J,
Mäki-Tanila A, Tixier-Boichard M, Vignal A, Wimmers K, Weigend S (2001)
Empirical evaluation of genetic clustering methods using multilocus
genotypes from 20 chicken breeds. Genetics, 159(2), 699-713.

Verdu P, Jewett EM, Pemberton TJ, Rosenberg NA, Baptista M (2017) Parallel
trajectories of genetic and linguistic admixture in a genetically admixed
creole population. Current Biology, 27(16), 2529-2535.

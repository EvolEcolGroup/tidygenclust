# data-raw/openadmixture

Scripts used to validate `gt_openadmixture()` against the original
implementations, and to create the reference results (fixtures) used by the
test suite. Run every script from the root of the package. None of this is
needed to use the package or to run its tests.

## Test fixtures

Neither OpenADMIXTURE.jl nor ADMIXTURE is under active development, so their
results are computed once and stored in `tests/testthat/fixtures/`. The tests
compare `gt_openadmixture()` with these stored results, so they need neither
Julia nor ADMIXTURE. The fixtures only need regenerating if a comparison
changes (e.g. different data, values of K or seeds).

### OpenADMIXTURE.jl: `fixtures/openadmixture_julia_reference.rds`

Used by `tests/testthat/test_openadmixture_vs_julia.R`. On a small simulated dataset (150
individuals, 300 loci, K = 4, 3% missing genotypes), and starting from fixed
values of Q and P, it holds the results of OpenADMIXTURE.jl for: the
log-likelihood, one EM step, one round of Newton updates of Q and P, and the
full quasi-Newton algorithm (with the log-likelihood at each iteration). The
simulated genotypes and the starting values are stored too, so the test does
not depend on R's random number generator. It also records the versions of
Julia and OpenADMIXTURE.jl that were used.

To regenerate it:

1. Install Julia and OpenADMIXTURE.jl in a dedicated conda environment
   (`openadmixturer_julia`; Julia's packages are kept inside the
   environment, not in `~/.julia`):

   ```bash
   bash data-raw/openadmixture/setup_julia_env.sh
   ```

   OpenADMIXTURE.jl is pinned to the commit used for the port. Julia is
   pinned to 1.12.7: conda-forge's Julia 1.10.4 crashes when downloading
   packages on recent Linux kernels, and conda-forge has no 1.11 builds.

2. Run OpenADMIXTURE.jl and save the results:

   ```bash
   Rscript data-raw/openadmixture/regenerate_julia_fixtures.R
   ```

   This writes the data as a PLINK bed file and calls
   `data-raw/openadmixture/julia/openadmixture_reference.jl`, which uses the internal functions of
   OpenADMIXTURE.jl. Julia is run through `conda run`; set `OPENADMIXTURE_JULIA_ENV`
   to use a differently named conda environment, or `OPENADMIXTURE_JULIA` to the path
   of a julia executable whose default environment includes
   OpenADMIXTURE.jl.

Two details matter when comparing the two implementations, and are handled
by these scripts:

- SnpArrays (used by OpenADMIXTURE.jl) counts copies of the second allele of
  the bim file, which for a `gen_tibble` is `allele_ref`, whereas
  `gt_openadmixture()` counts `allele_alt`. P matrices are flipped (1 - P) when
  passed to and from Julia; Q is unaffected.
- OpenADMIXTURE.jl indexes its per-thread buffers with `threadid()`, which
  breaks with the interactive thread added by default since Julia 1.12, so
  Julia is started with `--threads=1,0`.

### ADMIXTURE: `fixtures/openadmixture_admixture_eur.rds`

Used by `tests/testthat/test_openadmixture_vs_admixture.R`. It holds the estimates of
ADMIXTURE 1.3.0 (seed 1) for K = 2 to 4 on the `EUR_subset` data of
SnpArrays.jl (379 individuals, 54,051 SNPs), the data used in the
OpenADMIXTURE.jl tests: Q and P (stored exactly, as integers in millionths,
since ADMIXTURE writes 6 decimal places), the log-likelihood, the 5-fold
cross-validation error (`--cv`) and the log of each run, together with MD5
checksums of the data files. The data
themselves are not stored: the test downloads them (and checks the
checksums), so it is skipped when offline.

To regenerate it:

1. Install ADMIXTURE in the conda environment created by
   `tgc_tools_install()` (or put the `admixture` binary on the PATH):

   ```bash
   Rscript data-raw/openadmixture/setup_admixture.R
   ```

2. Download the data:

   ```bash
   Rscript data-raw/openadmixture/download_eur_subset.R
   ```

3. Run ADMIXTURE and save the results (this takes about 30-40 minutes,
   mostly because of the cross-validation):

   ```bash
   Rscript data-raw/openadmixture/regenerate_admixture_fixtures.R
   ```

ADMIXTURE's folds are random and cannot be reproduced (its source code is
not available), so the test only expects the cross-validation errors of
`gt_openadmixture()` to agree closely with ADMIXTURE's, and to rank K the same
way.

ADMIXTURE reports the frequency of the second allele of the bim file
(`allele_ref`), so the test compares 1 - P with the P matrices of
`tidygenclust`. ADMIXTURE also writes P with loci sorted by chromosome and
position; this is the order of the input for `EUR_subset`, which the test
checks.

## Other scripts

These scripts use the `EUR_subset` data downloaded by
`data-raw/openadmixture/download_eur_subset.R` (into `data-raw/openadmixture/data/`, which is not under
version control), and need `tidygenclust` to be installed (timings from
`devtools::load_all()` are misleading, as it compiles without optimisation):

- `validate_eur.R`: fits K = 4 with 1 and 4 cores, reporting times and the
  log-likelihood (compare with -15,552,193.9 reported by OpenADMIXTURE.jl);
- `compare_admixture.R`: compares a run of `gt_openadmixture()` with a live
  run of ADMIXTURE (via `tidypopgen::gt_admixture()`);
- `profile_eur.R`: times the individual C++ kernels;
- `benchmark_admixture.R`: times `gt_openadmixture()` against ADMIXTURE
  (single-threaded, K = 2 to 4, about 5 minutes of ADMIXTURE runs), writing
  the results to `benchmark_admixture.csv`.

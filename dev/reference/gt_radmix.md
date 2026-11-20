# Estimate ancestry proportions with the ADMIXTURE method

This function estimates individual ancestry proportions (Q) and the
allele frequencies of the ancestral populations (P) with the method of
ADMIXTURE (Alexander et al. 2009), i.e. by maximum likelihood under the
ADMIXTURE model. It uses the optimisation strategy of ADMIXTURE: a few
EM iterations to prime a random start, then block-wise Newton updates of
Q and P, each solved as a constrained quadratic program, accelerated by
the quasi-Newton scheme of Zhou, Alexander and Lange (2011). It also
offers the cross-validation of ADMIXTURE to choose K (Alexander and
Lange 2011), a start from a truncated singular value decomposition (SVD)
of the genotypes, supervised analyses with individuals of known
ancestry, and support for pseudohaploid data. With genotype likelihoods
(given by another package, see
[`gt_radmix_data()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix_data.md)),
it fits the model of NGSadmix (Skotte et al. 2013), the ADMIXTURE model
for low-depth sequencing data, with the same algorithm. If you use this
function, please cite Alexander et al. (2009) for the method, and, if
you use cross-validation, Alexander and Lange (2011); with genotype
likelihoods, also cite Skotte et al. (2013).

## Usage

``` r
gt_radmix(
  x,
  k,
  n_runs = 1,
  seed = NULL,
  max_iter = 1000,
  rtol = 1e-08,
  n_em_iter = 5,
  n_qn = 3,
  init = c("random", "svd"),
  svd_power = 3,
  init_q = NULL,
  init_p = NULL,
  fix_q = FALSE,
  fix_p = FALSE,
  supervised = NULL,
  use_imputed = FALSE,
  crossval = FALSE,
  n_folds = 5,
  cv_max_iter = max_iter,
  cv_rtol = max(rtol, 1e-07),
  n_cores = 1,
  verbose = FALSE
)
```

## Arguments

- x:

  a
  [`tidypopgen::gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.html),
  or genotype likelihoods from another package (see "Genotype
  likelihoods" in the details), possibly grouped

- k:

  an integer vector giving the number(s) of ancestral populations

- n_runs:

  the number of runs for each value of `k` (defaults to 1)

- seed:

  the seed for the random number generator (defaults to `NULL`). Either
  a vector of length `n_runs` (the same seeds are used for each `k`) or
  of length `n_runs * length(k)`. The global random number generator
  state is restored on exit.

- max_iter:

  the maximum number of iterations (defaults to 1000)

- rtol:

  the convergence threshold on the relative change in log-likelihood
  between iterations (defaults to 1e-8). A threshold of 1e-7 takes
  10-20% less time, but can stop a few log-likelihood units short of
  convergence, with ancestry proportions up to about 0.05 away from the
  converged ones

- n_em_iter:

  the number of EM iterations used to prime the estimates (defaults to
  5, as in ADMIXTURE); not used with `init = "svd"`

- n_qn:

  the number of secant pairs used by the quasi-Newton acceleration
  (defaults to 3, as in ADMIXTURE)

- init:

  how the starting values of each run are chosen: `"random"` (the
  default), random values followed by `n_em_iter` EM iterations, as in
  ADMIXTURE; or `"svd"`, from a truncated SVD of the genotypes, which
  needs fewer iterations but is not suited to exploring the modes of the
  likelihood with several runs (see "Starting values from the SVD" in
  the details). `"svd"` cannot be combined with `init_q` or `init_p`.

- svd_power:

  the number of power iterations of the randomised SVD of genotype
  likelihoods used by `init = "svd"` (defaults to 3; see "Starting
  values from the SVD" in the details). Not used with genotypes, whose
  truncated SVD is computed by
  [`tidypopgen::gt_pca_partialSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_partialSVD.html)

- init_q:

  optional matrix (individuals x k) of starting ancestry proportions;
  only allowed with a single value of `k`

- init_p:

  optional matrix (loci x k) of starting allele frequencies (of the
  alternate allele, or of `allele2` for genotype likelihoods, as in the
  P matrices returned by this function; those of
  [`tidypopgen::gt_admixture()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_admixture.html)
  give the frequency of the reference allele, so use `1 - P`); only
  allowed with a single value of `k`. If it has row names (as the P
  matrices returned by this function), they are used to match its rows
  to the loci of `x` (it can then include other loci, which are
  ignored); otherwise its rows must be in the order of the loci of `x`

- fix_q, fix_p:

  booleans; if `TRUE`, Q (respectively P) is kept at the value given in
  `init_q` (respectively `init_p`)

- supervised:

  optional vector (character or factor) with one element per individual
  of `x`, for a supervised analysis (as ADMIXTURE's `--supervised`): the
  population of each reference individual, whose ancestry is then fixed,
  and `NA` for the individuals to estimate. `k` must be a single value,
  equal to the number of populations. See "Supervised analyses" in the
  details.

- use_imputed:

  boolean; if `TRUE`, imputed genotypes are used as if they were
  observed; otherwise they are treated as missing (the default). Only
  for genotypes.

- crossval:

  boolean; should the cross-validation error be computed (defaults to
  `FALSE`)? See the details. Only for genotypes.

- n_folds:

  the number of folds used for cross-validation (defaults to 5, as in
  ADMIXTURE)

- cv_max_iter, cv_rtol:

  the maximum number of iterations and the convergence threshold of the
  refit of each fold in cross-validation (default to `max_iter`, as for
  the full data, and to 1e-7, or `rtol` if it is larger). Smaller
  `cv_max_iter` or larger `cv_rtol` make cross-validation faster but
  bias its errors downwards (see "Cross-validation" in the details).

- n_cores:

  the number of cores to use (defaults to 1); the results do not depend
  on it

- verbose:

  boolean; should progress be printed (defaults to `FALSE`)? With
  `verbose = TRUE`, the log-likelihood is also computed (and logged)
  every 5 EM steps during the initial EM iterations; otherwise only
  their number is logged, as these values are not needed by the
  algorithm.

## Value

an object of class `gt_admix`, a list with elements:

- `k` the number of clusters of each run

- `Q` a list of
  [`tidypopgen::q_matrix()`](https://evolecolgroup.github.io/tidypopgen/reference/q_matrix.html)
  objects with the ancestry proportions

- `P` a list of matrices (loci x k) with the alternate allele
  frequencies (frequencies of `allele2` for genotype likelihoods)

- `log` a list with the log of each run

- `loglik` the log-likelihood of each run

- `converged` whether each run converged

- `cv` the cross-validation error of each run (if `crossval` is `TRUE`)

- `id` the id column of `x`

- `group` the group column of `x` (if `x` is grouped)

- `ref_pops` the reference population of each column of Q and P (only
  for a supervised analysis)

- `algorithm` "radmix"

## Details

Genotypes are read directly from the backing file of the `gen_tibble`
(respecting any filtering of individuals and selection of loci), so no
intermediate PLINK file is written; the same holds for genotype
likelihoods, which the package that stores them gives to `gt_radmix()`
through a method of
[`gt_radmix_data()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix_data.md).
Computation is in C++, multithreaded with OpenMP. The C++ code started
as a port of the Julia package
[OpenADMIXTURE.jl](https://github.com/OpenMendel/OpenADMIXTURE.jl) (Ko
et al. 2023), and has since been extensively rewritten; its solver of
the quadratic programs, its projections, its quasi-Newton step and the
structure of its fitting loop are still derived from OpenADMIXTURE.jl,
and the alternating least squares of `init = "svd"` are adapted from
[ADAMIXTURE](https://github.com/AI-sandbox/ADAMIXTURE) (see the file
`COPYRIGHTS` of the package,
`system.file("COPYRIGHTS", package = "tidygenclust")`).

Missing genotypes are ignored in the likelihood. Pseudohaploid
individuals (see
[`tidypopgen::gt_pseudohaploid()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pseudohaploid.html))
contribute a single allele per locus. Only diploid and pseudohaploid
`gen_tibble`s are supported.

`P` gives the frequency of the **alternate** allele (`allele_alt` in
[`tidypopgen::show_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/show_loci.html))
in each ancestral population (for genotype likelihoods, of `allele2`,
often the minor allele; see "Genotype likelihoods" below).

Supplying `init_p` together with `fix_p = TRUE` estimates the ancestry
of new individuals given the allele frequencies of a reference panel (a
"projection" analysis); the loci of `x` must then match the rows of
`init_p`.

### Algorithm

As in ADMIXTURE (Alexander et al. 2009; Alexander and Lange 2011), each
run consists of:

1.  a random initialisation followed by a few EM iterations
    (`n_em_iter`; the EM algorithm of FRAPPE, Tang et al. 2005), or,
    with `init = "svd"`, a start from the SVD of the genotypes (see
    below);

2.  block-wise Newton updates (block relaxation): all of Q given P (one
    individual at a time, on the simplex), then all of P given the new Q
    (one locus at a time, in the unit interval). Each Newton step is
    solved exactly as a small quadratic program, by an active-set method
    based on the sweep operator (Lange 2010);

3.  quasi-Newton acceleration (Zhou, Alexander and Lange 2011): each
    iteration applies two rounds of block updates, and the changes that
    they made in the last `n_qn` iterations (secant pairs) give a
    quasi-Newton step, which is accepted only when it improves the
    log-likelihood of the previous iteration (otherwise the iteration
    ends with the plain block updates).

Iterations stop when the relative change in log-likelihood falls below
`rtol`. All estimates are kept within \[1e-5, 1 - 1e-5\], as in
ADMIXTURE. The log of each run (element `log` of the result) reports the
log-likelihood at each iteration, together with that of the quasi-Newton
point ("QN") and that of the plain block updates ("basic"). The latter
is only computed when the quasi-Newton point is rejected (otherwise it
is reported as "not computed"), which saves a pass over the genotypes in
most iterations.

A few numerical safeguards handle degenerate cases, without changing the
results when the computations are well behaved: for example, a Newton
step that is not finite (because a genotype is impossible under the
current estimates) is discarded, and the individual or locus keeps its
current value for that iteration.

### Genotype likelihoods

With genotype likelihoods, the genotypes are not observed: each
individual has the likelihoods of its data (its reads) given each of the
three genotypes at a locus, with 0, 1 and 2 copies of `allele2`.
tidygenclust does not store them: a package that does can provide a
method of
[`gt_radmix_data()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_radmix_data.md)
for its class. The model is that of NGSadmix (Skotte et al. 2013): the
genotypes follow the ADMIXTURE model, with `P` the frequency of
`allele2`, and the log-likelihood is \$\$\sum\_{ij} \log \sum\_{g=0}^{2}
P(D\_{ij} \| g) P(g \| x\_{ij}),\$\$ where \\x\_{ij} = \sum_k q\_{ik}
p\_{kj}\\, \\P(g \| x)\\ is binomial (Hardy-Weinberg proportions), and
the likelihoods \\P(D\_{ij} \| g)\\ of each individual are normalised to
sum to 1, as in the beagle files read by NGSadmix. Individuals without
data at a locus do not contribute to it, whereas NGSadmix counts
\\\log(1/3)\\ for each of them (they have three equal likelihoods):
otherwise, the log-likelihoods of the two programs are the same.

The algorithm is the same as for genotypes: every step uses the expected
numbers of copies of each allele given the likelihoods and the current
estimates in place of the observed ones. EM steps are then those of
NGSadmix, and the gradients of the block updates are exact; as the
log-likelihood of genotype likelihoods is not concave in each block, the
Newton steps use the curvature of the expected log-likelihood of the
genotypes (that of the EM algorithm), which keeps every quadratic
program convex. The quasi-Newton steps are accepted on the
log-likelihood above. Cross-validation (`crossval`) is available (see
"Cross-validation"), and so are starts from the SVD (`init = "svd"`):
the SVD is then that of the posterior dosages of the genotypes (their
expected numbers of copies of `allele2`) under Hardy-Weinberg
proportions with the allele frequencies estimated from the likelihoods
(emMAF, Kim et al. 2011), as in the first step of PCAngsd (Meisner and
Albrechtsen 2018), computed by a randomised SVD (Halko et al. 2011)
whose products are passes over the likelihoods; individuals without data
at a locus get the mean dosage. `use_imputed` is not available for
genotype likelihoods. Unlike NGSadmix, which by default drops loci with
a minor allele frequency below 0.05 (its option `-minMaf`),
`gt_radmix()` uses all the loci of `x`: select them beforehand if
needed.

### Supervised analyses

As with ADMIXTURE's `--supervised` option, the ancestry of some
individuals can be given rather than estimated: `supervised` gives the
population of each of these reference individuals (and `NA` for the
individuals to estimate), and each population is one of the `k`
ancestral populations, so `k` must equal the number of populations in
`supervised`. The ancestry proportions of each reference individual are
fixed at 1 for its population and 0 for the others (moved inside the
feasible set like every estimate, i.e. 1 - (k - 1) 1e-5 and 1e-5, as
reported by ADMIXTURE); the allele frequencies of the ancestral
populations (P) are estimated from all the individuals, and the ancestry
of the other individuals is estimated given them. Each column of Q and P
is one reference population: in the order of the levels of `supervised`
if it is a factor, or of their first appearance otherwise (the order
used by ADMIXTURE for the populations of its `.pop` file). The
populations are listed in the element `ref_pops` of the result (the
columns of Q keep the names `.Q1`, `.Q2`, ..., which tidypopgen's
methods expect).

The fit is the same as in an unsupervised analysis (random start, EM
priming, block updates and quasi-Newton steps), with the reference
individuals fixed from the start; their ancestry proportions are never
updated, and the passes over the genotypes skip the work that only
serves them, so a supervised fit with many reference individuals is
faster than an unsupervised one. Cross-validation (`crossval = TRUE`) is
available, as in ADMIXTURE: the reference individuals stay fixed in the
fit of each fold, and their genotypes belong to the folds and are
predicted like the others (ADMIXTURE does not document this detail).
`init_q` (whose rows for the reference individuals are ignored) and
`init_p` can give starting values; `fix_q`, `fix_p` and `init = "svd"`
cannot be combined with `supervised`. Supervised analyses are also
available for genotype likelihoods.

When the reference individuals are many, starting P from the allele
frequencies of the reference populations, as fastmixture does in its
supervised mode, can save iterations: pass them as `init_p` (the
frequency of the alternate allele, e.g. from
[`tidypopgen::loci_alt_freq()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_alt_freq.html)
on the individuals of each population, in the order of the columns),
with `n_em_iter = 0` (the individuals to estimate then start from random
values). On the HGDP data (K = 5, 660 of 1043 individuals as
references), this took 7 iterations instead of 11 and 5 EM steps, to the
same log-likelihood; with fewer references (`EUR_subset`, 102 of 379),
it made no difference.

### Starting values from the SVD

With `init = "svd"`, each run starts from the principal components of
the genotypes rather than from random values, as in ADAMIXTURE
(Saurina-i-Ricos et al. 2026) and fastmixture (Santander et al. 2024):

1.  a partial SVD of the genotypes, centred and scaled as by
    [`tidypopgen::gt_pca_partialSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_partialSVD.html),
    with `max(k) + 5` components, is computed once and shared by all the
    runs (monomorphic loci, which the default scaling of tidypopgen
    rejects, are left unscaled);

2.  for each run, alternating least squares (ALS) find Q and P whose
    product approximates the allele frequencies predicted by the first
    `k + 5` components, starting from random allele frequencies. As in
    ADAMIXTURE, the least squares are bounded, the loci are weighted by
    the inverse of their binomial variance when solving for Q, and Q is
    extrapolated (momentum) when solving for P. The ALS only uses the
    SVD, not the genotypes, so it is fast;

3.  the best of 3 such solutions, by log-likelihood, starts the fit,
    without EM priming (`n_em_iter` is not used).

The SVD cannot use missing genotypes, so they must be imputed first
(e.g. with
[`tidypopgen::gt_impute_simple()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_impute_simple.html));
the imputed genotypes are only used for the SVD, and the fit ignores
them (unless `use_imputed = TRUE`). The SVD also needs more individuals
and loci than `max(k) + 5`. With genotype likelihoods, the SVD is that
of the posterior dosages (see "Genotype likelihoods"), and needs no
imputation; its random test matrix is drawn after setting the first
value of `seed`, and its number of power iterations is `svd_power` (3 by
default, enough for the leading components, those of the structure; each
further iteration makes the trailing components more accurate at the
cost of two more passes over the likelihoods).

Fits from the SVD need fewer iterations (often less than half as many as
from random starts in our tests), and usually reach a good fit (though,
as from any start, not always the best one). They are therefore useful
for a single good fit of a large dataset, where many random runs would
be costly, or as a first exploration before a more extensive search with
random starts. However, the starts of different runs vary little, and
favour some solutions over others: to explore the modes of the
likelihood with several runs (e.g. with
[`gt_clumppling()`](https://evolecolgroup.github.io/tidygenclust/dev/reference/gt_clumppling.md)),
use random starts (the default), and do not mix runs started in the two
ways, as the runs started from the SVD would over-represent the
solutions that they favour.

### Cross-validation

With `crossval = TRUE`, the cross-validation error of ADMIXTURE
(Alexander and Lange 2011) is computed for each run, and can be used to
choose `k` (lower is better; see `autoplot(type = "cv")`). The observed
genotypes are split into `n_folds` folds of roughly equal size; for each
fold, the model is refitted with the genotypes of that fold treated as
missing, and each of them is predicted by its expected value \\\mu = n
\sum_k q\_{ik} p\_{kj}\\, where \\n\\ is the number of allele copies of
the individual. The error is the binomial deviance \$\$d(g, \mu) = g
\log(g / \mu) + (n - g) \log((n - g) / (n - \mu))\$\$ averaged over all
the masked genotypes of all folds. Where ADMIXTURE does not document the
details, we made the following choices:

- each fold is fitted starting from the estimates for the full data,
  with at most `max_iter` iterations and a convergence threshold of 1e-7
  (or `rtol`, if larger) by default (see below for `cv_max_iter` and
  `cv_rtol`);

- genotypes are assigned to folds by a hash of their position and of a
  seed drawn from the random number generator (reproducible with
  `seed`), so the folds differ from those of ADMIXTURE and the errors
  agree with ADMIXTURE's closely but not exactly;

- pseudohaploid individuals use \\n = 1\\.

With genotype likelihoods, the folds are drawn over the individuals and
loci with data, in the same way, and the error of a masked individual
and locus is the deviance of its likelihoods, \$\$d = \log \max_x f(x) -
\log f(x\_{ij}),\$\$ where \\f(x) = \sum\_{g} P(D\_{ij} \| g) P(g \|
x)\\ is the likelihood of its data given a frequency \\x\\ of `allele2`
(see "Genotype likelihoods"), maximised over \\x\\ in \\\[0, 1\]\\, and
\\x\_{ij}\\ is the prediction of the refitted model. For the likelihoods
of a hard call, this is the binomial deviance above (with \\n = 2\\);
with uncertain genotypes, the errors are smaller, as the data say less
about each genotype, so they can be compared across values of `k` but
not with those of genotypes.

Cross-validation does not change the estimates of Q and P, but requires
`n_folds` additional fits for each run (each usually converging in a few
iterations).

The refits of the folds stop earlier than the fit to the full data
(`cv_rtol = 1e-7`, against `rtol = 1e-8`), and can be shortened further
with `cv_max_iter` (at most that many iterations per fold; ADAMIXTURE
refits each fold with 5) or a larger `cv_rtol`. This saves time, but
biases the errors downwards: each refit starts from the estimates for
the full data, which have seen the masked genotypes, and a refit stopped
early still predicts them too well, the more so at larger `k`. In our
tests, the default halved the cost of cross-validation on EUR_subset
(379 individuals, 54,051 loci, `k` = 1-6, 5 seeds), where refits run to
`rtol` took 24 iterations on average, and lowered the errors by at most
1.5e-5 up to `k` = 4, and by 1.2e-4 and 2e-4 at `k` = 5 and 6, while the
errors of successive values of `k` differed by 5e-4 to 3e-3. On the HGDP
(1,043 individuals, 130,943 loci, `k` = 2-10, 3 seeds), where refits
converge in about 7 iterations, it saved 11% of the passes and changed
the errors by less than 1e-6. Shorter refits save more on EUR_subset
(`cv_rtol = 1e-6`: 74%; `cv_max_iter = 5`: 80%) with larger biases (up
to 6.7e-4 and 1e-3 at `k` = 6). None of these choices changed the value
of `k` with the lowest error, for any seed. To refit each fold to
convergence, set `cv_rtol = rtol`.

### Reproducibility

The random starting values (and the folds of cross-validation) are drawn
from R's random number generator, so runs are reproducible with `seed`
(or [`set.seed()`](https://rdrr.io/r/base/Random.html)). The work is
divided between threads in a way that depends only on the dimensions of
the data, so the results are identical for any `n_cores`.

On x86-64 Linux, when the package is compiled with GCC 12 or later, the
main computations are also compiled for processors with the AVX2 and FMA
instructions (most processors made since 2015), and the faster version
is chosen when the package is loaded. Results are reproducible on a
given machine, but can differ in the last digits between machines with
and without these instructions. Compiling the package with
`-DRADMIX_NO_MULTIVERSION` (e.g. added to `CXXFLAGS` in `~/.R/Makevars`)
builds only the portable version, whose results are the same on every
x86-64 machine.

### Comparison with ADMIXTURE

The model, the log-likelihood and the optimisation strategy are those of
ADMIXTURE (Alexander et al. 2009; Alexander and Lange 2011), so both
programs reach the same maximum likelihood estimates, up to local optima
and the convergence criterion (ADMIXTURE stops when the log-likelihood
increases by less than 1e-4 between iterations). The test suite checks
this on the `EUR_subset` data of SnpArrays.jl. Note that:

- ADMIXTURE reports the frequency of the second allele in the PLINK bim
  file, which for a `gen_tibble` is the reference allele, so its P
  matrices are 1 - P of those returned here;

- cross-validation follows ADMIXTURE (see above), and the test suite
  checks that the errors agree with those of ADMIXTURE on `EUR_subset`;

- ADMIXTURE treats the autosomal genotypes of every individual as
  diploid, whereas `gt_radmix()` counts a single allele per locus for
  pseudohaploid individuals (see above);

- supervised analyses (`supervised`) follow ADMIXTURE's `--supervised`
  option: the same fixed ancestry of the reference individuals, the same
  order of the populations, and the same requirement that `k` equals
  their number; `supervised` gives the populations as a vector (with
  `NA` for the individuals to estimate) instead of a `.pop` file (with
  `-`). The test suite checks the estimates against those of ADMIXTURE
  on `EUR_subset`;

- bootstrap standard errors, available in ADMIXTURE, are not
  implemented.

## References

Alexander DH, Novembre J, Lange K (2009) Fast model-based estimation of
ancestry in unrelated individuals. Genome Research, 19(9), 1655-1664.

Alexander DH, Lange K (2011) Enhancements to the ADMIXTURE algorithm for
individual ancestry estimation. BMC Bioinformatics, 12, 246.

Zhou H, Alexander D, Lange K (2011) A quasi-Newton acceleration for
high-dimensional optimization algorithms. Statistics and Computing,
21(2), 261-273.

Tang H, Peng J, Wang P, Risch NJ (2005) Estimation of individual
admixture: analytical and study design considerations. Genetic
Epidemiology, 28(4), 289-301.

Lange K (2010) Numerical Analysis for Statisticians, 2nd edition.
Springer, New York.

Ko S, Chu BB, Peterson D, Okenwa C, Papp JC, Alexander DH, Sobel EM,
Zhou H, Lange KL (2023) Unsupervised discovery of ancestry-informative
markers and genetic admixture proportions in biobank-scale datasets. The
American Journal of Human Genetics, 110(2), 314-325.

Saurina-i-Ricos J, Mas Montserrat D, Ioannidis AG (2026) ADAMIXTURE:
adaptive first-order optimization for biobank-scale genetic clustering.
Bioinformatics, 42(Supplement_1). doi:
https://doi.org/10.1093/bioinformatics/btag236

Santander CG, Refoyo Martinez A, Meisner J (2024) Faster model-based
estimation of ancestry proportions. Peer Community Journal, 4, e115.
doi: https://doi.org/10.24072/pcjournal.503

Skotte L, Korneliussen TS, Albrechtsen A (2013) Estimating individual
admixture proportions from next generation sequencing data. Genetics,
195(3), 693-702.

Kim SY, Lohmueller KE, Albrechtsen A, Li Y, Korneliussen T, Tian G,
Grarup N, Jiang T, Andersen G, Witte D, Jorgensen T, Hansen T, Pedersen
O, Wang J, Nielsen R (2011) Estimation of allele frequency and
association mapping using next-generation sequencing data. BMC
Bioinformatics, 12, 231.

Meisner J, Albrechtsen A (2018) Inferring population structure and
admixture proportions in low-depth NGS data. Genetics, 210(2), 719-731.

Halko N, Martinsson PG, Tropp JA (2011) Finding structure with
randomness: probabilistic algorithms for constructing approximate matrix
decompositions. SIAM Review, 53(2), 217-288.

## Examples

``` r
bed_file <-
  system.file("extdata", "pop_a.bed", package = "tidypopgen")
pop_a <- tidypopgen::gen_tibble(bed_file,
  backingfile = tempfile("pop_a"),
  quiet = TRUE
)
adm <- gt_radmix(pop_a, k = 2, seed = 1)
adm$Q[[1]]
#>            .Q1       .Q2
#> [1,] 0.5051543 0.4948457
#> [2,] 0.0000100 0.9999900
#> [3,] 0.9999900 0.0000100
#> [4,] 0.0000100 0.9999900
#> [5,] 0.9999900 0.0000100
#> attr(,"class")
#> [1] "q_matrix" "matrix"   "array"   

# choose k by cross-validation
adm_cv <- gt_radmix(pop_a, k = 1:3, seed = 1, crossval = TRUE)
adm_cv$cv
#> [1] 0.8224648 1.1840311 0.7096221

# supervised analysis: the first two individuals are references of
# population "A" and the next two of population "B"; the ancestry of the
# others is estimated
pops <- rep(NA, nrow(pop_a))
pops[1:4] <- c("A", "A", "B", "B")
adm_sup <- gt_radmix(pop_a, k = 2, seed = 1, supervised = pops)
adm_sup$ref_pops
#> [1] "A" "B"
adm_sup$Q[[1]]
#>          .Q1     .Q2
#> [1,] 0.99999 0.00001
#> [2,] 0.99999 0.00001
#> [3,] 0.00001 0.99999
#> [4,] 0.00001 0.99999
#> [5,] 0.00001 0.99999
#> attr(,"class")
#> [1] "q_matrix" "matrix"   "array"   
```

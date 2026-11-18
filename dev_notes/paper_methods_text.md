# Text for the methods of the paper

Text on how the native implementations of tidygenclust work and how they
compare with the programs they replace, written on 18 Nov 2026 (manica01)
from the code of `dev` (54ce6cd) and the sources of the other programs, with
the sources of each statement. The code derived from other software, and
its licences, are listed in `inst/COPYRIGHTS` (and, for `gt_radmix()`, in
"Origins of the code" of `admixture_implementation.md`).

# `gt_radmix()` (radmix)

Two paragraphs on how radmix works and how it compares with the programs
of the benchmark (ADMIXTURE, fastmixture, ADAMIXTURE and OpenADMIXTURE.jl),
from the sources of OpenADMIXTURE.jl (commit 401ac5c), ADAMIXTURE 1.8.9 and
fastmixture 1.3.0. The design is in `admixture_implementation.md`.

## The paragraphs

radmix (the function `gt_radmix()` of tidygenclust) estimates individual
ancestry proportions (Q) and ancestral allele frequencies (P) by maximum
likelihood under the model of ADMIXTURE (Alexander et al. 2009): the
genotype of individual *i* at locus *j* is binomial, with *n_i* draws (two
for diploids, one for pseudohaploids) and success probability
Σ_k q_ik p_kj, where each row of Q lies on the simplex and each p_kj in
[0, 1]. It follows the optimisation strategy of ADMIXTURE: a random start
is refined by a few EM steps (Tang et al. 2005), and the parameters are
then updated by block relaxation, in which Q given P and then P given Q are
updated, each individual and each locus by its own constrained Newton step;
every two rounds of block updates are accelerated by the quasi-Newton
extrapolation of Zhou, Alexander and Lange (2011), built from the last three
secant pairs. Because the log-likelihood is concave in each row of Q and of
P, each Newton step is a small convex quadratic program, which is solved
exactly by Lange's active-set method based on the sweep operator (Lange
2010). radmix began as a C++ port of OpenADMIXTURE.jl (Ko et al. 2023), a
Julia implementation of this algorithm, but has since been extensively
rewritten; the parts that remain translations of OpenADMIXTURE.jl (the
solver of the quadratic programs, the projections onto the feasible set,
the quasi-Newton step and the structure of the fitting loop) are labelled
as such in the source code, which reproduces their MIT licence. Among the
changes, the derivatives with respect to Q are those of the log-likelihood
itself (OpenADMIXTURE.jl differentiates a homogeneous form of it, which
doubles the work of each Hessian), so that the two blocks share the same
per-genotype coefficients; the genotypes are read directly from the
memory-mapped backing file of a `gen_tibble`, in tiles of 64 individuals by
64 loci, with branch-free vectorised arithmetic; the quadratic program of
each individual or locus is solved as soon as its Hessian has been
accumulated, so that the Hessians are never stored; the log-likelihood of a
candidate point is computed in the same pass over the genotypes as the next
update of Q, so that an accepted iteration takes four passes over the data
instead of six; and all sums are added in a fixed order, so that the
estimates are identical for any number of threads. radmix also adds
pseudohaploid data, the cross-validation of ADMIXTURE (the binomial
deviance of held-out genotypes; Alexander and Lange 2011), supervised
analyses, and genotype likelihoods, for which it fits the model of NGSadmix
(Skotte et al. 2013) by replacing each genotype with its expected allele
counts; the selection of ancestry-informative markers by sparse K-means and
the GPU support of OpenADMIXTURE.jl were not ported.

ADMIXTURE (distributed as closed-source binaries) and OpenADMIXTURE.jl thus
maximise the same likelihood with essentially the same algorithm as radmix,
and differ from it in their implementation and stopping rules: ADMIXTURE
stops when the log-likelihood increases by less than 10^-4 between
iterations, OpenADMIXTURE.jl when its relative change falls below 10^-7,
and radmix below 10^-8. The test suite of radmix checks its estimates,
log-likelihoods and cross-validation errors against stored results of
ADMIXTURE 1.3.0, and its updates and log-likelihoods against those of
OpenADMIXTURE.jl. fastmixture (Santander et al. 2024) and ADAMIXTURE
(Saurina-i-Ricos et al. 2026) differ in their starting values and, for
fastmixture, in the optimiser. Both start by default from a truncated
singular value decomposition (SVD) of the genotypes, from which alternating
least squares (ALS) extract non-negative factors that serve as initial Q and
P. fastmixture then maximises the likelihood by EM, with a quasi-Newton
extrapolation after every two EM steps, using mini-batches of loci in the
early iterations. ADAMIXTURE was introduced with an EM algorithm combined
with the Adam adaptive optimiser, and also offers ADMIXTURE's block
relaxation with quasi-Newton acceleration, which is its default solver in
version 1.8.9. radmix offers an SVD start as an option (`init = "svd"`): a
partial SVD of the centred and scaled genotypes with K + 5 components,
followed by ALS adapted from those of ADAMIXTURE (bounded least squares,
loci weighted by the inverse of their binomial variance when solving for
Q, and momentum extrapolation of Q; BSD 3-clause licence, labelled in the
source code), keeping the best of three ALS solutions by log-likelihood. An
SVD start typically halves the number of iterations, but repeated runs from
it are nearly deterministic and so explore few modes of the likelihood.
radmix therefore keeps random starts as its default, as the alignment of
replicate runs across modes (e.g. with `gt_clumppling()`) depends on that
exploration.

## A sentence for the acknowledgements or the code availability

radmix includes code translated from OpenADMIXTURE.jl (MIT licence;
Copyright (c) 2022 Seyoon Ko, Kenneth Lange and contributors) and adapted
from ADAMIXTURE (BSD 3-clause licence; Copyright (c) 2024 Ioannidis Lab);
the files concerned are listed, with the licences, in the file `COPYRIGHTS`
of tidygenclust.

## Sources of the statements

- The algorithm and its defaults (5 EM steps, 3 secant pairs, `rtol =
  1e-8`, K + 5 components and 3 ALS solutions): `admixture_implementation.md`
  and `src/radmix_fit.cpp`, `R/radmix_svd.R`.
- Four passes per accepted iteration: `radmix_fit.cpp` (three block updates
  and the fused pass); OpenADMIXTURE.jl's `admixture_qn!()` makes four block
  updates and two log-likelihood passes per iteration.
- The homogeneous form: `update_q_loop!()` in OpenADMIXTURE.jl's
  `src/loops.jl`, whose Hessian adds `p p'` and `(1 - p)(1 - p)'` terms for
  each genotype (radmix adds one).
- Stopping rules: "How the benchmark was run" in
  `vignettes/articles/admixture_benchmark.Rmd`.
- fastmixture: `functions.py` (`batQuasi()`: two EM steps, then `jumpP()`,
  a step length from one secant pair) and `utils.py` (`factorALS()`).
- ADAMIXTURE: `api.py` ("default: brqn + als"), `model/br_qn.py` (SQP block
  updates with ZAL quasi-Newton), `model/als.py`.
- SVD start halving the iterations, and why it is not the default:
  `admixture_implementation.md`, section 3.

## Open questions

- Speed figures are left out: add them from the benchmark article
  (`admixture_benchmark.Rmd`) once the paper's runs are fixed.
- The references in the text are those of `?gt_radmix`; ADMIXTURE's
  version and the stopping rule should be checked against the version used
  in the paper's runs.

# `gt_clumppling()`

One short paragraph (Andrea asked for a short one), from the code of `dev`
and the sources of clumppling 2.0.10, python-louvain 0.16,
markov_clustering 0.0.6.dev0 and cdlib 0.4.1. The design is in
`clumppling_implementation.md`.

## The paragraph

`gt_clumppling()` implements the method of clumppling (Liu, Kopelman and
Rosenberg 2024) in R and C++, as a port, with modifications, of its python
package (version 2.0.10). Every pair of replicate runs with the same K is
aligned by the matching of clusters that minimises the squared Euclidean
distance between their membership matrices; the replicates are grouped
into modes by community detection (Louvain by default) on the network of
alignment costs; and the modes are aligned across K, merging two clusters
when two values of K differ by one, so that clusters are labelled
consistently across K. Instead of solving an integer linear program with
GLPK for each pair, `gt_clumppling()` solves the equivalent linear
assignment problem exactly with the Hungarian algorithm (with extra "free"
columns when the two runs have different K), and community detection uses
C++ and R ports of the implementations that clumppling calls through cdlib
(python-louvain and markov_clustering) or the igraph package; the ported
code is labelled in the source, with its licences. On the Cape Verde
example of clumppling (200 ADMIXTURE runs, K = 2-5), `gt_clumppling()`
finds the same modes, representatives and alignments as the python package
(to within 2 x 10^-16) in 0.08 s instead of 17 s; the results differ only
where the community structure is ambiguous, and there the python results
themselves change with the random seed.

## A sentence for the acknowledgements or the code availability

`gt_clumppling()` is ported from clumppling (MIT licence; Xiran Liu), and
includes code ported from python-louvain (BSD 3-clause licence; Thomas
Aynaud) and markov_clustering (MIT licence; Guy Allard), and adapted from
the implementation of the Hungarian algorithm by Andrey Lopatin published
on cp-algorithms (CC BY-SA 4.0); the files concerned are listed, with the
licences, in the file `COPYRIGHTS` of tidygenclust.

## Sources of the statements

- The algorithm and the reduction of the ILP to an assignment problem:
  `src/clumppling_assignment.cpp` (top of the file) and `?gt_clumppling`;
  defaults `cd_method = "louvain"`, `merge = TRUE`, `use_rep = TRUE`.
- Same results and times: `vignettes/articles/data/clumppling_benchmark/`
  (`agreement.csv`: capeverde_50, 200 runs, same modes and
  representatives, largest difference of the aligned modes 2.2e-16;
  `timings.csv`: medians 0.077 s in R and 16.6 s in python, both with the
  Q matrices in memory; laptop i7-10510U, 31 Oct 2026, tidygenclust
  0.1.5.9003). The chicken example (K = 17, 19, 21) differs at K = 17 and
  19, where python's own modes change with the seed
  (`clumppling_implementation.md`). Also `validate_vs_python.csv` (60
  random scenarios, same modes in all).

## Open questions

- If the paper quotes the times, rerun the benchmark with the version of
  the paper (the numbers above are from 31 Oct 2026, before the
  simplification of #34).

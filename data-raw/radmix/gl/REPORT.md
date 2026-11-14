# `gt_radmix()` on genotype likelihoods: benchmarks

`gt_radmix()` fits the model of NGSadmix (Skotte et al. 2013) to genotype
likelihoods, with the algorithm of ADMIXTURE. tidygenclust does not store
genotype likelihoods: a package that does gives them to `gt_radmix()`
through a method of `gt_radmix_data()`, as a "source" (a C function that
fills the likelihoods of a block of individuals at one locus).
`dev_notes/admixture_implementation.md` (section "Genotype likelihoods")
describes how; `tests/testthat/test_radmix_vs_ngsadmix.R` checks the
estimates and the log-likelihood against NGSadmix. This directory holds the
benchmarks, on simulated data, with the example source of tidygenclust
(`src/radmix_gl_source.cpp`, which keeps three doubles per individual and
locus in memory; a source reading a compact encoding from disk has other
costs).

All the measurements use builds of tidygenclust installed with
`R CMD INSTALL` (`-O2`, with the AVX2 versions of the accumulation
functions), never `devtools::load_all()`, on a laptop with an Intel Core
i7-10510U (4 cores, AVX2) and 15 GB of RAM: `dev` (9b5f046) and `gl` (this
branch).

## Scripts

| file | what it does |
|---|---|
| `benchmark_gl.R` | the time of each pass over the data, and of a full fit, with genotypes and with likelihoods (and of NGSadmix on the same likelihoods, with `NGSADMIX` set), on simulated data; appends to `results/benchmark_gl.csv` |
| `perf_gl.R` | one pass (or a fit) of `benchmark_gl.R`, on genotypes or likelihoods, counted by `perf stat` (user-space instructions and cycles of the measured part only), to compare builds when the machine is not quiet; `results/perf_gl.csv` |

The data: 200 individuals and 50,000 loci from K = 3 ancestral populations
(`radmix_simulate_admixed()`, 1% missing genotypes), and genotype
likelihoods from reads at a mean depth of 4 with an error rate of 0.01
(`radmix_simulate_gl()`; 2.8% of the individuals have no reads at a
locus). Fits start from the same random start (seed 1), with the default
`rtol = 1e-8`.

## The passes over genotypes, and the cost of likelihoods

The machine was in use (wall-clock times of the same build varied by up to
80% between runs), so the builds are compared by the instructions and
cycles that `perf_gl.R` counts: two rounds, the builds interleaved, one
thread, 20 repetitions of each pass after a warm-up (load 1.2 to 1.8).
Instruction counts do not depend on the load, and cycles are only counted
while the process runs. Millions of cycles per pass (instructions in
parentheses), and billions for a fit, the range of the two rounds:

| pass (one thread) | `dev`, genotypes | `gl`, genotypes | `gl`, likelihoods |
|---|---|---|---|
| log-likelihood | 258-263 (882) | 259 (880) | 261 (802) |
| EM step | 263 (686) | 266 (682) | 453-456 (926) |
| update of Q | 148-149 (403) | 151 (401) | 311-313 (696) |
| update of P | 204-205 (553) | 205-206 (551) | 282-283 (693) |
| fit | 10.24-10.25 (29.7), 9 iterations | 10.35-10.37 (29.6), 9 iterations | 26.0-26.1 (64.1), 16 iterations |

- The passes over genotypes execute the same instructions as on `dev`
  (0.2 to 0.6% fewer), and take at most 1.5% more cycles (1.1% for a fit).
  The kind of data is a template parameter of the functions that process
  slices, chosen once per pass, so the passes over genotypes compile as
  before.
- A pass over likelihoods takes 1.0 (log-likelihood) to 2.1 (update of Q)
  times the cycles of a pass over genotypes: the source reads 24 bytes per
  individual and locus instead of 1, and the expected allele counts cost a
  division and a few multiplications more than the decoding of a genotype.
  A fit takes 2.5 times as many cycles, as it also needs more iterations:
  the Newton steps use the curvature of the expected log-likelihood of the
  genotypes, which is larger than that of the log-likelihood when genotypes
  are uncertain, so the steps are shorter.

## NGSadmix on the same likelihoods

`benchmark_gl.R` with 4 threads (`BENCH_CORES=4`), with NGSadmix (ANGSD
0.940; 4 threads, `-minMaf 0` so that both programs fit the same loci, and
its default convergence criteria) on the same likelihoods written to a
beagle file (normalised, 8 significant digits):

| program | seconds | iterations | log-likelihood (NGSadmix's convention) |
|---|---|---|---|
| `gt_radmix()` | 4.6 | 16 | -8148016.2 |
| NGSadmix | 182 | 153 | -8148007.4 |

- `gt_radmix()` is about 40 times faster (one run each, on a machine in
  use: the ratio, not the times, is what matters).
- Its maximum is 8.9 log-likelihood units lower (1e-6 relative):
  `gt_radmix()` keeps its estimates within [1e-5, 1 - 1e-5], as ADMIXTURE
  does, and NGSadmix within [1e-9, 1 - 1e-9]. Many estimates are at the
  bound here (the simulated allele frequencies are drawn from
  Beta(0.5, 0.5), many close to 0 or 1, and most individuals have most of
  their ancestry from one population). In the test against NGSadmix, our
  maximum is at least the log-likelihood of NGSadmix's estimates moved
  inside our bounds.
- NGSadmix counts `log(1/3)` in its log-likelihood for each individual
  without data at a locus, which `gt_radmix()` leaves out; the table adds it
  to `gt_radmix()`'s log-likelihood.

# Implementation of `gt_radmix()`

`gt_radmix()` estimates ancestry proportions (Q) and ancestral allele
frequencies (P) with the method of ADMIXTURE (Alexander et al. 2009). It
uses a modified version of the algorithm of OpenADMIXTURE (Ko et al. 2023),
and much of its code is ported to C++ (with OpenMP) from the Julia
implementation, OpenADMIXTURE.jl. It reads the genotypes directly from the
backing file of a `gen_tibble`. It is not called OpenADMIXTURE because it
differs from it: it does not include the selection of ancestry informative
markers, we implemented the cross-validation of ADMIXTURE (which is not part
of OpenADMIXTURE), and it changes the algorithm in the ways described
below.

This document describes how the algorithm works, step by step. It then lists
the changes from OpenADMIXTURE.jl that make it faster (and use less memory),
and finally the changes made for other reasons. The port was made from
OpenADMIXTURE.jl at commit `401ac5c4056f0f6f2a7ea6ffc3af69c2b1e56fdc`, and
the source files of that version are cited below by name (e.g. `loops.jl`).

Starting from the same values, `gt_radmix()` reproduces the
computations of OpenADMIXTURE.jl to numerical precision. The test suite checks
this step by step, against stored results of the Julia code
(`tests/testthat/test_radmix_vs_julia.R`). It also checks that the
estimates and cross-validation errors agree with those of ADMIXTURE 1.3.0
(`tests/testthat/test_radmix_vs_admixture.R`).

## Where the code is

| File | Role |
|---|---|
| `R/gt_radmix.R` | User interface: checks the arguments, prepares the data once, loops over K and runs, assembles the `gt_admix` object |
| `R/radmix_data.R` | Builds the description of the genotypes passed to C++ (backing file, selected rows and columns, decoding table, ploidy) and the cross-validation masks |
| `R/radmix_fit.R` | Thin wrapper around the C++ fitting loop |
| `R/radmix_cv.R` | Cross-validation, as in ADMIXTURE |
| `src/radmix_fit.cpp` | The fitting loop: random start, EM priming, block updates, quasi-Newton acceleration, convergence |
| `src/radmix_kernels.cpp` | Passes over the genotypes: log-likelihood, EM step, gradients and Hessians, Newton updates of Q and P, deviance |
| `src/radmix_qp_sweep.cpp` | Solver of the small quadratic program of each individual and locus (sweep operator, active set) |
| `src/radmix_quasi_newton.cpp` | Quasi-Newton step from the secant pairs |
| `src/radmix_projections.cpp` | Projections onto the simplex (Q) and the box (P) |
| `src/radmix_geno_access.h` | Read-only view of the genotypes in the backing file, with cross-validation masks |

The C++ files are commented in detail, and each change from the Julia code
is documented where it is implemented.

## The model

For individual `i` and locus `j`, the genotype `g_ij` is the number of copies
of the alternate allele. It is out of `n_i = 2` for diploids, or `n_i = 1`
for pseudohaploids. It is binomial:

```
g_ij ~ Binomial(n_i, qp_ij),    qp_ij = sum_k q_ki p_kj
```

where `q_ki` is the proportion of the genome of individual `i` that comes from
ancestral population `k`, and `p_kj` the frequency of the alternate allele at
locus `j` in population `k`. The log-likelihood, over the observed genotypes
(missing genotypes do not contribute), is

```
L = sum_ij g_ij log(qp_ij) + (n_i - g_ij) log(1 - qp_ij)
```

It is maximised subject to each column of Q lying on the simplex
(`q_ki >= 0`, `sum_k q_ki = 1`) and each `p_kj` lying in `[0, 1]`. All
estimates are kept within `[1e-5, 1 - 1e-5]`, as in ADMIXTURE and
OpenADMIXTURE.jl.

Internally, Q is stored as a K x I matrix and P as a K x J matrix,
column-major. The K values of one individual (or one locus), which are
updated together, are therefore contiguous.

## The algorithm, step by step

### 1. Preparing the data (once per call)

`radmix_prepare()` gathers what the C++ code needs, without copying
the genotypes:

- the file-backed matrix (FBM) of the `gen_tibble`, one byte per genotype,
  memory-mapped;
- the FBM row of each selected individual and the FBM column of each selected
  locus. Filtering, reordering and `select_loci()` only change these indices,
  never the backing file;
- a 256-entry decoding table of the bytes (0, 1, 2 are dosages, everything
  else is missing, and the imputed codes 4-6 are used only with
  `use_imputed = TRUE`). This makes the result independent of the `code256`
  currently attached to the R object;
- the number of allele copies `n_i` of each individual (2, or 1 for
  pseudohaploids), and the factor that turns a stored dosage into a count out
  of `n_i` (pseudohaploids are stored as 0/2).

This is shared by all the runs and values of K. In C++, a `GenoView`
(`radmix_geno_access.h`) reads genotype `(i, j)` as a memory read, a
table lookup and a multiplication.

### 2. Runs

For each value of K and each run, `gt_radmix()` sets the seed (if
given) and calls the fitting loop (`radmix_fit_cpp()`). The random
start is drawn from R's random number generator, so `set.seed()` makes fits
reproducible. The user's random number state is restored on exit.

### 3. Random start

Q and P are drawn from uniform distributions, and each column of Q is rescaled
to sum to one (as in the `AdmixData` constructor of the Julia code). Starting
values given by the user (`init_q`, `init_p`) replace the corresponding block
after the EM priming.

### 4. EM priming

`n_em_iter` (5 by default) iterations of the EM algorithm of FRAPPE move the
random start into a sensible region of the parameter space. They are skipped
when both `init_q` and `init_p` are given. Writing, for one genotype,

```
a = g / qp,   b = (n - g) / (1 - qp)
```

each EM step updates Q and P simultaneously from the current values:

```
q_new_ki = q_ki sum_j [a p_kj + b (1 - p_kj)] / (n_i J_i)
p_new_kj = A / (A + B),   A = p_kj sum_i a q_ki,   B = (1 - p_kj) sum_i b q_ki
```

where `J_i` is the number of loci observed in individual `i`. It takes two
passes over the genotypes: one over blocks of individuals, for Q, and one
over blocks of loci, for P.

### 5. Main loop

The log-likelihood of the starting point is computed. Then each iteration:

1. **Two block sweeps.** One sweep `F` updates all of Q given P, then all of
   P given the new Q. It is applied twice, giving `x0 = (Q, P)`,
   `x1 = F(x0)` and `x2 = F(x1)`. Given P, the log-likelihood separates into
   one term per individual, and given Q into one term per locus, so every
   column is updated independently, by a constrained Newton step:
   - **Gradient and negative Hessian** of each column, accumulated over the
     genotypes. With `wa = g / qp^2` and `wb = (n - g) / (1 - qp)^2`:
     ```
     Q:  grad_k = sum_j a p_kj + b (1 - p_kj)
         H_kk'  = sum_j wa p_kj p_k'j + wb (1 - p_kj)(1 - p_k'j)
     P:  grad_k = sum_i (a - b) q_ki
         H_kk'  = sum_i (wa + wb) q_ki q_k'i
     ```
     (The Q-side formulas come from the homogeneous form of the likelihood
     used by the Julia code. On the simplex they give the same Newton step
     as the plain derivatives.)
   - **Quadratic program.** The step `delta` maximises the local quadratic
     model `grad' delta - 0.5 delta' H delta` subject to
     `0 <= x + delta <= 1`, plus `sum_k delta_k = 1 - sum_k q_ki` for Q. It
     is solved exactly with K. Lange's active-set method based on the sweep
     operator (`radmix_qp_sweep.cpp`, ported from `qp.jl`). The
     method:
     - builds a tableau that bundles H, the gradient and (for Q) the sum
       constraint;
     - sweeps all the parameters to obtain the Newton step;
     - walks along that step until a parameter reaches a bound;
     - un-sweeps (fixes) that parameter and recomputes the step for the
       others;
     - sweeps back in any fixed parameter whose gradient points into the
       box, until the KKT conditions hold.

     For Q, a multiple `mu 11'` of the matrix of ones is added to H, with
     `mu` from Gerschgorin's theorem, so that the tableau can be swept even
     when H is singular along the direction of the ones. This does not change
     the solution on the simplex.
   - **Projection.** The new column `x + delta` is projected onto
     `{sum = 1, q >= 1e-5}` for Q (sort-and-threshold algorithm), or clamped
     to `[1e-5, 1 - 1e-5]` for P (`radmix_projections.cpp`).

2. **Quasi-Newton acceleration** (Zhou, Alexander and Lange 2011). The
   maximum likelihood estimate is a fixed point of `F`. Each iteration gives
   a secant pair `u = x1 - x0`, `v = x2 - x1`. The last `n_qn` pairs (3 by
   default, kept in a circular buffer) are the columns of U and V, and the
   step jumps towards the fixed point:
   ```
   x_qn = x1 - V (U'(U - V))^-1 U'(x0 - x1)
   ```
   Only an `n_qn` x `n_qn` system is solved. A system that is singular to
   machine precision is rejected, and the step is then skipped. The point is
   projected back onto the feasible set.

3. **Accept or reject.** The quasi-Newton point is kept if its
   log-likelihood is higher than that of `x0`. Otherwise the iteration ends
   at `x2`. As in the Julia code, the comparison is with `x0`, not with
   `x2`.

4. **Convergence.** The loop stops when the relative change in
   log-likelihood over an iteration falls below `rtol` (1e-7 by default), or
   after `max_iter` iterations (1000).

### 6. Cross-validation (optional)

With `crossval = TRUE`, the cross-validation error of ADMIXTURE (Alexander and
Lange 2011) is computed after the main fit, so the main estimates are the
same with or without it. The observed genotypes are split into `n_folds`
folds (5 by default). For each fold:

1. The model is refitted without the genotypes of that fold. This uses a
   "training" view of the genotypes, which hides them from every kernel
   exactly as if they were missing.
2. The genotypes of the fold are predicted by `mu = n qp`, through a
   "held-out" view that shows only them.
3. The binomial deviance of those genotypes is summed:
   ```
   d(g, mu) = g log(g / mu) + (n - g) log((n - g) / (n - mu))
   ```

The error is the deviance averaged over all the held-out genotypes of all
folds. The folds get their own seed, drawn from the random number generator
after the main fit, which is then restored. Two details that ADMIXTURE does
not document were decided here:

- each fold is fitted starting from the estimates for the full data;
- genotypes are assigned to folds by a hash of their position.

## Changes from OpenADMIXTURE.jl that speed up the algorithm

None of these changes alters the estimates of Q and P, except through
rounding. The sums are accumulated in a different order from the Julia code,
but always in the same order for any number of threads.

| OpenADMIXTURE.jl | `gt_radmix()` | Effect |
|---|---|---|
| Log-likelihood of `x2` in every iteration | Only when the quasi-Newton step is rejected | 5 passes over the genotypes per iteration instead of 6, in most iterations |
| Log-likelihood after every EM step | Only with `verbose = TRUE`, every 5 steps | Saves `n_em_iter` passes per run |
| Two logarithms per genotype | Products of 8 likelihood terms, one logarithm per product | 8 logarithms per 64 genotypes instead of up to 128 |
| Coefficients recomputed inside the loops over k and k' | Computed once per genotype, branch-free, vectorised | Fewer divisions; SIMD loops |
| Full K x K Hessians | Upper half only, with products of p hoisted out of the loop over individuals | About half the work for the Hessians |
| Hessians of all individuals (K x K x I) and all loci (K x K x J) stored, then the QPs solved in a separate pass | Accumulation and QP solved block by block | The K x K x J array (K times the size of P) is never allocated; the Hessians stay in cache |
| Q side parallelised over individuals only | Two-dimensional blocking (individuals x segments of loci) when there are few individuals | All threads busy with tens to hundreds of individuals |
| Active-set loop runs to its limit of 1000 passes when it cycles at a vertex | Cycling detected and stopped | `update_q` about 10% faster near convergence |
| Estimates copied on acceptance | Buffers swapped | No copy of the parameters per iteration |
| PLINK bed file read with SnpArrays | Backing file of the `gen_tibble` read directly, through row and column indices | No export to PLINK; any subset of individuals and loci at no cost |

### Fewer passes over the genotypes

Every pass over the genotypes reads the whole genotype matrix, so the number
of passes per iteration largely determines the run time on large datasets.

1. **Log-likelihood of `x2` only when needed** (`radmix_fit.cpp`).
   Julia computes `L(x2)` (`ll_basic`) in every iteration, before the
   quasi-Newton step. But the decision compares `L(x_qn)` with `L(x0)`, so
   `L(x2)` is only needed when the step is rejected: `x2` then becomes the
   new estimate, and its log-likelihood is used for the convergence test.
   The step is accepted in most iterations (typically 85-100%). An iteration
   therefore makes 5 passes instead of 6 (two updates of Q, two of P, one
   log-likelihood), with identical results. The log reports
   `basic = not computed` when `L(x2)` was not needed.

2. **No log-likelihood during EM priming.** Julia computes and prints the
   log-likelihood after every EM step. These values are not used by the
   algorithm (EM runs a fixed number of steps), so they are only computed
   with `verbose = TRUE`, and then every 5 steps. The log-likelihood after
   EM is the initial log-likelihood of the main loop, which is always
   computed.

3. **Batched logarithms** (`loglik_slice()` in `radmix_kernels.cpp`).
   For hard calls, the contribution of a genotype is `log(z)` with
   `z = qp^g (1 - qp)^(n - g)`, computed without logarithms or branches.
   Within a slice of 64 individuals, the z of 8 genotypes are multiplied
   together, in 8 interleaved products that vectorise, before taking a single
   logarithm. Each factor is at least 1e-30, so a product of 8 cannot
   underflow. Genotypes for which this does not apply go through the direct
   formula: dosages, and the rare z below 1e-30 (including z = 0, which
   correctly gives -Inf). The log-likelihood is computed at least once per
   iteration, so this matters.

### Faster kernels

The innermost unit of work is a *slice*: one locus and up to 64 individuals
(adjacent bytes in the column-major backing file). For each slice:

4. **Coefficients computed once, without branches, and vectorised.** `qp` is
   computed for the whole slice with unit-stride loops over a transposed copy
   of Q (made once per kernel call). The genotypes are decoded in a separate
   loop, the only one with a branch (missing or not). Missing genotypes get
   `g = n - g = 0`, so they drop out of every sum with no further tests.
   The coefficients `a`, `b`, `wa`, `wb` are then computed with arithmetic
   masks instead of tests, which also avoids `0/0` when `qp` is exactly 0 or
   1. All these loops are SIMD. Julia instead evaluates
   `gij / qp^2 * p[k] * p[k2]` (and its counterpart) inside the loops over
   `k` and `k2`, for every individual, locus, `k` and `k2`.

5. **Half Hessians and hoisted products.** The Hessians are symmetric, so
   only the `K (K + 1) / 2` entries with `k' >= k` are accumulated. The
   other half is filled in once at the end. On the Q side, the products
   `p_kj p_k'j` and `(1 - p_kj)(1 - p_k'j)` do not depend on the individual,
   so they are computed once per locus, outside the loop over the
   individuals of the slice. `1 - P` is computed once per kernel call and
   shared by all threads. On the P side, the coefficients only enter as
   `a - b` and `wa + wb`, which are formed first, halving the work.
   Accumulators are laid out so that each update is a unit-stride loop:
   - on the Q side, one vector of 64 values per `k`, or per `(k, k')` pair;
   - on the P side, a dot product over the slice, while the 64 x K corner
     of Q stays in L1 cache.

6. **Accumulation and QP solved block by block.** Julia first fills arrays
   with the Hessian of every individual (K x K x I) and of every locus
   (K x K x J), then solves the QPs in a separate pass. The array for the
   loci alone is K times the size of P. Here each thread accumulates the
   gradients and Hessians of one block of 64 individuals (or loci) into
   thread-local buffers, and solves their QPs straight away. The large
   arrays are never allocated, and each Hessian is used while it is still in
   cache.

### Parallelism

7. **Fixed tiles, no atomics, deterministic results.** Work is split into
   tiles of 64 individuals x 64 loci:
   - Q-side work (sums over loci for each individual) is parallelised over
     blocks of individuals;
   - P-side work (sums over individuals for each locus) over blocks of loci.

   Each thread owns the output columns of its blocks, so there are no races
   and no reductions across threads. The log-likelihood, the only scalar
   sum, is stored as one partial sum per block, and these are added up in
   order. Julia spawns nested tasks and adds the partial log-likelihoods
   with an atomic add, whose order can vary. The block boundaries depend
   only on the dimensions of the data, never on the number of threads, so
   the estimates are bit-identical for any `n_cores` (the HGDP benchmark in
   `vignettes/articles/admixture_benchmark.Rmd` confirms this).

8. **Two-dimensional blocking of the Q side.** With few individuals there
   are few blocks of 64. For example, 379 individuals give 6 blocks, which
   balance 4 threads poorly, and 50 individuals give a single block, so the
   Q side would run on one thread. When there are fewer than 32 blocks, the
   loci are also split into segments (of at least 256 loci), and each
   (block, segment) pair is a separate task. The partial sums of the
   segments are then added up in segment order. The number of segments
   depends only on I and J, so results remain independent of the number of
   threads. Population genetic datasets often have tens to hundreds of
   individuals and many loci, so this matters.

9. **Parallel quasi-Newton step.** The step involves vectors of length
   K (I + J). It is computed directly from the Q and P blocks, without
   stacking them into one vector or creating temporaries. Its long inner
   products are split into fixed chunks of 65,536 values, processed in
   parallel and added up in order.

### Quadratic programs

10. **Early exit when the active-set loop cycles** (`qp_solve_column()`).
    The Q solution can be a vertex of the simplex, where every parameter is
    on a bound. There, the active-set loop can enter a cycle:
    - it frees a fixed parameter whose gradient points into the box;
    - it takes a step of exactly 0, because the last free parameter,
      determined by the sum constraint, cannot move;
    - it fixes the parameter again;
    - and so on, up to its limit of 1000 passes, with `delta` unchanged.

    About 4% of the columns of Q do this near convergence, at about 200
    times the cost of a typical column. The loop now stops as soon as it is
    about to fix the parameter it freed in the previous pass, after a step
    of exactly 0. Fits are bit-identical, and `update_q` is about 10%
    faster near convergence.

### Memory

11. **The whole loop in C++, with buffers allocated once.** The fitting loop
    was first written in R, calling the C++ kernels for each pass. Every call
    returned a new R matrix, and the matrices it replaced (each a copy of P)
    accumulated until R's garbage collector ran. The peak memory reached
    about 30-34 copies of the parameters. In C++, every estimate lives in a
    buffer allocated at the start of the fit:
    - `x0`, `x1`, `x2` and the quasi-Newton point;
    - the `2 n_qn` secant vectors.

    That is 10 copies of K (I + J) values with the default `n_qn = 3`.
    Accepting an estimate swaps buffers in constant time, where Julia copies
    the parameters (`d.x .= d.x_tmp`). Less memory also means fewer cache
    misses and page faults on large datasets.

### Data access

12. **Reading the backing file of the `gen_tibble` directly.** Julia reads a
    PLINK bed file through SnpArrays. Running ADMIXTURE from R
    (`tidypopgen::gt_admixture()`) also requires exporting the `gen_tibble`
    to a bed file first. `gt_radmix()` instead reads the memory-mapped
    backing file of the `gen_tibble`, through the row and column indices of
    the current selection, so filtering individuals or selecting loci costs
    nothing. When the rows are in their original order, the 64 individuals
    of a slice are adjacent bytes.

13. **Cross-validation without stored masks, and with warm starts.**
    - The fold of a genotype is not stored: it is computed on the fly from a
      hash (the splitmix64 finaliser) of the seed and the position of the
      genotype. This needs no memory (a stored mask would be as large as the
      genotypes), is thread safe, and gives the same folds for the same seed.
    - The hash runs for every genotype in every pass of a fold, so it is
      kept cheap: the remainder by the number of folds uses a compile-time
      constant for 2 to 10 folds (a multiplication instead of a 64-bit
      division), and hidden genotypes are dropped without branches (they are
      spread at random, so a branch was mispredicted about once in
      `n_folds`). The folds and the results are unchanged; a pass with a
      fold hidden went from about 1.9 to 1.5 times the cost of a plain pass
      (quick and full HGDP data of the benchmark), and a fit with 5-fold
      cross-validation became about 30% faster.
    - Each fold is fitted starting from the estimates for the full data,
      without EM priming. They are close to the optimum without the fold, so
      each fold converges in a few iterations (fastmixture also starts its
      folds from the full-data estimates).

## Other differences from OpenADMIXTURE.jl

These are not about speed. Some extend what the Julia code does, and the
others are numerical safeguards that do not change the results when the
computations are well behaved.

- **Pseudohaploids.** The number of trials is `n_i` (2, or 1 for
  pseudohaploid individuals) instead of the constant 2. Julia assumes that
  every individual is diploid.
- **EM update of Q with missing data.** The EM update of Q is normalised by
  `n_i J_i` (the allele copies actually observed in individual `i`), which is
  also the sum of its numerators over k. Julia divides by `2J`, which equals
  that sum only when nothing is missing, so with missing data Julia's Q
  leaves the simplex.
- **Cross-validation**, which the Julia code does not offer (see above).
- **The penalty `mu` of the Q problems is floored at 0.** Julia does not
  floor it, so a negative `mu` could make the tableau indefinite (Lange's
  original code, kept commented out in `qp.jl`, does floor it). In practice
  Julia's `mu` is 0 anyway, because `mul!(tmp_k, v, tmp_k)` aliases its input
  and output, and BLAS `gemv` zeroes the output first. On the simplex `mu`
  does not change the solution.
- **No sweep on a vanishing pivot.** Julia sweeps a parameter back into the
  active set without checking its pivot, which could divide by zero and fill
  the step with Inf/NaN. The same `1e-10` guard as for the inverse sweep is
  applied.
- **Non-finite Newton steps are discarded.** If the step of a column is not
  finite (e.g. a genotype that is impossible under the current estimates),
  the column keeps its value for that sweep, instead of propagating NaN to
  the log-likelihood.
- **Projections tolerate non-finite values.** A column of Q with a
  non-finite value is reset to the centre of the simplex, and a NaN in P to
  0.5. This should never trigger, given the guards above.
- **P is the frequency of the alternate allele.** SnpArrays (and ADMIXTURE)
  report the frequency of the second allele of the bim file, which for a
  `gen_tibble` is the reference allele, so their P is `1 - P`.
- **Not ported:** GPU (CUDA) support, the sparse K-means selection of
  ancestry informative markers (SKFR), and the `:LBQN` variant of the
  quasi-Newton step (only the default `:ZAL` is implemented).

## Performance in practice

On the `EUR_subset` data of SnpArrays.jl (379 individuals, 54,051 loci), on a
laptop with an Intel i7-10510U
(`data-raw/radmix/benchmark_admixture.csv`), times are in seconds,
with the number of iterations in brackets:

| K | ADMIXTURE 1.3.0, 1 core | `gt_radmix()`, 1 core | `gt_radmix()`, 1 core, ADMIXTURE's threshold | `gt_radmix()`, 4 cores |
|---|---|---|---|---|
| 2 | 37.3 (20) | 13.2 (14) | 16.8 (17) | 3.8 (14) |
| 3 | 109.3 (39) | 23.7 (19) | 39.1 (29) | 8.5 (19) |
| 4 | 121.5 (33) | 32.7 (22) | 46.0 (31) | 12.6 (22) |

With its default threshold (relative change in log-likelihood below 1e-7),
`gt_radmix()` stops earlier than ADMIXTURE (absolute change below
1e-4). The column "ADMIXTURE's threshold" uses ADMIXTURE's criterion, for a
like-for-like comparison. On this laptop the 4-core speed-up is limited by
the CPU, which slows down after 10-20 s with 4 busy cores. The benchmark on
the HGDP (`vignettes/articles/admixture_benchmark.Rmd`) also compares
`gt_radmix()` with fastmixture, on 1 to 8 cores.

## Possible further improvements

From `data-raw/radmix/Notes_to_do`:

- store the secant pairs in single precision, saving 3 of the 10 copies of
  the parameters;
- keep `1 - P` in a buffer reused across kernel calls, rather than computing
  it in each call;
- tune the tile sizes and the OpenMP schedule, which have not been tuned
  since the port;
- for very large datasets, fit on a subset of loci and then estimate P for
  the others in chunks, to bound memory.

## References

Alexander DH, Novembre J, Lange K (2009) Fast model-based estimation of
ancestry in unrelated individuals. *Genome Research*, 19, 1655-1664.

Alexander DH, Lange K (2011) Enhancements to the ADMIXTURE algorithm for
individual ancestry estimation. *BMC Bioinformatics*, 12, 246.

Ko S, Chu BB, Peterson D, Okenwa C, Papp JC, Alexander DH, Sobel EM, Zhou H,
Lange KL (2023) Unsupervised discovery of ancestry-informative markers and
genetic admixture proportions in biobank-scale datasets. *The American
Journal of Human Genetics*, 110, 314-325.

Zhou H, Alexander D, Lange K (2011) A quasi-Newton acceleration for
high-dimensional optimization algorithms. *Statistics and Computing*, 21,
261-273.

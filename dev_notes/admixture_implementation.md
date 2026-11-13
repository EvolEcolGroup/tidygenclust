# Implementation of `gt_radmix()`

`gt_radmix()` estimates ancestry proportions (Q) and ancestral allele
frequencies (P) by maximum likelihood under the model of ADMIXTURE
(Alexander, Novembre and Lange 2009). It uses the optimisation strategy of
ADMIXTURE: block relaxation, in which Q and P are updated in turn by
constrained Newton steps, accelerated by the quasi-Newton scheme of Zhou,
Alexander and Lange (2011), after a few EM steps that prime a random start.
It offers the cross-validation of ADMIXTURE to choose K (Alexander and
Lange 2011) and, as an option, a start from a truncated singular value
decomposition (SVD) of the genotypes, as in fastmixture and ADAMIXTURE.

The computations are in C++, multithreaded with OpenMP, and read the
genotypes directly from the backing file of a `gen_tibble`. This document
derives the algorithm from first principles, describes how it is
implemented, and ends with the origins of the code.

## Where the code is

| File | Role |
|---|---|
| `R/gt_radmix.R` | User interface: checks the arguments, prepares the data once, loops over K and runs, assembles the `gt_admix` object |
| `R/radmix_data.R` | Builds the description of the genotypes passed to C++ (backing file, selected rows and columns, decoding table, ploidy) and the cross-validation masks |
| `R/radmix_svd.R` | The SVD start (`init = "svd"`): the factors of the SVD, and the best of several ALS solutions |
| `R/radmix_fit.R` | Thin wrapper around the C++ fitting loop |
| `R/radmix_cv.R` | Cross-validation |
| `src/radmix_fit.cpp` | The fitting loop: random start, EM priming, block updates, quasi-Newton steps, convergence |
| `src/radmix_kernels.cpp` | The passes over the genotypes: log-likelihood, deviance, EM step, gradients and Hessians, block updates |
| `src/radmix_qp_sweep.cpp` | The solver of the quadratic program of each column (sweep operator, active set) |
| `src/radmix_quasi_newton.cpp` | The quasi-Newton step |
| `src/radmix_projections.cpp` | Projections onto the simplex (Q) and the box (P) |
| `src/radmix_als.cpp` | Alternating least squares on the factors of the SVD (`init = "svd"`) |
| `src/radmix_geno_access.h` | Read-only view of the genotypes in the backing file, with the cross-validation masks |
| `src/radmix_exports.cpp` | The R interface of the building blocks (used by the tests, and by R for the log-likelihood, the deviance and the projections) |

The C++ files are commented in detail; this document gives the overall
picture and the reasons for the main choices.

## The model

For individual `i` and locus `j`, the genotype `g_ij` is the number of copies
of the alternate allele, out of `n_i = 2` for diploids or `n_i = 1` for
pseudohaploids. It is binomial,

```
g_ij ~ Binomial(n_i, qp_ij),    qp_ij = sum_k q_ki p_kj,
```

where `q_ki` is the proportion of the genome of individual `i` that comes from
ancestral population `k`, and `p_kj` the frequency of the alternate allele at
locus `j` in population `k`. Writing `h_ij = n_i - g_ij` for the copies of the
reference allele, the log-likelihood over the observed genotypes (missing
genotypes do not contribute, and the binomial coefficients, which do not
depend on Q or P, are dropped) is

```
L(Q, P) = sum_ij l_ij(qp_ij),    l_ij(x) = g_ij log(x) + h_ij log(1 - x).
```

It is maximised subject to each column of Q lying on the simplex
(`q_ki >= 0`, `sum_k q_ki = 1`) and each `p_kj` lying in `[0, 1]`. All
estimates are kept within `[1e-5, 1 - 1e-5]`, as in ADMIXTURE, so that `qp`
and `1 - qp` are strictly positive and the log-likelihood is always finite.

Internally, Q is stored as a K x I matrix and P as a K x J matrix,
column-major: the K values of one individual (or one locus), which are
updated together, are contiguous. Q and P are returned to R as I x K and
J x K matrices.

### Derivatives

Everything the algorithm needs comes from three per-genotype quantities. For
one genotype (dropping the indices), let

```
a  = g / qp,                     b = h / (1 - qp),
d1 = l'(qp)  = a - b,            d2 = -l''(qp) = g / qp^2 + h / (1 - qp)^2.
```

`qp_ij` is linear in the column `q_.i` (with coefficients `p_.j`) and in the
column `p_.j` (with coefficients `q_.i`), so by the chain rule

```
dL / dq_ki = sum_j d1_ij p_kj,      -d2L / dq_ki dq_k'i = sum_j d2_ij p_kj p_k'j,
dL / dp_kj = sum_i d1_ij q_ki,      -d2L / dp_kj dp_k'j = sum_i d2_ij q_ki q_k'i.
```

Three properties shape the algorithm:

- **Separability.** Given P, `L` is a sum of one term per individual, and
  given Q, of one term per locus. So with one block fixed, each column of the
  other block can be optimised on its own, with its own K x K Hessian.
- **Concavity.** `d2 >= 0`, so each negative Hessian is a sum of terms
  `d2 v v'` and is positive semi-definite: `L` is concave in each column,
  and a constrained Newton step for a column is a convex quadratic program.
- **Symmetry.** The Q and P sides have the same form, with the roles of the
  individuals and the loci swapped. The implementation shares the
  per-genotype coefficients between them; only the direction of the sums
  differs.

## The algorithm, step by step

### 1. Preparing the data (once per call)

`radmix_prepare()` gathers what the C++ code needs, without copying the
genotypes:

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
(`radmix_geno_access.h`) reads genotype `(i, j)` as a memory read, a table
lookup and a multiplication.

### 2. Runs

For each value of K and each run, `gt_radmix()` sets the seed (if given) and
calls the fitting loop (`radmix_fit_cpp()`). Random starting values are drawn
from R's random number generator, so `set.seed()` makes fits reproducible;
when `seed` is given, the user's random number state is restored on exit.

### 3. Starting values

There are three ways to start a fit.

**Random start (the default).** Q and P are drawn from uniform
distributions, and each column of Q is rescaled to sum to one. This is far
from the optimum, so the fit first runs `n_em_iter` (5 by default) EM steps
(section 4), which move it cheaply into a sensible region of the parameter
space. Random starts explore the modes of the likelihood, which is what the
multi-run analyses of `gt_clumppling()` need.

**Starting values given by the user** (`init_q`, `init_p`). They replace
the corresponding block after the EM priming, so that when only one block is
given the other still starts from a sensible value; with both, EM is
skipped. With `fix_q` or `fix_p`, the given block is kept fixed throughout
(e.g. to project new individuals onto the allele frequencies of a reference
panel).

**SVD start** (`init = "svd"`). Under the model, the expected allele
frequencies `Pi = P' Q` (loci x individuals) form a matrix of rank at most K.
A truncated SVD of the centred and scaled genotypes gives a low-rank
approximation of it, and alternating least squares (ALS) then looks for
non-negative factors of that approximation, which are a good start:

- the SVD (`radmix_svd_factors()`) is computed once per call, with
  `max(k) + 5` components, by `tidypopgen::gt_pca_partialSVD()`, with each
  locus centred by `2 f` and scaled by `sqrt(2 f (1 - f))`, `f` being the
  frequency of its alternate allele (monomorphic loci are left unscaled, and
  contribute nothing). It needs a `gen_tibble` without missing genotypes, or
  with imputed ones; the imputed genotypes are only used for the start. The
  factors give `Pi ~ Za Ua'`, with K + 5 components plus a column of allele
  frequencies (in `Za`) and of ones (in `Ua`);
- ALS (`radmix_als.cpp`) alternates a Q step, which for each individual
  solves a least squares problem bounded in `[0, 1]`, weighting each locus by
  the inverse of its binomial variance (`1 / max(f (1 - f), 1e-3)`,
  normalised to mean one), and a P step, which does the same for each locus
  without weights, from Q extrapolated along its last change (momentum 0.8).
  Each step only involves products with the factors, never the genotypes.
  The bounded problems share their matrix, so their unconstrained solutions
  are computed at once and used where they lie inside the box; the others
  are solved by coordinate descent. Q is then clipped to
  `[1e-5, 1 - 1e-5]` and normalised, and P clipped. ALS stops when the root
  mean square change of Q falls below 1e-4;
- three ALS solutions, from random allele frequencies, are computed, and the
  one with the highest log-likelihood (one pass over the genotypes each) is
  the start. The fit then skips EM.

K + 5 components rather than K, and the best of three ALS solutions, were
chosen in `data-raw/radmix/explore_svd_start_refined.R`: they avoid the poor
starts that the plain version gave on weakly structured data. An SVD start
typically halves the number of iterations, but it is nearly deterministic,
so repeated runs explore few modes: it is an option, for a fast single fit,
not the default. R's BLAS is limited to one thread while ALS runs: its
products are tall and skinny, and BLAS threads only competed with the OpenMP
threads.

### 4. The EM step

Each allele copy of individual `i` at locus `j` comes from one of the K
ancestral populations. Given an alternate allele, the probability that it
came from population `k` is `q_ki p_kj / qp_ij`; given a reference allele, it
is `q_ki (1 - p_kj) / (1 - qp_ij)`. The E step assigns to population `k`

```
g q_ki p_kj / qp + h q_ki (1 - p_kj) / (1 - qp) = q_ki [a p_kj + b (1 - p_kj)]
```

expected allele copies, of which `q_ki p_kj a` are alternate alleles. The M
step sets `q_ki` to the share of the allele copies of individual `i` assigned
to `k`, and `p_kj` to the share of alternate alleles among the copies at
locus `j` assigned to `k`:

```
q_new_ki = q_ki sum_j [a p_kj + b (1 - p_kj)] / (n_i J_i)
p_new_kj = A / (A + B),   A = p_kj sum_i a q_ki,   B = (1 - p_kj) sum_i b q_ki
```

This is the EM algorithm of FRAPPE (Tang et al. 2005). Summed over `k`, the
copies assigned at one locus are `a qp + b (1 - qp) = g + h = n_i`, so the
denominator `n_i J_i` (`J_i` being the number of loci observed in individual
`i`) is the sum of the numerators over `k`: normalising by that sum keeps Q
on the simplex exactly, with any pattern of missing data. Both blocks are
updated from the same Q and P, in two passes over the genotypes (one over
blocks of individuals for Q, one over blocks of loci for P). EM never
decreases the likelihood, but converges slowly, which is why it is only used
to prime random starts.

### 5. Block updates

Each iteration of the main loop applies the block-update map `F` twice: one
application updates all of Q given P, then all of P given the new Q. By
separability, every column is updated on its own, by a constrained Newton
step.

**The quadratic program.** For a column `x` (the K proportions of an
individual, or the K frequencies of a locus), with gradient `grad` and
negative Hessian `H` from the formulas above, the step `delta` maximises the
quadratic approximation of the log-likelihood,

```
maximise  grad' delta - delta' H delta / 2
subject to  0 <= x + delta <= 1    (and sum_k delta_k = 1 - sum_k x_k for Q).
```

Without the constraints the solution is the Newton step `H^-1 grad`.

**The solver** (`radmix_qp_sweep.cpp`). K is small (typically below 20)
and there are I + J problems per sweep, so each is solved exactly by a small
dense active-set method based on the sweep operator (Lange 2010), with
preallocated scratch memory and no library calls:

- a symmetric tableau bundles H, the gradient and, for Q, the sum
  constraint. Sweeping a diagonal entry brings a parameter into the Newton
  system (Gauss-Jordan elimination on one pivot, keeping the matrix
  symmetric), and the inverse sweep takes it out again, each at O(K^2) cost;
- sweeping all the parameters, and then the sum constraint, gives the Newton
  step of the problem with the equality constraint (by Lagrange
  multipliers) but without the bounds;
- the solver walks along that step until a parameter reaches a bound, fixes
  that parameter (inverse sweep), recomputes the step of the others, and
  frees again any fixed parameter whose gradient points back into the box,
  until the KKT conditions hold;
- a pivot that is not positive, or that has shrunk below 1e-8 of its
  original value while sweeping, means that H is not numerically positive
  definite: the step is then 0.

For Q, sweeping the parameters requires the K x K block to be positive
definite, whereas the constrained problem only needs H to be positive
definite on the subspace `sum_k delta_k = 0`. Since H is positive
semi-definite, `H + mu 11'` is positive definite for any `mu > 0` as soon as H
is positive definite on that subspace (for `x` with `1'x != 0`,
`x'(H + mu 11')x >= mu (1'x)^2 > 0`), and the shift does not change the
solution, as every feasible step has the same `1' delta`. The tableau uses
`mu = trace(H) / K^2`, which adds the mean eigenvalue of H along the direction
of the ones and so keeps the matrix as well conditioned as H.

A solution at a vertex of the simplex (every parameter of a column of Q on a
bound) makes the active set cycle: a parameter is freed, cannot move, is
fixed again, and so on. This is detected (a parameter fixed again right
after being freed, after a step of exactly 0), and the solver stops at once
with the right step. Its return value (0 in that case, or on failure) is
therefore not a failure signal, and the updates use the step it found
regardless.

**Projection.** The new column `x + delta` is projected onto
`{sum = 1, q >= 1e-5}` for Q (by the sort-and-threshold algorithm for the
Euclidean projection onto a simplex), or clamped to `[1e-5, 1 - 1e-5]` for P
(`radmix_projections.cpp`). This removes rounding errors and keeps every
estimate strictly inside `(0, 1)`.

**Impossible genotypes.** If a genotype is impossible under the current
estimates (`qp` of exactly 0 with `g > 0`, or 1 with `h > 0`), its
coefficients are infinite, and so is the step of its column. Such a step is
discarded: the column keeps its value for this update, the other columns
are updated, and it gets another chance at the next update, once the other
block has moved. Without this, the NaN would spread to the whole column
and, through the next updates, to the log-likelihood. This is a safeguard
only: the estimates stay within `[1e-5, 1 - 1e-5]`, and the exact 0 or 1
that an EM step can give at a monomorphic locus agrees with its genotypes.

### 6. Quasi-Newton acceleration

The maximum likelihood estimate is a fixed point of the block-update map,
`x = F(x)`, with `x = (Q, P)`. Each iteration gives `x0` (the current
estimate), `x1 = F(x0)` and `x2 = F(x1)`, hence a secant pair `u = x1 - x0`,
`v = x2 - x1`, which tells how `F` changes along `u`. With the last `n_qn`
pairs (3 by default) as the columns of U and V, the smallest matrix `M` that
satisfies the secant conditions `M U = V` is `M = V (U'U)^-1 U'`. Solving
`x = F(x0) + M (x - x0)` for the fixed point of the linear approximation of
`F`, by the Woodbury identity, gives the quasi-Newton point (Zhou, Alexander
and Lange 2011)

```
x_qn = x1 - V (U'(U - V))^-1 U'(x0 - x1).
```

Only an `n_qn` x `n_qn` system is solved, so the step costs two passes over
the parameters (`radmix_quasi_newton.cpp`). A system that is singular to
machine precision (secant pairs become nearly collinear close to
convergence), or that is not finite, is rejected, and the step is skipped.
The point ignores the constraints, so it is projected onto the feasible set.

**Accept or reject.** The quasi-Newton point is kept if its log-likelihood
is higher than that of `x0`; otherwise the iteration ends at `x2`. The
comparison is with `x0` rather than with `x2`: a point slightly worse than
`x2` still makes progress, and tends to pay off in the following
iterations. It also means that `L(x2)` is only needed when the point is
rejected; the point is accepted in most iterations.

**Evaluating a point and starting from it in one pass.** The first block
update of an iteration, the update of Q given `x0 = (q, p)`, forms `qp` for
every genotype, which is all that the log-likelihood of `x0` needs. So
whenever the loop evaluates a point from which the next iteration may start
(the start, the quasi-Newton point, and `x2` after a rejection), it also
computes the update of Q from that point in the same pass
(`update_q_loglik()`), and the next iteration starts from it. An iteration
then costs four passes over the genotypes (three block updates and the
fused pass at the quasi-Newton point) instead of five, and five instead of
six after a rejection. The update is identical to that of a separate pass,
and the log-likelihood differs from that of `loglik()` only by rounding (the
same terms, added up in another order, still independent of the number of
threads). The update computed with the last log-likelihood of a fit (at
convergence) or with a rejected point is wasted; it is not computed after
the last iteration allowed (`max_iter`), nor when Q is fixed.

### 7. Convergence

The loop stops when the relative change in log-likelihood over an iteration
falls below `rtol` (1e-8 by default), or after `max_iter` iterations (1000).
A NaN change (e.g. with no observed genotype) never counts as convergence.

### 8. Cross-validation (optional)

With `crossval = TRUE`, the cross-validation error of ADMIXTURE (Alexander and
Lange 2011) is computed after the main fit, so the main estimates are the
same with or without it. The observed genotypes are split into `n_folds`
folds (5 by default). For each fold:

1. The model is refitted without the genotypes of that fold, through a
   "training" view of the genotypes, which hides them from every pass
   exactly as if they were missing. The fit starts from the estimates for
   the full data (a warm start, as in fastmixture), without EM: they are
   close to the optimum without the fold, so each fold converges in a few
   iterations. The convergence criteria are those of the main fit, unless
   `cv_max_iter` or `cv_rtol` shorten the refits (faster, but the errors
   are then biased downwards, as the start has seen the held-out
   genotypes; measured in `data-raw/radmix/newspeed_algo/`).
2. The genotypes of the fold are predicted by `mu = n qp`, through a
   "held-out" view that shows only them.
3. The binomial deviance of those genotypes is summed:
   ```
   d(g, mu) = g log(g / mu) + (n - g) log((n - g) / (n - mu))
   ```
   with `0 log 0 = 0` (ADMIXTURE gives this formula for diploids).

The error is the deviance averaged over all the held-out genotypes of all the
folds. The fold of a genotype is computed from a hash (the splitmix64
finaliser) of a seed and the position of the genotype, modulo the number of
folds. This is thread safe and gives the same folds for the same seed;
folds are of roughly equal size, and genotypes that are already missing
belong to no fold. The fit of a fold evaluates the hash once per genotype,
before its first pass, and keeps which genotypes are visible as a bit mask
(one bit per genotype, 1/8 of the size of the genotypes; e.g. 18 MB for the
HGDP of 1043 x 130,943, 87 MB for 1043 x 642,794), so that its passes read
one 64-bit word per slice of 64 individuals instead of hashing every
genotype. Beyond 1 GiB of mask (8.6e9 genotypes), the passes hash instead.
The results are the same either way; the held-out view, used once per fold,
always hashes. The seed is drawn from the random number
generator after the main fit, whose state is then restored. ADMIXTURE's own
random partition cannot be reproduced (its source code is not available), so
the errors agree with ADMIXTURE's statistically, not exactly.

## How the passes over the genotypes are organised

Every pass reads the whole genotype matrix, so the number of passes and
their cost determine the run time on large data. A fit makes `n_em_iter`
EM steps (two passes each), one pass for the starting log-likelihood (and
the first update of Q), and four passes per iteration (five after a
rejected quasi-Newton step), see "Evaluating a point and starting from it
in one pass" above.

### Data access

The backing file is memory-mapped, and read through the row and column
indices of the current selection, so filtering individuals or selecting loci
costs nothing, and nothing is exported to a PLINK file. The genotypes are
stored individuals x loci, column-major, so the genotypes of consecutive
individuals at one locus are adjacent bytes. When the individuals of a slice
(below) are adjacent rows of the backing file and decoded alike (no
pseudohaploids), their bytes are read directly and decoded with two tables
giving `g` and `h` (0 for missing); otherwise each genotype goes through the
row index, the decoding table and the ploidy scale. With a cross-validation
view, hidden genotypes are zeroed by multiplying rather than by a branch,
which would be mispredicted (they are spread at random). The fit of a fold
reads their visibility from its bit mask (see section 8), one word per
slice; without the mask, the fold hash runs for every genotype in every
pass, and is kept cheap (the remainder by the number of folds uses
compile-time constants for 2 to 10 folds, a multiplication instead of a
64-bit division), but such a pass costs about 1.5 times a plain pass.

### Slices and coefficients

The unit of work is a *slice*: one locus and up to 64 individuals. For each
slice, `qp` is computed for all its individuals with unit-stride loops over
a transposed copy of Q (made once per pass); the genotypes are decoded in a
separate loop (the only one with table lookups); and the coefficients (`a`
and `b` for EM, `d1` and `d2` for the Newton updates) are computed without
branches. A missing genotype has `g = h = 0`, so it drops out of every sum
without a test, and terms with a zero count are computed as 0 (by dividing
the count by 1 instead of `qp`), so that `qp` of exactly 0 or 1 does not give
`0/0`. These loops vectorise. `d1` and `d2` take two divisions rather than
four (`rx = 1 / qp`, `a = g rx`, `a rx = g / qp^2`).

### Accumulation

- **Q side** (one sum per individual): a block of 64 individuals is
  streamed over the loci, with one accumulator vector of length 64 per
  gradient entry and per entry of the Hessians, so that each update is a
  unit-stride loop over the individuals of the slice. The products
  `p_kj p_k'j` do not depend on the individual, so they are computed once
  per locus.
- **P side** (one sum per locus): for a block of 64 loci, each sum is a dot
  product over the individuals of the slice, while the 64 x K corner of Q
  stays in the L1 cache.
- **Two loci at a time.** The coefficients of two loci are computed first,
  and each accumulator is then read and written once for both, which halves
  the memory traffic of the accumulators (the accumulation is limited by
  loads and stores rather than arithmetic; 6-14% faster). A range with an
  odd number of loci is padded with a locus whose coefficients are zero.
- **Half Hessians.** Only the K (K + 1) / 2 entries of each symmetric
  Hessian with `k' >= k` are accumulated; the others are filled in at the
  end.
- **Hessians never stored.** Each thread accumulates the gradients and
  Hessians of one block of individuals (or loci) into thread-local buffers,
  and solves their quadratic programs straight away. The array of the
  Hessians of all the loci (K times the size of P) is never allocated, and
  each Hessian is used while it is in cache.

### The log-likelihood

For a hard call, the contribution of a genotype is `log(z)` with
`z = qp^g (1 - qp)^h`, which can be computed without logarithms or branches
(a polynomial in `g` and `h`). Within a slice of 64 individuals, the `z` of
8 genotypes are multiplied together, in 8 interleaved products that
vectorise, before taking a single logarithm: 8 logarithms per 64 genotypes
instead of up to 128. Each factor is at least 1e-30, so a product of 8 cannot
underflow; dosages and the rare `z` below 1e-30 (including `z = 0`, which
correctly gives `-Inf`) use the direct formula.

### Parallelism and determinism

Work is split into tiles of 64 individuals x 64 loci. Q-side work is
parallelised over blocks of individuals and P-side work over blocks of loci,
so each thread writes its own output columns: there are no races and no
reductions across threads. The scalar sums (log-likelihood, deviance) keep
one partial sum per block of loci, added up in order, and the long inner
products of the quasi-Newton step are split into chunks of fixed size, added
up in order.

With few individuals there are few blocks of 64 (379 individuals give 6, 50
individuals a single one), which cannot keep many threads busy, and
population genetic data often have tens to hundreds of individuals and many
loci. When there are fewer than 32 blocks, the Q side therefore also splits
the loci into segments (of at least 256 loci), and each (block, segment) pair
is a separate task; the partial sums of the segments are added up in order.

Block and segment boundaries depend only on the dimensions of the data,
never on the number of threads, and each block is processed start to finish
by one thread, so the estimates are bit-identical for any `n_cores`.

### Instruction sets

R compiles packages for the baseline x86-64 instruction set (SSE2, no fused
multiply-add). On x86-64 Linux with GCC 12 or later, the two accumulation
functions, which hold the hot loops, are also compiled for the x86-64-v3
level (AVX2 and FMA), and the version that suits the processor is chosen
when the package is loaded; this makes the block updates and EM steps
15-30% faster. The two versions can differ by rounding, so results are
reproducible on a given machine but not bit-identical across machines;
compiling with `-DRADMIX_NO_MULTIVERSION` builds only the baseline version.
The log-likelihood is not multiversioned (its AVX2 version was slower).

### Memory

The whole fitting loop runs in C++, and every estimate lives in a buffer
allocated at the start of the fit: `x0`, `x1`, `x2` and the `2 n_qn` secant
vectors, i.e. `3 + 2 n_qn` copies of the K (I + J) parameters (9 with the
default `n_qn = 3`). The quasi-Newton point overwrites `x0`, which is no
longer needed once the secant pair is stored, and accepting an estimate
swaps buffers rather than copying them. Each pass also makes a transposed
copy of Q (K x I), and each thread small buffers of the size of a tile.

## Numerical safeguards

None of these changes the results when the computations are well behaved:

- coefficients with a zero count are 0 even when `qp` is exactly 0 or 1;
- a Newton step that is not finite is discarded (section 5);
- the quadratic program gives up (step 0) on a pivot that is not positive or
  has shrunk below 1e-8 of its original value, never sweeps a vanishing pivot
  back in, and stops when the active set cycles at a vertex;
- the projections reset a column of Q with a non-finite value to the centre
  of the simplex, and a NaN in P to 0.5 (this should never happen, given the
  guards above);
- the quasi-Newton step is rejected if its small system is singular to
  machine precision or not finite;
- EM can return `p` of exactly 0 or 1 at monomorphic loci; the other passes
  handle it through the guards on `qp` of exactly 0 or 1.

## Validation

The test suite checks:

- every pass against a plain R implementation on dense matrices
  (`tests/testthat/helper-radmix-reference.R`), with missing data,
  pseudohaploids, unsorted rows and columns, tiles at the edges and the
  segmented Q side, and the gradients and Hessians against numerical
  derivatives;
- the quadratic program against `quadprog`, and the quasi-Newton step against
  its formula;
- that the results are identical for any number of threads;
- the estimates, log-likelihoods and cross-validation errors against stored
  results of ADMIXTURE 1.3.0 on the `EUR_subset` data of SnpArrays.jl
  (`tests/testthat/test_radmix_vs_admixture.R`);
- the log-likelihood, EM step, Newton updates and full path against stored
  results of OpenADMIXTURE.jl, an independent implementation of the same
  model and of a closely related algorithm
  (`tests/testthat/test_radmix_vs_julia.R`).

`data-raw/radmix/README.md` describes how the stored results were made.

## Performance

`vignettes/articles/admixture_benchmark.Rmd` compares `gt_radmix()` with
ADMIXTURE, fastmixture, ADAMIXTURE and OpenADMIXTURE.jl on the HGDP data,
on 1 to 8 cores, with and without cross-validation and with random and SVD
starts. `dev_notes/radmix_performance.md` records the changes that made the code
faster, with their measurements, the profile and the multi-core scaling,
the ideas that were tried and rejected, and the leads for further speed;
its last section explains how to read the full reports, scripts and
results of that work, archived under the git tag
`archive/radmix-perf-verify`.

## Possible further improvements

For speed, see "Leads not pursued yet" in `dev_notes/radmix_performance.md`. For
memory:

- Store the secant pairs in single precision (saving 3 of the 9 copies of
  the parameters).
- For very large data, fit on a subset of loci and then estimate P for the
  others, one locus at a time, to bound memory.

## Origins of the code

The C++ code of `gt_radmix()` started as a port of OpenADMIXTURE.jl (Ko et
al. 2023; commit `401ac5c4056f0f6f2a7ea6ffc3af69c2b1e56fdc`), a Julia
implementation of the ADMIXTURE model and algorithm. It has since been
rewritten extensively, and its computations now differ from those of
OpenADMIXTURE.jl in many places. The main changes were:

- the derivatives of the Q block are those of the log-likelihood itself
  (OpenADMIXTURE.jl differentiates a homogeneous form of it, which gives the
  same Newton steps but twice the work for the Hessians), so that both
  blocks share the same per-genotype coefficients; and the shift `mu` of the
  quadratic programs of Q is `trace(H) / K^2`;
- the passes over the genotypes were reorganised: quadratic programs solved
  block by block, without storing the Hessians; branch-free, vectorised
  coefficients; half Hessians; two loci at a time; batched logarithms;
  two-dimensional blocking of the Q side; deterministic sums; and
  multiversioning;
- the fitting loop runs entirely in C++, with a fixed set of buffers, and
  computes the log-likelihood of `x2` only when the quasi-Newton point is
  rejected (and that of the EM steps only when it is reported);
- the EM step normalises Q by the observed allele copies of each individual,
  which keeps it on the simplex with missing data;
- pseudohaploid individuals, the cross-validation of ADMIXTURE, the SVD
  start and the direct reading of the backing file of a `gen_tibble` were
  added;
- numerical safeguards were added (section "Numerical safeguards"), and the
  active-set solver stops when it cycles at a vertex;
- the selection of ancestry informative markers by sparse K-means, which is
  a central part of OpenADMIXTURE, and its GPU support were not ported.

The part that remains closest to the original is the active-set solver of
the quadratic programs, which follows K. Lange's sweep-operator code in
OpenADMIXTURE.jl; its MIT licence is reproduced in `inst/COPYRIGHTS`.

## References

Alexander DH, Novembre J, Lange K (2009) Fast model-based estimation of
ancestry in unrelated individuals. *Genome Research*, 19, 1655-1664.

Alexander DH, Lange K (2011) Enhancements to the ADMIXTURE algorithm for
individual ancestry estimation. *BMC Bioinformatics*, 12, 246.

Ko S, Chu BB, Peterson D, Okenwa C, Papp JC, Alexander DH, Sobel EM, Zhou H,
Lange KL (2023) Unsupervised discovery of ancestry-informative markers and
genetic admixture proportions in biobank-scale datasets. *The American
Journal of Human Genetics*, 110, 314-325.

Lange K (2010) *Numerical Analysis for Statisticians*, 2nd edition.
Springer.

Tang H, Peng J, Wang P, Risch NJ (2005) Estimation of individual admixture:
analytical and study design considerations. *Genetic Epidemiology*, 28,
289-301.

Zhou H, Alexander D, Lange K (2011) A quasi-Newton acceleration for
high-dimensional optimization algorithms. *Statistics and Computing*, 21,
261-273.

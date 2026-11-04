# Simplification plan for tidygenclust

Simplifications that would make the code shorter or clearer without losing
performance. The plan comes from the code review of October 2026
(`data-raw/code_review.md`) and was updated after the bug fixes that
followed it. Its `gt_radmix()` items were done on 3 November 2026, together
with a rewrite of the Q-side derivatives (see "Done on 3 November 2026"
below); the other items have not been implemented yet.

"Identical" means bit-identical outputs (`gt_radmix()`) or identical results
on the 60 scenarios of `data-raw/clumppling/validate_vs_python.R`
(clumppling). Line counts are approximate.

## Tier 1: prototyped, identical results and no loss of speed

- **clumppling alignment**: `clumppling_reorder_q_across_k()` as a single
  loop over K (about 20 lines); `clumppling_align_q_wrt_p()` written with
  `order()`/`duplicated()` without merging, and with
  `clumppling_move_columns()` with merging (about 25 lines).

## Tier 2: refactoring off the hot path

- **Plots**: `plot_major_modes()`, `plot_all_modes()` and
  `plot_modes_within_k()` (`R/autoplot_gt_clumppling.R`) differ only in the
  modes they select (about 50 lines); `R/plot_modes.R` repeats its inset and
  rescaling code (about 70 lines).
- **Argument checks**: `gt_clumppling()` checks its logical options with
  its own loop; it could use `check_flag()` (from `R/gt_radmix.R`,
  which `gt_fastmixture()` already shares; its messages say "should be"
  rather than "must be").
- **Tests**: shared helpers for pseudohaploid gen_tibbles, starting values
  and cluster matching in the tests of `gt_radmix()` (about 50 lines).
- **Installer**: resolve the yaml file of each environment with one helper,
  and always leave fastmixture out of the yaml (it is reinstalled from
  github in any case).
- **`R/utils.R`**: `factor(x, levels = unique(x))` instead of fetching the
  internal `fct_inorder_base()` of tidypopgen.
- **`clumppling_run()`** returns several elements that nothing uses (about
  10 lines).

## Not recommended

- Reading genotypes through bigstatsr's `SubBMCode256Acc`: about 19% slower
  than the custom genotype view.
- Simpler sums in clumppling (`colSums()` for the representatives, the
  within-K cost from the assignment objective): they change which of
  several near-identical replicates is the representative, breaking the
  agreement with python.
- `igraph::cluster_louvain()` instead of the louvain port: a different
  algorithm, so modes would differ from python in ambiguous cases.
- Merging the Q-side and P-side accumulation of the kernels: they need
  different memory layouts to be fast.

## Done on 3 November 2026

The `gt_radmix()` items of tiers 1 to 3. Fits are slightly faster (a few
percent), with the same numbers of iterations, and results equal to those
before up to rounding (the Q-side derivatives changed, see tier 3).

- Tier 1, radmix kernels: no copy of `1 - p` is made any more
  (`one_minus()` is gone; the EM step computes `1 - p` in place, and the
  Newton updates no longer need it).
- Tier 1, radmix fit: the quasi-Newton point is written into the buffers of
  x0 (`qq`/`pq` are gone), keeping the restore of a fixed block after the
  quasi-Newton step (5cf546f). A fit now holds 3 + 2 `n_qn` copies of the
  parameters (9 by default, instead of 10).
- Tier 2, radmix kernels: the segmented and unsegmented Q paths are
  unified in `q_side()` (with `p_side()` as its P-side counterpart, both
  taking a per-block callback); the per-column Newton step is shared in
  `newton_column()`; one driver, `sum_slices()`, serves `loglik()` and
  `deviance()`. The R wrappers of the building blocks are now all in
  `src/radmix_exports.cpp` (`src/radmix_qp_exports.cpp` and
  `src/radmix_geno_debug.cpp` were removed).
- Tier 2, test-only C++: the external-pointer API of the quasi-Newton state
  was replaced by a stateless `radmix_qn_step_cpp()`, and the two projection
  exports were merged into `radmix_project_cpp(x, simplex)`.
- Tier 3, the mu term of the Q tableau: not dropped, but replaced, together
  with the derivatives of the Q block. These are now the plain derivatives
  of the log-likelihood, with the same per-genotype coefficients as the P
  block (OpenADMIXTURE.jl, from which the code started, differentiated a
  homogeneous form that needed `1 - p` and twice the work for the
  Hessians), and the shift is `mu = trace(H) / K^2`, which keeps the
  tableau as well conditioned as H. Results change at the level of
  rounding errors.
- Tier 3, replacing RcppArmadillo with R's LAPACK: no longer worthwhile, as
  the alternating least squares of the SVD start (`src/radmix_als.cpp`)
  also use Armadillo.

## Done since the review

These items were part of the plan, and were done with the bug fixes:

- QP solver: early exit when a vertex solution starts cycling (0850d80).
- `gt_fastmixture()`: arguments checked before python is started
  (ce0bc07).
- `check_count()`/`check_flag()` helpers for the arguments of
  `gt_radmix()` and `gt_fastmixture()` (59eba8d, ce0bc07).
- `gt_single_group()`, shared by `gt_radmix()` and
  `gt_fastmixture()` to get the group of each individual (b3af5e6).

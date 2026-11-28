# Simplification plan for tidygenclust

Simplifications that would make the code shorter or clearer without losing
performance. The plan comes from the code review of October 2026
(`data-raw/code_review.md`) and was updated after the bug fixes that
followed it. Its `gt_radmix()` items were done on 3 November 2026, together
with a rewrite of the Q-side derivatives (see "Done on 3 November 2026"
below), and the other items on 11 November 2026, on the `simplify` branch
(see "Done on 11 November 2026" below). The plan has no open items left.

"Identical" means bit-identical outputs (`gt_radmix()`) or identical results
on the 60 scenarios of `data-raw/clumppling/validate_vs_python.R`
(clumppling). Line counts are approximate.

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

## Done on 11 November 2026

The clumppling, plot, argument-check, test, installer and `R/utils.R` items
of tiers 1 and 2 (about 260 lines fewer in `R/`). Each item was checked to
give the same results as before:

- Tier 1, clumppling alignment: `clumppling_reorder_q_across_k()` is a
  single loop over K (the anchor mode of the smallest K keeps its order,
  the anchor of each larger K follows the mode of the previous K in its
  anchor pair, and the other modes of a K follow its anchor), and
  `clumppling_align_q_wrt_p()` finds the closest column of Q for each
  column of P with `order()`/`duplicated()`, merging columns with
  `clumppling_move_columns()` (which takes the number of columns of the
  result). Checked with `identical()` against the previous code on 10000
  random alignments (with tied distances), on 300 random scenarios like
  those of `data-raw/clumppling/validate_vs_python.R` (also with values of
  K that are not consecutive, and K up to 8) and on the capeverde and
  chicken examples with each option, with the same run times. (The python
  implementation itself was not available to run
  `validate_vs_python.R` again, but the previous code gave the results
  recorded in `validate_vs_python.csv`.)
- Tier 2, `clumppling_run()`: no longer returns `alignment_acrossK`,
  `best_acrossK`, `major_acrossK`, `cost_withinK` and `communities`.
- Tier 2, plots: `autoplot()` checks `k` and warns about subsets of
  individuals once, and the 'major_modes', 'all_modes' and
  'modes_within_k' types share one code path (`plot_mode_list()`), only
  selecting different modes and labels. `plot_modes()` computes all its
  positions with `mode_cells()` (replacing `get_top_xy()`,
  `get_bottom_xy()`, `get_inset_corners()` and the rescaling in
  `make_connection_df()`), its edges with `mode_edges()`, adds every inset
  with `add_inset()`, and builds the structure plots directly instead of
  calling `autoplot(type = "all_modes")`. `mode_k()` and `mode_rank()`
  parse mode labels for `tidy()`, `subset_gt_clumppling()` and the plots.
  The four types of plot, with and without groups, on the capeverde
  example, its subsets by K and by individuals, a single K and the chicken
  example (40 plots) render to byte-identical png files.
- Tier 2, argument checks: `gt_clumppling()` uses `check_flag()`, so its
  message is now "should be TRUE or FALSE"; `check_count()`,
  `check_flag()` and `gt_single_group()` moved to `R/utils.R`.
- Tier 2, `R/utils.R`: `factor(x, levels = unique(x))` replaces the copy of
  the internal `fct_inorder_base()` of tidypopgen.
- Tier 2, installer: `tgc_env_yaml()` resolves the yaml file of one
  environment, for `tgc_resolve_conda_yaml()` and
  `tgc_admixture_install()` (the results of `tgc_resolve_conda_yaml()` are
  the same as before for all the forms of `conda_yaml`, on linux and osx);
  fastmixture is always left out of the yaml.
- Tier 2, tests: `radmix_pseudohaploid_gt()`, `radmix_ref_start()` (with
  an optional seed and range of p, giving the same draws as before) and
  `radmix_match_clusters()` are shared by the tests of `gt_radmix()`.

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

# Simplification plan for tidygenclust

Simplifications that would make the code shorter or clearer without losing
performance. The plan comes from the code review of October 2026
(`data-raw/code_review.md`) and was updated after the bug fixes that
followed it. It has not been implemented yet.

"Identical" means bit-identical outputs (`gt_radmix()`) or identical results
on the 60 scenarios of `data-raw/clumppling/validate_vs_python.R`
(clumppling). Line counts are approximate.

## Tier 1: prototyped, identical results and no loss of speed

- **radmix kernels** (`src/radmix_kernels.cpp`): compute
  `1 - p` in place instead of the full copy made by `one_minus()` in the EM
  and update kernels. Identical, same speed, one K x J matrix less per call;
  about 15 lines.
- **radmix fit** (`src/radmix_fit.cpp`): write the
  quasi-Newton point into the buffers of x0 and drop `qq`/`pq`. Identical,
  and one copy of K(I + J) doubles less. The prototype predates the fix
  that copies a fixed block back into `qq`/`pq` after a quasi-Newton step
  (5cf546f), which it would need to keep.
- **clumppling alignment**: `clumppling_reorder_q_across_k()` as a single
  loop over K (about 20 lines); `clumppling_align_q_wrt_p()` written with
  `order()`/`duplicated()` without merging, and with
  `clumppling_move_columns()` with merging (about 25 lines).

## Tier 2: refactoring off the hot path

- **`gt_fastmixture()`**: the two argparse namespaces differ only in the
  cross-validation arguments, and the four branches that collect the results
  only in which elements they keep. One call and one result format would
  remove about 60 lines.
- **Plots**: `plot_major_modes()`, `plot_all_modes()` and
  `plot_modes_within_k()` (`R/autoplot_gt_clumppling.R`) differ only in the
  modes they select (about 50 lines); `R/plot_modes.R` repeats its inset and
  rescaling code (about 70 lines).
- **`inst/python/py_rfastmixture.py`**: about 120 lines of commented-out
  code (the command line parser and the code that writes files).
- **Test-only C++**: the external-pointer API of the quasi-Newton state
  (`radmix_qn_new_cpp()`, `_push_cpp()`, `_step_cpp()`) can be one
  stateless function (about 75 lines, including generated code). The two
  projection exports, `radmix_project_q_cpp()` (also used for the
  starting values of Q) and `radmix_project_p_cpp()` (only used by
  the tests), can be merged.
- **radmix kernels**: let the segmented Q path handle the
  unsegmented case (about 20 lines), share the per-column Newton step
  between the update kernels (about 10 lines), and use one driver for
  `loglik()` and `deviance()` (about 25 lines).
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

## Tier 3: needs more work

- Dropping the mu term from the Q tableau of the QP solver is closer to
  OpenADMIXTURE.jl, but changes results at the level of rounding errors.
- Replacing RcppArmadillo (used for one small solve in
  `src/radmix_quasi_newton.cpp`) with R's LAPACK would remove a
  dependency; it needs checking against the quasi-Newton tests.

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

## Done since the review

These items were part of the plan, and were done with the bug fixes:

- QP solver: early exit when a vertex solution starts cycling (0850d80).
- `gt_fastmixture()`: arguments checked before python is started
  (ce0bc07).
- `check_count()`/`check_flag()` helpers for the arguments of
  `gt_radmix()` and `gt_fastmixture()` (59eba8d, ce0bc07).
- `gt_single_group()`, shared by `gt_radmix()` and
  `gt_fastmixture()` to get the group of each individual (b3af5e6).

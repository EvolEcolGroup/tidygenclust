# Code review of tidygenclust (October 2026)

A review of the whole package on the `rclump` branch, after the native
clumppling port, the optional conda tools and the simplification of the
OpenADMIXTURE genotype view. It covers correctness, the test suite and ways
to simplify the code without losing performance. Every bug listed under
"High" and "Medium" was reproduced; "Low" items were found by reading the
code and checked where noted.

The review also added inline comments throughout (comment-only changes,
checked mechanically: the code of every edited file is token-for-token
identical), raising the share of comment lines most in the files that had
few: `R/plot_modes.R` (12% to 51% of code lines),
`src/clumppling_louvain.cpp` (17% to 54%), `src/clumppling_assignment.cpp`
(24% to 50%), `R/clumppling_across_k.R` (6% to 35%), `R/clumppling_run.R`
(4% to 30%) and `R/autoplot_gt_clumppling.R` (7% to 30%). Comments that no
longer matched the code, and errors in the roxygen documentation, were
corrected.

## Bugs

### High

1. **Clusters change colour across K when K >= 10** (`plot_q()` in
   `R/autoplot_gt_clumppling.R`, `tidy_q()` in `R/tidy_clumppling.R`). The
   cluster column `q` is a character vector, so its levels sort as "1", "10",
   "11", ..., "2", and `scale_fill_distruct()` assigns colours (and the
   stacking order) by level. The same aligned cluster is then drawn in a
   different colour in each K: in the chicken example of the vignette, the
   dominant cluster of the first individual (cluster 5 in K17M1, K19M1 and
   K21M1) is drawn in three different colours.
2. **`gt_fastmixture()` ignores the selection of individuals and loci of a
   gen_tibble** (`R/gt_fastmixture.R`, `inst/python/py_rfastmixture.py`).
   The backing file is read in full and assumed to have as many
   individuals and loci as the gen_tibble, but filtering, reordering and `select_loci()` only
   change the indices of a gen_tibble, not its backing file. A filtered
   gen_tibble fails with "AssertionError: Backing file doesn't match number
   of samples and SNPs!"; a reordered one runs silently, but the rows of Q
   are in the order of the backing file while `id` follows the gen_tibble,
   so each individual gets someone else's ancestry (P and supervised labels
   have the same problem).
3. **`gt_fastmixture(threads = )` has no effect, and changes R's
   environment** (`inst/python/py_rfastmixture.py`). The wrapper sets
   `OMP_NUM_THREADS` (and similar variables) after R has already started
   the OpenMP runtime, so it is ignored: `threads = 1` used about three
   cores. As python shares the process environment with R, the variables
   also stay set in the R session after the call.

### Medium

4. **`gt_clumppling()` hangs with a non-finite `cd_res`** (`R/gt_clumppling.R`,
   `src/clumppling_louvain.cpp`). With `cd_res = NA`, `NaN` or `Inf` the
   modularity is NaN, the stopping test of louvain is never true, and the
   C++ loop (which does not check for interrupts) runs forever. Other
   options are not validated either (`merge = "no"` gives "invalid 'x' type
   in 'x && y'"), and `clumppling_louvain_cpp()` accepts infinite weights.
5. **`subset_gt_clumppling()`** (`R/gt_clumppling.R`):
   - with `k` covering all of `K_range`, it fails with "invalid argument to
     unary operator" (`x$cost_acrossK[-NULL]`);
   - modes and costs are selected with `grep(paste0("K", k))`, so `k = 2`
     also matches K20 to K29;
   - `k` must be consecutive integers, so values of K with gaps (e.g. 17,
     19, 21 in the chicken example) cannot be selected, while values of `k`
     that are not in `K_range` silently give an empty object;
   - `indivs` is not checked against `N` (`indivs = 390:420` with N = 399
     returns rows of NA).
6. **`gt_openadmixture()`** (`R/gt_openadmixture.R`):
   - it refuses `k >= n_indiv` even with `fix_p = TRUE`, where each
     individual is estimated independently, so a small set of samples
     (e.g. one new sample) cannot be projected onto a reference panel;
   - `init_p` is only checked for its dimensions: a P with the loci in a
     different order (e.g. from another dataset) is accepted silently and
     gives wrong results;
   - without `seed`, turning on `crossval` changes the estimates of later
     runs and values of K (cross-validation draws from the global random
     number generator), although the documentation says it does not change
     the estimates.
7. **Installation of the optional tools** (`R/tgc_tools_install.R`):
   - `tgc_tools_install(tools = "admixture", reset = TRUE)` removes
     `ctidygenclust` (and the fastmixture installed in it), and on macOS
     never recreates it;
   - fastmixture is installed with `pip --force-reinstall`, which also
     reinstalls its dependencies and overwrites the numpy pinned by the
     conda snapshot;
   - on macOS, `cadmixture86` is updated by calling `conda` from the PATH
     rather than the conda binary used by reticulate (often missing in
     RStudio), with an unquoted path;
   - a single `conda_yaml` file is used for both environments on macOS
     (only "auto" and NULL should be recycled).
8. **Plots of gt_clumppling objects** (`R/autoplot_gt_clumppling.R`,
   `R/plot_modes.R`):
   - on objects subset with `subset_gt_clumppling(indivs = )`, the group
     lines and labels are drawn at the wrong positions (bars are drawn at
     the original indices, groups from 0), and the label of the number of
     runs in `type = "modes"` stretches the axis;
   - a factor `group` fails (`rle()` rejects factors);
   - `type = "modes_within_k"` does not check `k`: a value that is not in
     `K_range` gives an empty plot that fails when printed;
   - with a single K, `type = "modes"` warns ("no non-missing arguments to
     max");
   - the vignette and the documentation said that edges of the 'modes' plot
     are labelled with their cost, but they are only coloured.
9. **The QP solver of Q cycles until its iteration limit at vertex
   solutions** (`src/openadmixture_qp_sweep.cpp`). When the solution of a
   column is at a vertex, a parameter is freed, takes a step of exactly 0,
   and is fixed again, up to the limit of 1000 passes; about 4% of columns
   do this near convergence, making `update_q` 6-21% slower. An early exit
   gives identical results.

### Low

- **Argument validation in `gt_openadmixture()`**: `fix_q = NA` is treated
  as TRUE by C++; `rtol = NA` (or negative) silently disables convergence;
  `max_iter`, `n_qn` and `n_cores` give a generic C++ error; non-integer `k`,
  `n_runs` and `n_folds` are truncated; seeds are only checked by
  `set.seed()`, after the first fit.
- **A gen_tibble grouped by more than one variable** makes
  `gt_openadmixture()` fail only after all the fits have run.
- **Argument validation in `gt_fastmixture()`**: the error message names
  `seeds` and `repeats` rather than `seed` and `n_runs`; the default `seed`
  fails with `n_runs > 1`; arguments are only checked after python is
  started; `cv = NA` gives a cryptic error. The `safety` argument is never
  used by the python wrapper (fastmixture 1.3.0 has no such option).
- **`message()` is looked up in the global environment** by the C++ fit
  (`Function message("message")` in `src/openadmixture_fit.cpp`), so a
  user-defined `message()` receives the log of verbose fits.
- **`fix_q = TRUE` does not keep Q exactly fixed**: the quasi-Newton step
  re-projects the fixed block, which changes it at the level of rounding
  errors (P stays exact).
- **Internal C++ functions crash R on degenerate inputs** (not reachable
  from the exported functions): `openadmixture_q_segments_cpp()` and the Q
  kernels divide by zero with no individuals (SIGFPE); `project_q()` with
  K = 0; `openadmixture_qn_step_cpp()` writes past the end of its output if
  q and p have different numbers of rows; `clumppling_align_ilp_cpp()`
  with an empty P; the cross-validation seed accepts `Inf`/very large
  values (undefined behaviour in the cast to uint64_t).
- **Cosmetic errors**: missing space in the Windows error message
  ("subsystemfor"); unbalanced parenthesis, and 'groups' instead of 'group',
  in the error about non-adjacent groups; `ci_install = TRUE` from outside
  the package root gives "replacement has length zero"; the homebrew path
  is hard-coded for Apple silicon.
- **DESCRIPTION**: `plot_layout(axes = )` needs patchwork >= 1.2.0 and
  `linewidth` needs ggplot2 >= 3.4.0, which are not declared;
  `R_unif_index()` (used by louvain) is only part of R's C API since R
  3.6.0, while the package declares R >= 3.5.0.

## Test suite

- The plotting tests only check the class of the objects: nothing is
  rendered, so errors at printing time pass; some `autoplot()` calls have no
  expectation. Missing cases: factor `group`, invalid `k`, a single K, one
  mode per K, K >= 10, subset objects with `group`.
- `test_fastmixture.R`: `expect_true(all(names(x) == names(x)))` is always
  true, and `inherits(x, c("list", "gt_admix"))` is true for any list. No
  test uses a filtered or reordered gen_tibble, checks `id`/`group`, or the
  argument checks.
- None of the bugs above is covered by a test.
- Fragile: the chicken test relies on R's seed 7 reproducing python's modes
  for K = 17 and 19 (an implementation detail of the random number
  generator); the yaml filtering tests depend on the number of lines of the
  latest conda snapshot.
- Weak or dead: `n_cores = 1` vs `2` compared with `expect_equal` although
  results are bit-identical; the fixed-Q test uses `tolerance = 1e-4`;
  `openadmixture_has_real_qp()` skips are dead (the solver is always
  compiled in); a kernel test exercises its `p = 0` branch only if the
  simulation happens to produce such loci.
- The installer's decisions (reset, existing environments, operating
  systems, yaml recycling) are untested; they could be tested by mocking
  reticulate's conda functions with `testthat::local_mocked_bindings()`.

## Simplification plan

Simplifications that do not lose performance, by tier. Tier 1 was
prototyped; "identical" means bit-identical outputs (OpenADMIXTURE) or
identical results on the 60 scenarios of
`data-raw/clumppling/validate_vs_python.R` (clumppling).

### Tier 1: prototyped, identical results and no loss of speed

- OpenADMIXTURE kernels: compute `1 - p` in place instead of the full copy
  made by `one_minus()` (identical, same speed, one K x J matrix less per
  call; about 15 lines).
- OpenADMIXTURE fit: write the quasi-Newton point into the buffers of x0
  and drop `qq`/`pq` (identical, one copy of K(I + J) doubles less).
- QP solver: early exit when a vertex solution starts cycling (bug 9).
- clumppling: `clumppling_reorder_q_across_k()` as a single loop over K
  (about 20 lines); `clumppling_align_q_wrt_p()` without merging written
  with `order()`/`duplicated()`, and with merging using
  `clumppling_move_columns()` (about 25 lines).

### Tier 2: refactoring off the hot path

- `gt_fastmixture()`: the four branches that collect the results differ
  only in which elements they keep; one call and one result format would
  remove about 60 lines, and the arguments could be checked before python
  is started.
- Plots: `plot_major_modes()`, `plot_all_modes()` and
  `plot_modes_within_k()` differ only in the modes they select (about 50
  lines); `R/plot_modes.R` repeats its inset and rescaling code (about 70
  lines).
- `inst/python/py_rfastmixture.py`: about 120 lines of commented-out code.
- Test-only C++: the external-pointer API of the quasi-Newton state
  (`openadmixture_qn_new_cpp()`, `_push_cpp()`, `_step_cpp()`) can be one
  stateless function (about 75 lines including generated code); the two
  projection exports can be merged.
- OpenADMIXTURE kernels: let the segmented Q path handle the unsegmented
  case (about 20 lines), share the per-column Newton step between the update
  kernels (about 10 lines), and use one driver for `loglik()` and
  `deviance()` (about 25 lines).
- Argument checks: small `check_count()`/`check_flag()` helpers (which also
  add the missing checks).
- Tests: shared helpers for pseudohaploid gen_tibbles, starting values and
  cluster matching (about 50 lines).
- Installer: resolve the yaml file per environment with one helper, and
  always leave fastmixture out of the yaml (it is reinstalled from github).
- `R/utils.R`: `factor(x, levels = unique(x))` instead of copying an
  internal function of tidypopgen.
- `clumppling_run()` returns several elements that nothing uses (about 10
  lines).

### Tier 3: needs more work

- Dropping the mu term from the Q tableau of the QP solver is closer to
  OpenADMIXTURE.jl but changes results at the level of rounding errors.
- Replacing RcppArmadillo (used for one 3 x 3 solve) with R's LAPACK would
  remove a dependency; it needs checking against the quasi-Newton tests.

### Not recommended

- Reading genotypes through bigstatsr's `SubBMCode256Acc`: about 19%
  slower than the custom genotype view.
- Simpler sums in clumppling (`colSums()` for the representatives, the
  within-K cost from the assignment objective): they change which of
  several near-identical replicates is the representative, breaking the
  agreement with python.
- `igraph::cluster_louvain()` instead of the louvain port: a different
  algorithm, so modes would differ from python in ambiguous cases.
- Merging the Q-side and P-side accumulation of the kernels: they need
  different memory layouts to be fast.

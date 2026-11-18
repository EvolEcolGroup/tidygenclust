# tidygenclust: summary

Updated 2026-11-18 (manica01, overnight). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py` of the dev-notes skill).

## tidygenclust: dramanica/claude_jobs#3

tidygenclust: `gt_radmix()`, the native implementation of ADMIXTURE's method (port, speed-ups, consolidation), and its extension to genotype likelihoods. Earlier log of the likelihoods: dramanica/claude_jobs#12.

**Status.** PR #38 (`radmix_extensions`, squash-merged into `dev` as df747aa on 18 Nov, CI green) added a fused pass for supervised fits (fits 0.88-0.89 of the time, the same estimates), cross-validation and SVD starts for genotype likelihoods, and inlined likelihood helpers; with tidypopgl's faster source (its branch `radmix_source_speed`), fits of likelihoods take 0.80-0.83 of the time. tidypopgl's `radmix_source_speed` waits for review. Handover: `dev_notes/handover.md`.

`dev` holds PR #36 (merged 17 Nov: the speed-ups of `newspeed`, supervised
analyses, the benchmark article with the full HGDP run) and genotype
likelihoods (#35, 14 Nov), used by tidypopgl's `gl_workflow`, which
dropped NGSadmix. PR #38 added:
- **Supervised fits** evaluate each point and start the next iteration in
  one pass, as unsupervised fits do: the free individuals' Q side gives
  their update and log-likelihood, and a log-likelihood-only pass the
  references'. HGDP K = 5 (139 or 660 references): 0.88-0.89 of the time on
  1 and 8 threads, Q and P bit-identical.
- **Cross-validation of likelihoods:** the deviance
  `log max_x f(x) - log f(qp)` of each held-out individual and locus
  (ADMIXTURE's deviance for hard calls); it chose the simulated K in 12 of
  12 cases (depths 1-4 and genotypes).
- **SVD start of likelihoods:** a randomised SVD of the posterior dosages
  under emMAF frequencies, computed by passes over the source; at K = 5 on
  simulated data it converged 2.6-4x faster than random starts for two
  seeds of three, to higher maxima.
- **Speed of likelihood fits:** `slice_gl()` and `loglik_gl_slice()` now
  run in the AVX2 clones (they were 17% of a fit in baseline code), and
  tidypopgl's source decodes with one table lookup and AVX2 gathers (it
  was 23% of a fit).

The designs are in `admixture_implementation.md`, the timings in
`radmix_performance.md` and `radmix_extensions_2026-11/`.

**Next steps**
- [ ] **P2** [review] tidypopgl's `radmix_source_speed`: merge it; then update tidypopgl's vignette (it says cross-validation is not available for likelihoods).
- [ ] **P2** [claude] then [review] `articles_revised`: resolve its conflict in `pseudohaploids.Rmd`, then review and merge.
- [ ] **P3** [decide] Keep the 1e-5 bounds for likelihoods, or use NGSadmix's 1e-9? With 1e-5, the log-likelihood is about 1e-6 below NGSadmix's. Recommendation: keep 1e-5.
- [ ] **P3** [cleanup] On laptopX1, delete the untracked `data-raw/radmix/perf_verification/`; on manica01, delete the merged `radmix_extensions` (local and origin); delete `~/tgc_newspeed_work`, `~/tgc_bench_article`, `~/tgc_supfused_work` and `vignettes/articles/data/hgdp/runs` (about 7 GB of builds and raw outputs).

Backlog (P3, in `handover.md`): `-fprofile-update=single` under covr, the vignette's P comparison as code, P as a matrix or a data.frame (decide), bootstrap standard errors, the remaining speed leads, an EM-based PCA for the SVD start of likelihoods.

**Done**
- Extensions of `gt_radmix()`: fused pass for supervised fits, cross-validation and SVD starts for likelihoods, `svd_power`: PR #38 squash-merged into `dev` (df747aa), Andrea, 18 Nov 2026.
- laptopX1: 17 local branches without origin and the `coverage-speed` worktree deleted (tips in `/media/andrea/annex/tmp/claude/tidygenclust/deleted_branches_2026-11-18.txt`), 18 Nov 2026.
- Code derived from other software labelled in `gt_radmix()`, `gt_fastmixture()` and `gt_clumppling()`, licences completed in `COPYRIGHTS` and `cph` roles (ADAMIXTURE, cp-algorithms and ggpubr were uncredited); text for the paper in `dev_notes/paper_methods_text.md`: merged into `dev` (4c8464d), Andrea, 18 Nov 2026.

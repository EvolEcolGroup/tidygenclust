# tidygenclust: summary

Updated 2026-11-27 (laptopX1, on the local branch `ada196_article`). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py` of the dev-notes skill).

## tidygenclust: dramanica/claude_jobs#3

tidygenclust: `gt_radmix()`, the native implementation of ADMIXTURE's method (port, speed-ups, consolidation), and its extension to genotype likelihoods. Earlier log of the likelihoods: dramanica/claude_jobs#12.

**Status.** `gt_radmix()` streams data larger than the memory: merged into the local branch `ada196_article` (version 0.1.5.9007, not pushed, by Andrea's choice), measured on manica01, and shown in the benchmark article's new section: 50 GB of genotypes under an 8 GiB cap took 34% longer than in memory on an NVMe SSD, with identical results. The Q side now takes chunks in memory too (Andrea, 26 Nov), 15-21% faster with 8,192 individuals or more; the article's fits in memory were rerun with it, and the section waits for Andrea's review. Handover: `dev_notes/handover.md`.

PR #38 (`radmix_extensions`, squash-merged into `dev` as df747aa on 18 Nov, CI green) added a fused pass for supervised fits (fits 0.88-0.89 of the time, the same estimates), cross-validation and SVD starts for genotype likelihoods, and inlined likelihood helpers; with tidypopgl's faster source (its branch `radmix_source_speed`), fits of likelihoods take 0.80-0.83 of the time. `dev` is now 0.1.5.9006, which tidypopgl's `gl_workflow` requires since its PR #47.

Streaming (`radmix_streaming.md`, "Desktop"): with 10,000 individuals and a third of a 1.5 GB file cached, a short fit takes 9.1 s on the NVMe (6.6 s in memory; `dev` 105.5 s, reading the file 100 times), 25 s at 500 MB/s (`dev` 318 s). The prefetch only read 1/16 of each chunk on the desktop; fixed (madvise in pieces of 128 KB), it saves 20-27% with many individuals.

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
- [ ] **P1** [review] the section "Data larger than the memory" of the benchmark article (`ada196_article`, `vignettes/articles/admixture_benchmark.Rmd`; in-memory fits rerun with chunks on 27 Nov): then the PR into `dev` when Andrea says so.
- [ ] **P3** [cleanup] delete the article's data and the builds on manica01 (`~/project_temp/claude/tidygenclust/stream/`, about 75 GB) once the article is final.
- [ ] **P2** [claude] then [review] `articles_revised`: resolve its conflict in `pseudohaploids.Rmd`, then review and merge.

Backlog (P3, in `handover.md`): `-fprofile-update=single` under covr, the vignette's P comparison as code, P as a matrix or a data.frame (decide), bootstrap standard errors, the remaining speed leads, an EM-based PCA for the SVD start of likelihoods.

**Done**
- Streaming defaults kept (stream above 0.8 times the memory available, prefetch in pieces of 128 KB, alternation on), and the Q side in chunks for data in memory too (c521e36, 79e21c4). Andrea, 26 Nov 2026.
- Bounds of likelihood estimates: keep `[1e-5, 1 - 1e-5]`, as for genotypes, not NGSadmix's 1e-9 (`admixture_implementation.md`). Andrea, 20 Nov 2026.
- Clean-up: local `radmix_extensions` and untracked `data-raw/radmix/perf_verification/` deleted on laptopX1; `~/tgc_newspeed_work`, `~/tgc_bench_article`, `~/tgc_supfused_work` and `vignettes/articles/data/hgdp/runs` (about 6.4 GB) deleted on manica01. 20 Nov 2026.
- tidypopgl's PR #47 (`radmix_source_speed`: faster source, requires 0.1.5.9006, cross-validation in its vignette) squash-merged into `gl_workflow` (201f29f). Andrea, 19 Nov 2026.
- Version 0.1.5.9006 on `dev` (d889f63), the first with #38's cross-validation and SVD starts for likelihoods; required by tidypopgl's PR #47. Andrea, 19 Nov 2026.
- `radmix_extensions` (merged as #38) deleted on origin and manica01. Andrea, 19 Nov 2026.

# tidygenclust: summary

Updated 2026-11-16 (laptopX1, evening). One section: the package's issue in
dramanica/claude_jobs, with the single list of its tasks. Its first
paragraph, status and next steps are copied into the issue
(`scripts/sync_issue.py` of the dev-notes skill).

## tidygenclust: dramanica/claude_jobs#3

tidygenclust: `gt_radmix()`, the native implementation of ADMIXTURE's method (port, speed-ups, consolidation), and its extension to genotype likelihoods. Earlier log of the likelihoods: dramanica/claude_jobs#12.

**Status.** PR #36 (branch `newspeed_port`: the speed-ups of `newspeed`, the supervised analyses, and #37's test changes) passes R CMD check on the desktop, and its fits take 0.78-0.85 of `dev`'s time (0.92 at HGDP K = 10 on 1 thread; desktop, 16 Nov); the benchmark article now shows the full HGDP run (desktop, 16 Nov), so #36 waits only for review. `gt_radmix()` on genotype likelihoods is in `dev` (#35, 14 Nov) and used by tidypopgl's `gl_workflow`, which dropped NGSadmix. Handover: `dev_notes/handover.md`.

The port rebuilt each package commit of the two branches on top of #35,
which had made the kind of data (genotypes or genotype likelihoods) a
template parameter of the kernels, and applied the resolutions of the merge
`newspeed_supervised`. The fused pass and the supervised index also work on
genotype likelihoods (the NGSadmix fixture test passes); supervised
analyses of genotype likelihoods are new in the port, tested against the
same fits on hard-call genotypes. Decisions applied: refits of
cross-validation to `cv_rtol = 1e-7`, `MIN_Q_TASKS` = 128, no subset start;
FMA rounding, the inlined log-likelihood and the automatic bit mask kept;
the supervised API as built. The kits, reports and desktop results stay on
the old branches, archived under the tags `archive/newspeed` and
`archive/supervised`; their conclusions are in `radmix_performance.md`.

The design is in `admixture_implementation.md` ("Genotype likelihoods"),
the validation and benchmarks in `data-raw/radmix/gl/REPORT.md`: the same
log-likelihood as NGSadmix to 1e-9 on its fixture, 35-45 times faster on
simulated data. PR #36 (`newspeed_port`) extends the fused
pass and supervised analyses to genotype likelihoods.

**Next steps**
- [ ] **P1** [review] Merge PR #36 (checks, timing and the benchmark article done) and delete `newspeed`, `supervised`, `newspeed_supervised` and `coverage-speed`: #36 also adds supervised analyses of genotype likelihoods, which were not asked for.
- [ ] **P2** [claude] then [review] `articles_revised`: resolve its conflict in `pseudohaploids.Rmd`, then review and merge.
- [ ] **P3** [claude] Extensions of `gt_radmix()`: a fused pass for supervised fits (update every column of Q, then restore the references: 4 passes per iteration instead of 5); for likelihoods, cross-validation (deviance) and an SVD start (both stop with an error for now); optional, a faster source in tidypopgl (adjacent individuals read without the index).
- [ ] **P3** [decide] Keep the 1e-5 bounds for likelihoods, or use NGSadmix's 1e-9? With 1e-5, the log-likelihood is about 1e-6 below NGSadmix's.
- [ ] **P3** [cleanup] Delete origin `gl` and `simplify`, the local `gl_gl16`, `gl_private` and stale branches, and the untracked copy of `data-raw/radmix/perf_verification/`; on manica01, disable `tgc-desktop-checks.timer` and prune `~/tgc_newspeed_work`.

Backlog (P3, in `handover.md`): `-fprofile-update=single` under covr, the vignette's P comparison as code, P as a matrix or a data.frame (decide), bootstrap standard errors, the remaining speed leads.

**Done**
- Full HGDP benchmark of `admixture_benchmark.Rmd` with `newspeed_port`; article switched to it: `5eb3c73`, 16 Nov 2026.
- Timing of `newspeed_port` against `dev`: 0.78-0.85 of `dev`'s time (0.92 at HGDP K = 10, 1 thread), same iterations and log-likelihoods; `aa5795e`, `radmix_performance.md`, 16 Nov 2026.
- PR #37 closed, superseded by #36: Andrea, 16 Nov 2026.
- Desktop verification of `newspeed_port` (3f2a8e4): CRAN mode 0 ERROR/WARNING, all tests pass; PR #36 comment, 16 Nov 2026.
- Decisions on the newspeed defaults and the supervised API: Andrea, 15 Nov 2026.
- Port of `newspeed_supervised` onto `dev`: PR #36, 15 Nov 2026.
- Why covr was slow: gcov counters in one multithreaded test file, now skipped under covr; covr's flags unchanged (#37's `-O1` dropped); CRAN tests on 2 threads at most. Andrea, 16 Nov 2026; `covr_timing_2026-11.md`.
- Genotype likelihoods in `gt_radmix()` (the NGSadmix model with radmix's algorithm): in `dev`, #35, 14 Nov 2026; claude_jobs#12 merged into this issue, 16 Nov 2026.

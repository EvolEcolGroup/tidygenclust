# tidygenclust: summary

Updated 2026-11-15 (laptopX1). One section per job in dramanica/claude_jobs.
The first paragraph of each section and its next steps are copied into the
issue (`scripts/sync_issue.py`).

## tidygenclust radmix: dramanica/claude_jobs#3

Native implementation of ADMIXTURE's method in tidygenclust (`gt_radmix()`): port, speed-ups and consolidation.

**Status.** The speed-ups of `newspeed` and the supervised analyses of `supervised` are ported onto `dev` in PR #36 (branch `newspeed_port`, open), with Andrea's decisions of 15 Nov applied, and all CI checks pass, coverage included (an out-of-bounds index found by the coverage job's bounds checks was fixed, 3192e45). Fits take 0.75-0.85 of the time of the previous code (desktop, measured before the port). Handover: `dev_notes/handover.md`.

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

**Next steps**
- [ ] **P1** [review] Review and merge PR #36 (`newspeed_port` into `dev`), then delete `newspeed`, `supervised` and `newspeed_supervised` (archived under tags): it also adds supervised analyses of genotype likelihoods, which were not asked for.
- [ ] **P2** [review] Review and merge PR #37 (`coverage-speed`): covr at -O1 with libstdc++ bounds checks, and a thread test of `gt_adamixture()` that no longer fails when `OMP_NUM_THREADS` is set. (The `dev` coverage failures were Codecov's deactivated repository, reactivated by Andrea on 15 Nov.)
- [ ] **P2** [claude] then [review] `articles_revised`: resolve its conflict in `pseudohaploids.Rmd`, then review and merge.
- [ ] **P3** [desktop] Time the port against `dev` (a20eb08) and `newspeed` (13f89be): the `iter8` step of the newspeed kit, HGDP K = 5 and 10, 1 and 8 threads, about 30 min; #35's templating may shift the gains.
- [ ] **P3** [claude] A fused pass for supervised fits (update every column of Q, then restore the references): 4 passes per iteration instead of 5.
- [ ] **P3** [cleanup] Delete origin `gl` and `simplify`, the local `gl_gl16`, `gl_private` and stale branches, and the untracked copy of `data-raw/radmix/perf_verification/`; on manica01, disable `tgc-desktop-checks.timer` and prune `~/tgc_newspeed_work`.

Backlog (P3, in `handover.md`): where the covr overhead goes, the vignette's P comparison as code, P as a matrix or a data.frame (decide), bootstrap standard errors, the remaining speed leads.

**Done**
- Decisions on the newspeed defaults and the supervised API: Andrea, 15 Nov 2026.
- Port of `newspeed_supervised` onto `dev`: PR #36, 15 Nov 2026.

## tidygenclust gl: dramanica/claude_jobs#12

`gt_radmix()` on genotype likelihoods (the NGSadmix model with radmix's algorithm), with tidypopgl's `gl_tibble`s as a source.

**Status.** Done: in tidygenclust `dev` (#35, 14 Nov), and used by tidypopgl's `gl_workflow`, which dropped NGSadmix; CI passes on both. Only extensions remain; close this issue once they move to issue 3 or to tidygenclust's own issues. Handover: `dev_notes/handover.md`.

The design is in `admixture_implementation.md` ("Genotype likelihoods"),
the validation and benchmarks in `data-raw/radmix/gl/REPORT.md`: the same
log-likelihood as NGSadmix to 1e-9 on its fixture, 35-45 times faster on
simulated data. The port of `newspeed_port` (issue 3) extends the fused
pass and supervised analyses to genotype likelihoods.

**Next steps**
- [ ] **P3** [claude] Cross-validation (deviance) and an SVD start for likelihoods. Both stop with an error for now.
- [ ] **P3** [decide] Keep the 1e-5 bounds for likelihoods, or use NGSadmix's 1e-9? With 1e-5, the log-likelihood is about 1e-6 below NGSadmix's.
- [ ] **P3** [claude] Optional: a faster source in tidypopgl (adjacent individuals read without the index).
- [ ] **P3** [cleanup] Delete origin `gl`, and the local `gl_gl16` and `gl_private` (also listed in issue 3).

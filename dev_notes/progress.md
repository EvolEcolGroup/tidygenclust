# tidygenclust: progress

The log of major activity, newest first. Never rewrite old entries: add a
correction as a new entry. Earlier history (before 15 Nov 2026) is in the
comments of dramanica/claude_jobs#3 (radmix) and #12 (genotype
likelihoods), and in `radmix_performance.md`.

## 2026-11-15 (laptopX1): newspeed and supervised ported onto dev

- **Done:** rebuilt the package commits of `newspeed` (fused pass, CV bit
  mask, AVX2 inlining, cheaper log-likelihood, tile settings,
  `cv_max_iter`/`cv_rtol`) and `supervised` (supervised analyses, ADMIXTURE
  fixture, docs) on top of #35, with the resolutions of the merge
  `newspeed_supervised`; the kernels keep the GL template parameter, and
  `slice_gl()` reads the free individuals of a supervised fit through their
  index, so both features also work on genotype likelihoods. Then the two
  decided defaults, and the notes.
- **Where:** branch `newspeed_port`, commits 32aeab2 to the commit of this entry, on
  `dev` a20eb08; tags `archive/newspeed` (ec96315) and `archive/supervised`
  (fc43f9b) keep the old branches' records.
- **Result:** the full test suite passes on an -O2 build (3918
  expectations, none skipped, with the ADMIXTURE, OpenADMIXTURE.jl and
  NGSadmix fixtures); the new likelihood tests fail without the indexed
  read. Not timed yet (desktop item in `handover.md`).
- **Decisions:** (Andrea) newspeed defaults as recommended: `cv_rtol` 1e-7
  (implemented as `max(rtol, 1e-7)`), no subset start, FMA rounding
  accepted, inlined log-likelihood, `MIN_Q_TASKS` 128, automatic bit mask;
  the supervised API as built. (Claude) supervised analyses allowed for
  genotype likelihoods, for Andrea's review in the PR.
- **Notes changed:** `radmix_performance.md`, `admixture_implementation.md`;
  new `README.md`, `summary.md`, `handover.md`, `progress.md`.

# Developer notes of tidygenclust

Notes for developers and Claude sessions. Start with `summary.md` and
`handover.md`.

## Living documents

- `summary.md`: status and priorities, one section per job in
  dramanica/claude_jobs (for Andrea).
- `handover.md`: branches, next tasks, desktop queue and backlog (for the
  next session).
- `progress.md`: the log of major activity, newest first.

## `gt_radmix()` (native ADMIXTURE)

- `admixture_implementation.md`: the design, from first principles: the
  model, the algorithm, the kernels, cross-validation, supervised analyses
  and genotype likelihoods.
- `radmix_performance.md`: what made it faster, what was tried and
  rejected, the leads left, and where the full records are archived
  (tags `archive/radmix-perf-verify`, `archive/newspeed`,
  `archive/supervised`); the timing of the port against `dev` (16 Nov
  2026) is in `radmix_port_timing_2026-11/`.
- `../data-raw/radmix/README.md`: the test fixtures (ADMIXTURE,
  OpenADMIXTURE.jl, NGSadmix) and how to regenerate them.
- `../data-raw/radmix/gl/REPORT.md`: the validation and benchmarks of
  genotype likelihoods.
- `../data-raw/radmix/Notes_to_do`: the older backlog.

## `gt_clumppling()`

- `clumppling_implementation.md`: the native port of clumppling.

## Package-wide

- `simplification_plan.md`: the simplification plan (completed in #34).
- `covr_timing_2026-11.md`: why covr was slow (gcov counters in
  multithreaded tests) and what was decided; data in `covr_timing/`.

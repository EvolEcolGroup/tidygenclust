# Desktop runs: preparation (manica01, Friday 13 November 2026)

The desktop is manica01 (i9-13900, 24 cores, SMT off: P-cores 0-7, E-cores
8-23). The timed runs of `HANDOVER.md` are scheduled for Sunday 15
November, 09:00. A systemd user timer, `tgc-desktop-checks.timer`, starts
`run_desktop_all.sh`. The script waits until the machine is idle, then runs
the phases in order. When it finishes, it commits the summaries here, pushes
them, and comments on claude_jobs #3. Linger is enabled for the user, so the
timer also fires if the user is logged out.

The work done beforehand ran on the E-cores while the machine was in use,
so none of it was timed:

- **Builds.** All the labels of `README.md` ("Commands") are in
  `~/tgc_newspeed_work/libs`, built with `build_variants.sh`. Every build
  shows `-O2` and the AVX2 clones (5; `base` and `sup` have 4).
  - The newspeed labels: `base`, `fused`, `cvmask`, `inline`, `final`,
    `head`, `final_llnoinl`, `final_q64`, `final_q128`, `final_guided` and
    `final_tj32`.
  - `sup`: `supervised` at 7c576ab.
  - `merged`: `newspeed_supervised`, the two branches merged.
- **Data.** The kit's backing files for `hgdp`, `hgdp_full`, `hgdp_quick`
  and `eur` are in `~/tgc_newspeed_work/data`. The supervised HGDP
  gen_tibble is in `~/tgc_newspeed_work/sup`.
- **Smoke tests.** All of them pass:
  - the kit, with toy data at 1 and 2 threads;
  - `newspeed_algo/run_desktop.sh smoke`;
  - the supervised timing script and the comparison with ADMIXTURE
    `--supervised` (pops, -j4). For the 904 individuals to estimate, Q
    agrees within 4.3e-5 (mean 3.5e-6), and gt_radmix's log-likelihood is
    higher by 1.1.
- **Test suites.** Each suite ran with `NOT_CRAN=true` against the
  installed `-O2` build, and none had a failure:

  | branch | build | tests passed |
  |---|---|---|
  | `supervised` | `sup` | 3707 |
  | `newspeed` | `head` | 3700 |
  | `newspeed_supervised` | `merged` | 3760 |
- **Supervised thread invariance.** On HGDP pruned, K = 5, n_cores 1 vs 8
  give `identical()` Q, P and log-likelihood, for both pops and regions.
  This is item 5 of the supervised desktop tests
  (`data-raw/radmix/supervised/thread_invariance.R`).
- **CV trajectories.** These are step 3, `cv`, of `newspeed_algo` and give
  machine-independent results. They run on the E-cores as the unit
  `tgc-cvtraj`, which stops at 08:45 on Sunday. The Sunday run finishes
  whatever is left, at the end.

The Sunday run uses this order: kit step 1 (kernels, iter8, whole, seeds,
cv), CV-policy wall-clock (`cv_timing`), warm-ups (hgdp, eur), supervised
(timings at 1 and 8 threads, per-pass costs, ADMIXTURE -j8 against
gt_radmix at 8 cores), then tiles/schedule, then the rest of the CV
trajectories and the `newspeed_algo` summary.

- Tiles come after the supervised runs because `HANDOVER.md` makes them
  optional ("only if step 1 leaves time").
- `warmup hgdp_full` is not scheduled: `HANDOVER.md` says to run it only
  if the subset start is still a candidate.
- No phase starts after Monday 16 November, 07:00.

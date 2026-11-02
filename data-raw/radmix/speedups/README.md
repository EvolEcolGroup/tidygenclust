# data-raw/radmix/speedups

The record of the work that sped up `gt_radmix()` (October-November 2026):
the SVD start (`init = "svd"`), the faster kernels, and the default
convergence threshold of 1e-8. The work was done on a laptop (Intel
i7-10510U, 4 cores) and verified on a desktop (Intel i9-13900, 8 P-cores).

| file | what it holds |
|---|---|
| `speedup_plan.md` | the plan and its history: the comparison with ADAMIXTURE and admixer, the perf profile, the experiments in order of priority (including the SVD start, experiment 5), the protocol |
| `kernel_experiments_laptop.md` | the full log of the four kernel experiments on the laptop (divisions, decoding, compiler flags and multiversioning, Hessian accumulation) |
| `handover.md` | the summary of the laptop work, the decisions left to the user, and the plan of the verification on the desktop |
| `verification_desktop.md` | the report of the verification on the desktop, with the laptop's figures next to the desktop's, the multi-core scaling, the SVD start on the full HGDP, and the evidence for the decisions |
| `verification_desktop_scaling.png` | the scaling plot of that report |

Outcome. The decisions in section 2.4 of `handover.md` were taken as
follows; the changes are in the branch `radmix-speedups`:
- the SVD start is an option (`init = "svd"`); random starts stay the
  default. Its ALS runs with one BLAS thread, a fix found during the
  verification (section "New: the BLAS-thread fix" of
  `verification_desktop.md`);
- the kernel speed-ups are kept with multiversioning (an AVX2 + FMA version
  of the kernels on x86-64 Linux with GCC >= 12);
- the default `rtol` is 1e-8 instead of 1e-7.

The scripts, patches and results (CSV) these documents refer to are in
`data-raw/radmix/perf_verification/` on the branch `radmix-perf-verify`
(commit f423896), which is kept as the archive of this work; they are not
copied here. The scripts of the SVD start exploration are in the parent
directory (`explore_svd_start.R`, `explore_svd_start_refined.R`).

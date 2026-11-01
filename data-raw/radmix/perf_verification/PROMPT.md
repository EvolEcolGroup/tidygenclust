# Prompt for the desktop

Start Claude Code in the directory where the repository is (or should be)
cloned, and paste:

```text
I'm handing over performance work on the R package tidygenclust (GitHub: EvolEcolGroup/tidygenclust) from a Claude instance on my laptop. Everything is on the branch radmix-perf-verify.

1. If the repository isn't cloned here yet, clone git@github.com:EvolEcolGroup/tidygenclust.git. Then git fetch and git switch radmix-perf-verify.
2. Read data-raw/radmix/perf_verification/HANDOVER.md in full, and follow its sections 4-7. Section 2 summarises what was done on the laptop. The scripts are in data-raw/radmix/perf_verification/scripts and were tested on the laptop.

By the end of the day I want:
- verification on this machine of everything done on the laptop: the test suite on the new code, each kernel experiment on its own, the compiler-flag measurements, the combined speed-ups, and the stricter convergence threshold (rtol 1e-8);
- how gt_radmix() scales with the number of cores, for kernels and whole fits, baseline against the combined speed-ups, on the HGDP and on the unpruned HGDP; and why it stops scaling where it does;
- the SVD start (init = "svd") against random starts on the full HGDP (1043 x 130,943), K = 3, 5, 7 and 10, five seeds each: time, iterations, and how often each reaches the best solution;
- a report, desktop_results/REPORT.md, with this machine's numbers next to the laptop's, and the evidence for the three decisions in section 2.4 of the handover.

Ground rules:
- First tell me what the machine is (cores, AVX-512, BLAS, memory) and your plan with time estimates, then go ahead without waiting for me.
- Install R packages as needed (tidypopgen from its dev branch). Don't change system settings (CPU governor, turbo, perf permissions) without asking me.
- Never time anything while another heavy job runs. Run the long phases in the background, and give me a short update as each phase ends.
- Commit the results (CSV, summaries, the report; not the .rds files) to radmix-perf-verify as you go, and push that branch. Don't push any other branch, open PRs, or merge.
- Don't modify the commits on svd-start or kernel-speedups. If you find a bug, show me a reproducible example and put a fix on a new local branch.
```

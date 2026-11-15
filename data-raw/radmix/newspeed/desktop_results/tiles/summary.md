# Fits: /home/andrea/tgc_newspeed_work/results/tiles 

Reference build: `final`. Times are medians over the repeats (n) and seeds; ratio = median of the build / median of the reference (below 1 is faster). Load: the range of the 1-minute load average before the fits. Passes: over the genotypes in the main loop (the start included, EM priming not), counted from the log, before / with the fused pass.

| data | K | threads | max_iter | build | n | median s | ratio | s / iteration | iterations | rejected | passes before / fused | load | max rel. loglik diff vs ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hgdp | 5 | 1 | 8 | final | 3 | 17.58 | 1.000 | 2.198 | 8 | 0 | 41 / 33 | 0.4-0.8 | 0 |
| hgdp | 5 | 1 | 8 | final_q64 | 3 | 17.97 | 1.022 | 2.246 | 8 | 0 | 41 / 33 | 0.6-0.8 | 2.9e-15 |
| hgdp | 5 | 1 | 8 | final_q128 | 3 | 17.50 | 0.996 | 2.188 | 8 | 0 | 41 / 33 | 0.7-0.9 | 3.1e-15 |
| hgdp | 5 | 1 | 8 | final_guided | 3 | 17.88 | 1.017 | 2.235 | 8 | 0 | 41 / 33 | 0.7-0.9 | 0 |
| hgdp | 5 | 1 | 8 | final_tj32 | 3 | 17.95 | 1.021 | 2.244 | 8 | 0 | 41 / 33 | 0.7-0.8 | 8.8e-16 |
| hgdp | 5 | 8 | 8 | final | 3 | 2.38 | 1.000 | 0.298 | 8 | 0 | 41 / 33 | 1.1-2.3 | 0 |
| hgdp | 5 | 8 | 8 | final_q64 | 3 | 2.29 | 0.964 | 0.287 | 8 | 0 | 41 / 33 | 1.5-2.7 | 2.9e-15 |
| hgdp | 5 | 8 | 8 | final_q128 | 3 | 2.23 | 0.936 | 0.279 | 8 | 0 | 41 / 33 | 1.6-2.5 | 3.1e-15 |
| hgdp | 5 | 8 | 8 | final_guided | 3 | 2.53 | 1.062 | 0.316 | 8 | 0 | 41 / 33 | 2.1-2.8 | 0 |
| hgdp | 5 | 8 | 8 | final_tj32 | 3 | 2.39 | 1.003 | 0.298 | 8 | 0 | 41 / 33 | 2.3-3.0 | 8.8e-16 |
| hgdp | 10 | 1 | 8 | final | 3 | 36.81 | 1.000 | 4.601 | 8 | 0 | 41 / 33 | 0.9-1.1 | 0 |
| hgdp | 10 | 1 | 8 | final_q64 | 3 | 36.13 | 0.982 | 4.517 | 8 | 0 | 41 / 33 | 0.9-1.0 | 5.9e-16 |
| hgdp | 10 | 1 | 8 | final_q128 | 3 | 35.84 | 0.974 | 4.480 | 8 | 0 | 41 / 33 | 0.9-1.2 | 4.5e-16 |
| hgdp | 10 | 1 | 8 | final_guided | 3 | 36.58 | 0.994 | 4.572 | 8 | 0 | 41 / 33 | 0.9-1.0 | 0 |
| hgdp | 10 | 1 | 8 | final_tj32 | 3 | 36.59 | 0.994 | 4.574 | 8 | 0 | 41 / 33 | 1.1-1.4 | 4.5e-15 |
| hgdp | 10 | 8 | 8 | final | 3 | 4.96 | 1.000 | 0.620 | 8 | 0 | 41 / 33 | 2.7-3.6 | 0 |
| hgdp | 10 | 8 | 8 | final_q64 | 3 | 4.83 | 0.974 | 0.604 | 8 | 0 | 41 / 33 | 2.9-3.7 | 5.9e-16 |
| hgdp | 10 | 8 | 8 | final_q128 | 3 | 4.79 | 0.966 | 0.599 | 8 | 0 | 41 / 33 | 3.0-3.8 | 4.5e-16 |
| hgdp | 10 | 8 | 8 | final_guided | 3 | 5.31 | 1.070 | 0.663 | 8 | 0 | 41 / 33 | 2.6-3.9 | 0 |
| hgdp | 10 | 8 | 8 | final_tj32 | 3 | 4.98 | 1.004 | 0.622 | 8 | 0 | 41 / 33 | 2.8-4.1 | 4.5e-15 |
| hgdp_full | 5 | 1 | 8 | final | 3 | 85.56 | 1.000 | 10.695 | 8 | 0 | 41 / 33 | 1.1-1.1 | 0 |
| hgdp_full | 5 | 1 | 8 | final_q64 | 3 | 83.95 | 0.981 | 10.494 | 8 | 0 | 41 / 33 | 0.9-1.2 | 8e-16 |
| hgdp_full | 5 | 1 | 8 | final_q128 | 3 | 83.11 | 0.971 | 10.389 | 8 | 0 | 41 / 33 | 1.0-1.1 | 4e-16 |
| hgdp_full | 5 | 1 | 8 | final_guided | 3 | 83.82 | 0.980 | 10.477 | 8 | 0 | 41 / 33 | 1.0-1.1 | 0 |
| hgdp_full | 5 | 1 | 8 | final_tj32 | 3 | 83.88 | 0.980 | 10.485 | 8 | 0 | 41 / 33 | 1.0-1.1 | 8e-16 |
| hgdp_full | 5 | 8 | 8 | final | 3 | 12.35 | 1.000 | 1.543 | 8 | 0 | 41 / 33 | 4.1-5.6 | 0 |
| hgdp_full | 5 | 8 | 8 | final_q64 | 3 | 12.07 | 0.978 | 1.509 | 8 | 0 | 41 / 33 | 4.5-4.8 | 8e-16 |
| hgdp_full | 5 | 8 | 8 | final_q128 | 3 | 11.91 | 0.965 | 1.489 | 8 | 0 | 41 / 33 | 4.7-4.9 | 4e-16 |
| hgdp_full | 5 | 8 | 8 | final_guided | 3 | 12.98 | 1.051 | 1.622 | 8 | 0 | 41 / 33 | 4.8-5.6 | 0 |
| hgdp_full | 5 | 8 | 8 | final_tj32 | 3 | 12.38 | 1.003 | 1.548 | 8 | 0 | 41 / 33 | 5.1-5.3 | 8e-16 |

## Thread invariance (identical log-likelihood histories and Q)

- final hgdp K = 10 seed 1 max_iter 8, threads 1,8: identical
- final_guided hgdp K = 10 seed 1 max_iter 8, threads 1,8: identical
- final_q128 hgdp K = 10 seed 1 max_iter 8, threads 1,8: identical
- final_q64 hgdp K = 10 seed 1 max_iter 8, threads 1,8: identical
- final_tj32 hgdp K = 10 seed 1 max_iter 8, threads 1,8: identical
- final hgdp K = 5 seed 1 max_iter 8, threads 1,8: identical
- final_guided hgdp K = 5 seed 1 max_iter 8, threads 1,8: identical
- final_q128 hgdp K = 5 seed 1 max_iter 8, threads 1,8: identical
- final_q64 hgdp K = 5 seed 1 max_iter 8, threads 1,8: identical
- final_tj32 hgdp K = 5 seed 1 max_iter 8, threads 1,8: identical
- final hgdp_full K = 5 seed 1 max_iter 8, threads 1,8: identical
- final_guided hgdp_full K = 5 seed 1 max_iter 8, threads 1,8: identical
- final_q128 hgdp_full K = 5 seed 1 max_iter 8, threads 1,8: identical
- final_q64 hgdp_full K = 5 seed 1 max_iter 8, threads 1,8: identical
- final_tj32 hgdp_full K = 5 seed 1 max_iter 8, threads 1,8: identical

# Fits: /home/andrea/tgc_newspeed_work/results/iter8 

Reference build: `base`. Times are medians over the repeats (n) and seeds; ratio = median of the build / median of the reference (below 1 is faster). Load: the range of the 1-minute load average before the fits. Passes: over the genotypes in the main loop (the start included, EM priming not), counted from the log, before / with the fused pass.

| data | K | threads | max_iter | build | n | median s | ratio | s / iteration | iterations | rejected | passes before / fused | load | max rel. loglik diff vs ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hgdp | 5 | 1 | 8 | base | 3 | 22.95 | 1.000 | 2.869 | 8 | 0 | 41 / 33 | 0.3-0.9 | 0 |
| hgdp | 5 | 1 | 8 | fused | 3 | 21.41 | 0.933 | 2.676 | 8 | 0 | 41 / 33 | 0.6-0.8 | 0 |
| hgdp | 5 | 1 | 8 | final | 3 | 17.50 | 0.763 | 2.188 | 8 | 0 | 41 / 33 | 0.8-1.0 | 5.8e-16 |
| hgdp | 5 | 1 | 8 | head | 3 | 17.41 | 0.758 | 2.176 | 8 | 0 | 41 / 33 | 0.8-0.9 | 5.8e-16 |
| hgdp | 5 | 4 | 8 | base | 3 | 5.88 | 1.000 | 0.735 | 8 | 0 | 41 / 33 | 1.0-1.8 | 0 |
| hgdp | 5 | 4 | 8 | fused | 3 | 5.55 | 0.944 | 0.694 | 8 | 0 | 41 / 33 | 1.1-1.7 | 0 |
| hgdp | 5 | 4 | 8 | final | 3 | 4.57 | 0.778 | 0.572 | 8 | 0 | 41 / 33 | 1.2-1.7 | 5.8e-16 |
| hgdp | 5 | 4 | 8 | head | 3 | 4.50 | 0.765 | 0.562 | 8 | 0 | 41 / 33 | 1.3-1.8 | 5.8e-16 |
| hgdp | 5 | 8 | 8 | base | 3 | 3.11 | 1.000 | 0.389 | 8 | 0 | 41 / 33 | 2.4-3.0 | 0 |
| hgdp | 5 | 8 | 8 | fused | 3 | 2.93 | 0.941 | 0.366 | 8 | 0 | 41 / 33 | 2.2-3.2 | 0 |
| hgdp | 5 | 8 | 8 | final | 3 | 2.42 | 0.776 | 0.302 | 8 | 0 | 41 / 33 | 2.5-2.9 | 5.8e-16 |
| hgdp | 5 | 8 | 8 | head | 3 | 2.40 | 0.770 | 0.300 | 8 | 0 | 41 / 33 | 2.7-2.7 | 5.8e-16 |
| hgdp | 10 | 1 | 8 | base | 3 | 46.54 | 1.000 | 5.818 | 8 | 0 | 41 / 33 | 0.8-1.0 | 0 |
| hgdp | 10 | 1 | 8 | fused | 3 | 43.30 | 0.930 | 5.412 | 8 | 0 | 41 / 33 | 0.9-1.0 | 0 |
| hgdp | 10 | 1 | 8 | final | 3 | 37.04 | 0.796 | 4.630 | 8 | 0 | 41 / 33 | 0.9-1.0 | 4.5e-16 |
| hgdp | 10 | 1 | 8 | head | 3 | 35.48 | 0.762 | 4.435 | 8 | 0 | 41 / 33 | 0.9-0.9 | 4.5e-16 |
| hgdp | 10 | 4 | 8 | base | 3 | 11.66 | 1.000 | 1.457 | 8 | 0 | 41 / 33 | 1.6-2.4 | 0 |
| hgdp | 10 | 4 | 8 | fused | 3 | 11.06 | 0.949 | 1.383 | 8 | 0 | 41 / 33 | 1.8-2.3 | 0 |
| hgdp | 10 | 4 | 8 | final | 3 | 9.21 | 0.790 | 1.151 | 8 | 0 | 41 / 33 | 1.9-2.4 | 4.5e-16 |
| hgdp | 10 | 4 | 8 | head | 3 | 9.15 | 0.785 | 1.144 | 8 | 0 | 41 / 33 | 2.0-2.3 | 4.5e-16 |
| hgdp | 10 | 8 | 8 | base | 3 | 6.28 | 1.000 | 0.785 | 8 | 0 | 41 / 33 | 2.5-3.7 | 0 |
| hgdp | 10 | 8 | 8 | fused | 3 | 5.96 | 0.949 | 0.745 | 8 | 0 | 41 / 33 | 2.7-4.0 | 0 |
| hgdp | 10 | 8 | 8 | final | 3 | 5.01 | 0.798 | 0.626 | 8 | 0 | 41 / 33 | 2.9-4.0 | 4.5e-16 |
| hgdp | 10 | 8 | 8 | head | 3 | 4.97 | 0.792 | 0.622 | 8 | 0 | 41 / 33 | 3.0-4.0 | 4.5e-16 |
| hgdp_full | 5 | 1 | 8 | base | 3 | 112.41 | 1.000 | 14.051 | 8 | 0 | 41 / 33 | 0.9-1.3 | 0 |
| hgdp_full | 5 | 1 | 8 | fused | 3 | 107.18 | 0.954 | 13.398 | 8 | 0 | 41 / 33 | 0.9-1.0 | 0 |
| hgdp_full | 5 | 1 | 8 | final | 3 | 85.83 | 0.764 | 10.728 | 8 | 0 | 41 / 33 | 1.0-1.0 | 6e-16 |
| hgdp_full | 5 | 1 | 8 | head | 3 | 83.56 | 0.743 | 10.445 | 8 | 0 | 41 / 33 | 0.9-1.1 | 6e-16 |
| hgdp_full | 5 | 4 | 8 | base | 3 | 29.20 | 1.000 | 3.650 | 8 | 0 | 41 / 33 | 2.3-3.2 | 0 |
| hgdp_full | 5 | 4 | 8 | fused | 3 | 27.48 | 0.941 | 3.435 | 8 | 0 | 41 / 33 | 2.9-3.3 | 0 |
| hgdp_full | 5 | 4 | 8 | final | 3 | 22.40 | 0.767 | 2.800 | 8 | 0 | 41 / 33 | 2.9-3.2 | 6e-16 |
| hgdp_full | 5 | 4 | 8 | head | 3 | 22.60 | 0.774 | 2.825 | 8 | 0 | 41 / 33 | 3.0-3.3 | 6e-16 |
| hgdp_full | 5 | 8 | 8 | base | 3 | 16.00 | 1.000 | 2.001 | 8 | 0 | 41 / 33 | 4.0-5.5 | 0 |
| hgdp_full | 5 | 8 | 8 | fused | 3 | 15.07 | 0.942 | 1.884 | 8 | 0 | 41 / 33 | 4.7-5.7 | 0 |
| hgdp_full | 5 | 8 | 8 | final | 3 | 12.45 | 0.778 | 1.557 | 8 | 0 | 41 / 33 | 5.0-5.9 | 6e-16 |
| hgdp_full | 5 | 8 | 8 | head | 3 | 12.47 | 0.779 | 1.559 | 8 | 0 | 41 / 33 | 4.9-6.1 | 6e-16 |

## Thread invariance (identical log-likelihood histories and Q)

- base hgdp K = 10 seed 1 max_iter 8, threads 1,4,8: identical
- final hgdp K = 10 seed 1 max_iter 8, threads 1,4,8: identical
- fused hgdp K = 10 seed 1 max_iter 8, threads 1,4,8: identical
- head hgdp K = 10 seed 1 max_iter 8, threads 1,4,8: identical
- base hgdp K = 5 seed 1 max_iter 8, threads 1,4,8: identical
- final hgdp K = 5 seed 1 max_iter 8, threads 1,4,8: identical
- fused hgdp K = 5 seed 1 max_iter 8, threads 1,4,8: identical
- head hgdp K = 5 seed 1 max_iter 8, threads 1,4,8: identical
- base hgdp_full K = 5 seed 1 max_iter 8, threads 1,4,8: identical
- final hgdp_full K = 5 seed 1 max_iter 8, threads 1,4,8: identical
- fused hgdp_full K = 5 seed 1 max_iter 8, threads 1,4,8: identical
- head hgdp_full K = 5 seed 1 max_iter 8, threads 1,4,8: identical

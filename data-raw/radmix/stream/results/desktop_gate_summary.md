# Fits: results/gate_fits 

Reference build: `base`. Times are medians over the repeats (n) and seeds; ratio = median of the build / median of the reference (below 1 is faster). Load: the range of the 1-minute load average before the fits. Passes: over the genotypes in the main loop (the start included, EM priming not), counted from the log, before / with the fused pass.

| data | K | threads | max_iter | build | n | median s | ratio | s / iteration | iterations | rejected | passes before / fused | load | max rel. loglik diff vs ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eur | 3 | 1 | 1000 | base | 3 | 4.86 | 1.000 | 0.211 | 23 | 3 | 119 / 96 | 0.8-1.0 | 0 |
| eur | 3 | 1 | 1000 | new | 3 | 5.06 | 1.041 | 0.220 | 23 | 3 | 119 / 96 | 0.7-0.9 | 0 |
| eur | 3 | 8 | 1000 | base | 3 | 0.66 | 1.000 | 0.029 | 23 | 3 | 119 / 96 | 4.2-5.5 | 0 |
| eur | 3 | 8 | 1000 | new | 3 | 0.67 | 1.002 | 0.029 | 23 | 3 | 119 / 96 | 3.8-5.0 | 0 |
| eur | 4 | 1 | 1000 | base | 3 | 6.07 | 1.000 | 0.253 | 24 | 3 | 124 / 100 | 0.7-0.7 | 0 |
| eur | 4 | 1 | 1000 | new | 3 | 6.26 | 1.032 | 0.261 | 24 | 3 | 124 / 100 | 0.7-0.7 | 0 |
| eur | 4 | 8 | 1000 | base | 3 | 0.81 | 1.000 | 0.034 | 24 | 3 | 124 / 100 | 2.6-3.5 | 0 |
| eur | 4 | 8 | 1000 | new | 3 | 0.82 | 1.005 | 0.034 | 24 | 3 | 124 / 100 | 2.4-3.6 | 0 |
| hgdp | 5 | 1 | 1000 | base | 3 | 41.61 | 1.000 | 1.981 | 21 | 0 | 106 / 85 | 0.0-0.9 | 0 |
| hgdp | 5 | 1 | 1000 | new | 3 | 40.95 | 0.984 | 1.950 | 21 | 0 | 106 / 85 | 0.5-1.1 | 0 |
| hgdp | 5 | 8 | 1000 | base | 3 | 5.09 | 1.000 | 0.243 | 21 | 0 | 106 / 85 | 0.6-2.1 | 0 |
| hgdp | 5 | 8 | 1000 | new | 3 | 5.14 | 1.009 | 0.245 | 21 | 0 | 106 / 85 | 1.1-2.4 | 0 |
| hgdp | 10 | 1 | 1000 | base | 3 | 132.75 | 1.000 | 3.588 | 37 | 14 | 200 / 163 | 0.9-1.0 | 0 |
| hgdp | 10 | 1 | 1000 | new | 3 | 133.36 | 1.005 | 3.604 | 37 | 14 | 200 / 163 | 0.9-1.0 | 0 |
| hgdp | 10 | 8 | 1000 | base | 3 | 17.75 | 1.000 | 0.480 | 37 | 14 | 200 / 163 | 2.6-4.9 | 0 |
| hgdp | 10 | 8 | 1000 | new | 3 | 17.93 | 1.010 | 0.485 | 37 | 14 | 200 / 163 | 3.5-5.1 | 0 |

## Thread invariance (identical log-likelihood histories and Q)

- base hgdp K = 10 seed 1 max_iter 1000, threads 1,8: identical
- new hgdp K = 10 seed 1 max_iter 1000, threads 1,8: identical
- base eur K = 3 seed 1 max_iter 1000, threads 1,8: identical
- new eur K = 3 seed 1 max_iter 1000, threads 1,8: identical
- base eur K = 4 seed 1 max_iter 1000, threads 1,8: identical
- new eur K = 4 seed 1 max_iter 1000, threads 1,8: identical
- base hgdp K = 5 seed 1 max_iter 1000, threads 1,8: identical
- new hgdp K = 5 seed 1 max_iter 1000, threads 1,8: identical

# Passes over the genotypes: results/gate_kernels 

Median seconds of one pass (medians over runs of medians); iteration = the passes of an accepted iteration: before, 2 update_q + 2 update_p + loglik; fused, update_q + update_q_loglik + 2 update_p. Ratio: to the iteration of the reference build `base` (before).

| data | K | threads | view | build | loglik | update_q | update_p | update_q_loglik | iteration before | iteration fused | ratio fused |
|---|---|---|---|---|---|---|---|---|---|---|---|
| eur | 3 | 1 | none | base | 0.054 | 0.038 | 0.043 | 0.057 | 0.216 | 0.181 | 0.838 |
| eur | 3 | 1 | none | new | 0.058 | 0.037 | 0.045 | 0.058 | 0.222 | 0.185 | 0.856 |
| eur | 4 | 1 | none | base | 0.058 | 0.044 | 0.054 | 0.065 | 0.254 | 0.217 | 0.854 |
| eur | 4 | 1 | none | new | 0.062 | 0.045 | 0.057 | 0.066 | 0.266 | 0.225 | 0.886 |
| hgdp | 5 | 1 | none | base | 0.418 | 0.359 | 0.423 | 0.495 | 1.982 | 1.700 | 0.858 |
| hgdp | 5 | 1 | none | new | 0.414 | 0.355 | 0.404 | 0.502 | 1.932 | 1.665 | 0.840 |
| hgdp | 10 | 1 | none | base | 0.570 | 0.647 | 0.782 | 0.796 | 3.428 | 3.007 | 0.877 |
| hgdp | 10 | 1 | none | new | 0.570 | 0.642 | 0.797 | 0.803 | 3.448 | 3.039 | 0.887 |
| eur | 3 | 8 | none | base | 0.008 | 0.005 | 0.006 | 0.008 | 0.030 | 0.025 | 0.833 |
| eur | 3 | 8 | none | new | 0.008 | 0.005 | 0.006 | 0.008 | 0.030 | 0.025 | 0.833 |
| eur | 4 | 8 | none | base | 0.008 | 0.007 | 0.007 | 0.009 | 0.036 | 0.030 | 0.833 |
| eur | 4 | 8 | none | new | 0.008 | 0.006 | 0.007 | 0.009 | 0.034 | 0.029 | 0.806 |
| hgdp | 5 | 8 | none | base | 0.054 | 0.049 | 0.052 | 0.064 | 0.256 | 0.217 | 0.848 |
| hgdp | 5 | 8 | none | new | 0.052 | 0.047 | 0.053 | 0.065 | 0.252 | 0.218 | 0.852 |
| hgdp | 10 | 8 | none | base | 0.072 | 0.083 | 0.098 | 0.102 | 0.434 | 0.381 | 0.878 |
| hgdp | 10 | 8 | none | new | 0.069 | 0.083 | 0.100 | 0.104 | 0.435 | 0.387 | 0.892 |


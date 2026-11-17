# Fits: /home/andrea/tgc_bench_article/results/port_timing 

Reference build: `dev`. Times are medians over the repeats (n) and seeds; ratio = median of the build / median of the reference (below 1 is faster). Load: the range of the 1-minute load average before the fits. Passes: over the genotypes in the main loop (the start included, EM priming not), counted from the log, before / with the fused pass.

| data | K | threads | max_iter | build | n | median s | ratio | s / iteration | iterations | rejected | passes before / fused | load | max rel. loglik diff vs ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eur | 3 | 1 | 1000 | dev | 3 | 6.26 | 1.000 | 0.272 | 23 | 3 | 119 / 96 | 0.9-1.0 | 0 |
| eur | 3 | 1 | 1000 | port | 3 | 5.02 | 0.802 | 0.218 | 23 | 3 | 119 / 96 | 0.8-0.9 | 2.4e-16 |
| eur | 3 | 8 | 1000 | dev | 3 | 0.89 | 1.000 | 0.038 | 23 | 3 | 119 / 96 | 4.2-6.3 | 0 |
| eur | 3 | 8 | 1000 | port | 3 | 0.69 | 0.782 | 0.030 | 23 | 3 | 119 / 96 | 4.2-5.4 | 2.4e-16 |
| eur | 4 | 1 | 1000 | dev | 3 | 7.52 | 1.000 | 0.313 | 24 | 3 | 124 / 100 | 0.8-0.8 | 0 |
| eur | 4 | 1 | 1000 | port | 3 | 6.22 | 0.827 | 0.259 | 24 | 3 | 124 / 100 | 0.7-0.8 | 2.4e-16 |
| eur | 4 | 8 | 1000 | dev | 3 | 1.06 | 1.000 | 0.044 | 24 | 3 | 124 / 100 | 2.8-3.8 | 0 |
| eur | 4 | 8 | 1000 | port | 3 | 0.84 | 0.795 | 0.035 | 24 | 3 | 124 / 100 | 2.6-3.5 | 2.4e-16 |
| hgdp | 5 | 1 | 1000 | dev | 3 | 50.31 | 1.000 | 2.396 | 21 | 0 | 106 / 85 | 0.1-0.9 | 0 |
| hgdp | 5 | 1 | 1000 | port | 3 | 40.02 | 0.795 | 1.906 | 21 | 0 | 106 / 85 | 0.6-0.9 | 2.9e-16 |
| hgdp | 5 | 8 | 1000 | dev | 3 | 6.68 | 1.000 | 0.318 | 21 | 0 | 106 / 85 | 0.7-2.2 | 0 |
| hgdp | 5 | 8 | 1000 | port | 3 | 5.28 | 0.790 | 0.251 | 21 | 0 | 106 / 85 | 1.2-3.1 | 2.9e-16 |
| hgdp | 10 | 1 | 1000 | dev | 3 | 183.32 | 1.000 | 4.955 | 37 | 15 | 201 / 164 | 1.0-1.2 | 0 |
| hgdp | 10 | 1 | 1000 | port | 3 | 168.01 | 0.916 | 4.541 | 37 | 15 | 201 / 164 | 1.0-1.1 | 7.4e-16 |
| hgdp | 10 | 8 | 1000 | dev | 3 | 25.79 | 1.000 | 0.697 | 37 | 15 | 201 / 164 | 3.2-5.8 | 0 |
| hgdp | 10 | 8 | 1000 | port | 3 | 21.89 | 0.849 | 0.592 | 37 | 15 | 201 / 164 | 4.4-6.1 | 7.4e-16 |

## Thread invariance (identical log-likelihood histories and Q)

- dev hgdp K = 10 seed 1 max_iter 1000, threads 1,8: identical
- port hgdp K = 10 seed 1 max_iter 1000, threads 1,8: identical
- dev eur K = 3 seed 1 max_iter 1000, threads 1,8: identical
- port eur K = 3 seed 1 max_iter 1000, threads 1,8: identical
- dev eur K = 4 seed 1 max_iter 1000, threads 1,8: identical
- port eur K = 4 seed 1 max_iter 1000, threads 1,8: identical
- dev hgdp K = 5 seed 1 max_iter 1000, threads 1,8: identical
- port hgdp K = 5 seed 1 max_iter 1000, threads 1,8: identical

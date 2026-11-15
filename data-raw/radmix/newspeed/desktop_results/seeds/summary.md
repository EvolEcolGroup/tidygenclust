# Fits: /home/andrea/tgc_newspeed_work/results/seeds 

Reference build: `base`. Times are medians over the repeats (n) and seeds; ratio = median of the build / median of the reference (below 1 is faster). Load: the range of the 1-minute load average before the fits. Passes: over the genotypes in the main loop (the start included, EM priming not), counted from the log, before / with the fused pass.

| data | K | threads | max_iter | build | n | median s | ratio | s / iteration | iterations | rejected | passes before / fused | load | max rel. loglik diff vs ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eur | 3 | 8 | 1000 | base | 3 | 0.79 | 1.000 | 0.040 | 20,23 | 0,2,3 | 101,103,119 / 81,83,96 | 4.8-5.7 | 0 |
| eur | 3 | 8 | 1000 | final | 3 | 0.62 | 0.777 | 0.031 | 20,23 | 0,2,3 | 101,103,119 / 81,83,96 | 4.0-4.7 | 8.4e-16 |
| eur | 4 | 8 | 1000 | base | 3 | 1.18 | 1.000 | 0.048 | 24,53 | 2,3,23 | 123,124,289 / 99,100,236 | 3.3-3.6 | 0 |
| eur | 4 | 8 | 1000 | final | 3 | 0.90 | 0.765 | 0.038 | 24,53 | 2,3,23 | 123,124,289 / 99,100,236 | 3.1-3.4 | 7.2e-16 |
| hgdp | 5 | 8 | 1000 | base | 3 | 7.32 | 1.000 | 0.345 | 21,27 | 0 | 106,136 / 85,109 | 0.6-1.9 | 0 |
| hgdp | 5 | 8 | 1000 | final | 3 | 5.72 | 0.782 | 0.272 | 21,27 | 0 | 106,136 / 85,109 | 2.6-3.1 | 3.5e-15 |
| hgdp | 10 | 8 | 1000 | base | 3 | 26.67 | 1.000 | 0.756 | 24,35,37 | 0,8,15 | 121,184,201 / 97,149,164 | 2.6-5.0 | 0 |
| hgdp | 10 | 8 | 1000 | final | 3 | 21.29 | 0.798 | 0.608 | 24,35,37 | 0,8,15 | 121,184,201 / 97,149,164 | 5.4-5.7 | 2.7e-15 |

## Thread invariance (identical log-likelihood histories and Q)


# Fits: /media/andrea/annex/tmp/claude/tidygenclust/newspeed/out/kit_tiny 

Reference build: `base`. Times are medians over the repeats (n) and seeds; ratio = median of the build / median of the reference (below 1 is faster). Load: the range of the 1-minute load average before the fits. Passes: over the genotypes in the main loop (the start included, EM priming not), counted from the log, before / with the fused pass.

| data | K | threads | max_iter | build | n | median s | ratio | s / iteration | iterations | rejected | passes before / fused | load | max rel. loglik diff vs ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hgdp_quick | 3 | 1 | 1 | base | 2 | 0.63 | 1.000 | 0.633 | 1 | 0 | 6 / 5 | 5.7-5.9 | 0 |
| hgdp_quick | 3 | 1 | 1 | sg5ll | 2 | 0.53 | 0.840 | 0.531 | 1 | 0 | 6 / 5 | 5.7-5.9 | 0 |
| hgdp_quick | 3 | 2 | 1 | base | 2 | 0.63 | 1.000 | 0.628 | 1 | 0 | 6 / 5 | 5.3-5.5 | 0 |
| hgdp_quick | 3 | 2 | 1 | sg5ll | 2 | 0.54 | 0.856 | 0.538 | 1 | 0 | 6 / 5 | 5.3-5.5 | 0 |

## Thread invariance (identical log-likelihood histories and Q)

- base hgdp_quick K = 3 seed 1 max_iter 1, threads 1,2: identical
- sg5ll hgdp_quick K = 3 seed 1 max_iter 1, threads 1,2: identical

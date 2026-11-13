# Fits: /media/andrea/annex/tmp/claude/tidygenclust/newspeed/out/whole_eur 

Reference build: `base`. Times are medians over the repeats (n) and seeds; ratio = median of the build / median of the reference (below 1 is faster). Load: the range of the 1-minute load average before the fits. Passes: over the genotypes in the main loop (the start included, EM priming not), counted from the log, before / with the fused pass.

| data | K | threads | max_iter | build | n | median s | ratio | s / iteration | iterations | rejected | passes before / fused | load | max rel. loglik diff vs ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eur | 3 | 1 | 1000 | base | 2 | 26.31 | 1.000 | 1.144 | 23 | 3 | 119 / 96 | 3.9-4.7 | 0 |
| eur | 3 | 1 | 1000 | sg5ll | 2 | 21.59 | 0.820 | 0.939 | 23 | 3 | 119 / 96 | 4.0-4.6 | 1.2e-16 |
| eur | 4 | 1 | 1000 | base | 2 | 30.43 | 1.000 | 1.268 | 24 | 3 | 124 / 100 | 4.4-4.5 | 0 |
| eur | 4 | 1 | 1000 | sg5ll | 2 | 26.81 | 0.881 | 1.117 | 24 | 3 | 124 / 100 | 4.1-4.2 | 3.6e-16 |

## Thread invariance (identical log-likelihood histories and Q)


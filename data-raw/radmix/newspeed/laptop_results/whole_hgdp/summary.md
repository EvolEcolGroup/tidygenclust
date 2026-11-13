# Fits: /media/andrea/annex/tmp/claude/tidygenclust/newspeed/out/whole_hgdp 

Reference build: `base`. Times are medians over the repeats (n) and seeds; ratio = median of the build / median of the reference (below 1 is faster). Load: the range of the 1-minute load average before the fits. Passes: over the genotypes in the main loop (the start included, EM priming not), counted from the log, before / with the fused pass.

| data | K | threads | max_iter | build | n | median s | ratio | s / iteration | iterations | rejected | passes before / fused | load | max rel. loglik diff vs ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hgdp | 5 | 1 | 1000 | base | 3 | 225.89 | 1.000 | 10.449 | 21,27 | 0 | 106,136 / 85,109 | 4.4-5.1 | 0 |
| hgdp | 5 | 1 | 1000 | sg5ll | 3 | 175.56 | 0.777 | 7.625 | 21,27 | 0 | 106,136 / 85,109 | 2.9-4.4 | 3.5e-15 |

## Thread invariance (identical log-likelihood histories and Q)


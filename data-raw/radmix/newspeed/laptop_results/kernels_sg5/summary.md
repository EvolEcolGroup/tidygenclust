# Passes over the genotypes: /media/andrea/annex/tmp/claude/tidygenclust/newspeed/out/kernels_sg5 

Median seconds of one pass (medians over runs of medians); iteration = the passes of an accepted iteration: before, 2 update_q + 2 update_p + loglik; fused, update_q + update_q_loglik + 2 update_p. Ratio: to the iteration of the reference build `base` (before).

| data | K | threads | view | build | loglik | update_q | update_p | update_q_loglik | iteration before | iteration fused | ratio fused |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hgdp | 5 | 1 | none | base | 2.125 | 1.843 | 2.138 |  | 10.086 |  |  |
| hgdp | 5 | 1 | none | fused | 2.200 | 2.082 | 2.302 | 3.308 | 10.968 | 9.993 | 0.991 |
| hgdp | 5 | 1 | none | sg5 | 1.828 | 1.652 | 1.892 | 2.768 | 8.916 | 8.204 | 0.813 |
| hgdp | 5 | 1 | none | sg5ll | 1.752 | 1.683 | 1.938 | 2.269 | 8.995 | 7.829 | 0.776 |


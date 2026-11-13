# Passes over the genotypes: /media/andrea/annex/tmp/claude/tidygenclust/newspeed/out/kernels_ll 

Median seconds of one pass (medians over runs of medians); iteration = the passes of an accepted iteration: before, 2 update_q + 2 update_p + loglik; fused, update_q + update_q_loglik + 2 update_p. Ratio: to the iteration of the reference build `base` (before).

| data | K | threads | view | build | loglik | update_q | update_p | update_q_loglik | iteration before | iteration fused | ratio fused |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hgdp | 5 | 1 | none | base | 2.255 | 2.076 | 2.327 |  | 11.060 |  |  |
| hgdp | 5 | 1 | none | fused | 2.373 | 2.138 | 2.485 | 3.756 | 11.619 | 10.864 | 0.982 |
| hgdp | 5 | 1 | none | ll2 | 3.540 | 1.930 | 2.554 | 4.910 | 12.508 | 11.948 | 1.080 |
| hgdp | 5 | 1 | none | ll2inline | 3.486 | 2.197 | 3.181 | 4.894 | 14.242 | 13.454 | 1.216 |


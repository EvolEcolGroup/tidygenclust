# Passes over the genotypes: /home/andrea/tgc_newspeed_work/results/kernels 

Median seconds of one pass (medians over runs of medians); iteration = the passes of an accepted iteration: before, 2 update_q + 2 update_p + loglik; fused, update_q + update_q_loglik + 2 update_p. Ratio: to the iteration of the reference build `base` (before).

| data | K | threads | view | build | loglik | update_q | update_p | update_q_loglik | iteration before | iteration fused | ratio fused |
|---|---|---|---|---|---|---|---|---|---|---|---|
| eur | 4 | 1 | none | base | 0.072 | 0.051 | 0.065 |  | 0.304 |  |  |
| eur | 4 | 1 | none | final | 0.053 | 0.045 | 0.056 | 0.065 | 0.255 | 0.222 | 0.730 |
| eur | 4 | 1 | none | final_llnoinl | 0.053 | 0.044 | 0.057 | 0.083 | 0.255 | 0.241 | 0.793 |
| eur | 4 | 1 | none | fused | 0.073 | 0.051 | 0.063 | 0.100 | 0.301 | 0.277 | 0.911 |
| eur | 4 | 1 | none | inline | 0.072 | 0.045 | 0.060 | 0.096 | 0.282 | 0.261 | 0.859 |
| hgdp | 5 | 1 | none | base | 0.507 | 0.440 | 0.484 |  | 2.355 |  |  |
| hgdp | 5 | 1 | none | final | 0.368 | 0.387 | 0.432 | 0.521 | 2.006 | 1.772 | 0.752 |
| hgdp | 5 | 1 | none | final_llnoinl | 0.379 | 0.375 | 0.413 | 0.646 | 1.955 | 1.847 | 0.784 |
| hgdp | 5 | 1 | none | fused | 0.516 | 0.435 | 0.480 | 0.775 | 2.346 | 2.170 | 0.921 |
| hgdp | 5 | 1 | none | inline | 0.511 | 0.380 | 0.415 | 0.728 | 2.101 | 1.938 | 0.823 |
| hgdp | 10 | 1 | none | base | 0.658 | 0.998 | 1.144 |  | 4.942 |  |  |
| hgdp | 10 | 1 | none | final | 0.469 | 0.900 | 0.999 | 1.082 | 4.267 | 3.980 | 0.805 |
| hgdp | 10 | 1 | none | final_llnoinl | 0.469 | 0.884 | 0.945 | 1.197 | 4.127 | 3.971 | 0.804 |
| hgdp | 10 | 1 | none | fused | 0.662 | 0.997 | 1.125 | 1.344 | 4.906 | 4.591 | 0.929 |
| hgdp | 10 | 1 | none | inline | 0.649 | 0.903 | 0.954 | 1.280 | 4.363 | 4.091 | 0.828 |


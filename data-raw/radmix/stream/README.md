# Benchmark of streaming in gt_radmix()

Scripts and results of `dev_notes/radmix_streaming.md` (run as root in a
cloud container: they use cgroup v1 memory and blkio limits and
`drop_caches`).

- `make_data.R`: simulated genotypes (an FBM.code256).
- `bench_fit.R`: one fit of the C++ loop, with the bytes read and a hash of
  the estimates; modes `base` (the build of dev), `FALSE`, `TRUE`,
  `TRUE_nopf` (no prefetch), `TRUE_nozz` (no alternating directions).
- `bench_stream.sh`: `bench_fit.R` from a cold page cache under a memory
  cap and a read throttle.
- `run_ooc.sh`: the out-of-core matrix (`results/out_of_core.csv`; the two
  `base` rows of `big_a:60` come from a first run of the matrix).
- `bench_kernels.R`: the time of each pass, warm cache
  (`results/conventional_passes.txt`; full fits in
  `results/conventional_fits.csv`).
- `bench_passes.R`, `bench_passes.sh`, `run_passes.sh`, `summarise_passes.R`:
  step 1, the cost of each pass of dev (and of streaming) out of core
  (`results/passes.csv`).

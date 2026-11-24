# Benchmark of streaming in gt_radmix()

Scripts and results of `dev_notes/radmix_streaming.md`. The results were
measured as root in a cloud container (cgroup v1 memory and blkio limits).
`run_capped.sh` caps the memory (page cache included) and the read
throughput of the disk either that way or, on a machine with systemd and
cgroup v2 (the desktop), with `sudo systemd-run --scope -p MemoryMax=...
-p IOReadBandwidthMax=...` (untested: written in the container); both drop
the page cache first. Paths come from the environment: `BENCH_DATA` (the
data of `make_data.R`), `LIB_BASE` and `LIB_NEW` (the builds of dev and of
the branch), `N_CORES`, `FILE_MB`; the memory caps add a fraction of the
file to the peak memory of the R script, measured by `measure_peak.sh`
into `anon_<prefix>.txt` (fits) and `anonp_<prefix>.txt` (passes).

- `run_capped.sh`: runs a command under the caps, from a cold page cache.
- `measure_peak.sh`: the peak anonymous memory of an R script.
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
- `run_ooc_third.sh`: data about 3 times the memory
  (`results/third_of_memory.csv`; repeats of the runs that were disturbed,
  with the user and system CPU time, in `results/third_of_memory_repeats.csv`).
- The cap of the prefetch of the Q side: `results/prefetch_cap.csv` and
  `results/prefetch_cap_alternated.csv` (the capped and uncapped builds,
  alternated).

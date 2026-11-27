# Streaming data larger than the memory in `gt_radmix()`

How `gt_radmix()` reads a backing file that does not fit in memory, what
was measured, and what is left. The code is in `src/radmix_geno_access.h`
("Streaming larger-than-memory data": `stream()`, `chunk_loci()`,
`direction()`, `prefetch()`), `src/radmix_kernels.cpp` ("Streaming":
`Prefetcher`, `q_side_stream()`, the streaming branches of `p_side()` and
`sum_slices()`) and `R/radmix_stream.R` (when to stream). The benchmark
scripts are in `data-raw/radmix/stream/`.

## The problem

The genotypes (or likelihoods) are never copied into memory: every pass
reads the memory-mapped backing file, so data larger than the memory never
exhaust it, and the operating system reads the file again as needed. But
the usual layout of the work reads it many times per pass once it does not
fit in the page cache:

- Q side with one segment (8,192 individuals or more): each block of 64
  individuals is streamed over all the loci by one thread, and then the
  next block starts again at the first locus. With `nbi` blocks and `n`
  threads, a Q-side pass reads the file about `nbi / n` times (40 times for
  10,000 individuals on 4 threads): a fit of 10,000 x 150,000 genotypes
  (1.5 GB) with 67% of the file in the page cache read 294 GB from disk and
  took 518 s, against 25 s in memory.
- Q side with segments (fewer individuals): the same within a segment, so a
  segment that does not fit in the cache is read many times.
- Every pass is a forward scan: with a file a little larger than the cache,
  each pass evicts the pages that the next one reads first (LRU), and
  almost nothing is read from the cache.

## Step 1: what the passes of dev cost out of core

Before the design below was benchmarked as a whole, each pass of `dev` was
measured on its own (`bench_passes.R`, `run_passes.sh`; summary by
`summarise_passes.R`; data in `results/passes.csv`): the passes of a fit in
the order of the fitting loop (EM, the fused Q update, then P, Q, P, fused
Q per iteration), each with its time, the bytes read from disk, the read
throughput and the CPU used by all the threads, in memory and with 60% or
85% of the file cacheable (the cap is the peak memory of the script plus
that fraction), on the disk as it is and throttled to 500 MB/s. Medians
from the third pass on, 4 threads, K = 4:

| data | cache | disk | pass | dev: time (s), reads of the file, CPU | in memory (s) |
|---|---|---|---|---|---|
| 10,000 x 150,000 | 60% | 2 GB/s | Q update | 103, 40.0, 28% | 2.6 |
| | | | fused Q update | 105, 40.0, 29% | 3.5 |
| | | | P update | 3.2, 1.0, 82% | 2.0 |
| | 60% | 500 MB/s | Q update | 120, 40.0, 24% | |
| | | | P update | 3.0, 1.0, 80% | |
| | 85% | 2 GB/s | Q update | 103, 40.0, 28% | |
| | | | P update | 3.2, 1.0, 83% | |
| 1,000 x 1,500,000 | 60% | 2 GB/s | Q update | 3.2, 0.96, 95% | 2.4 |
| | | | fused Q update | 4.3, 1.0, 94% | 2.8 |
| | | | P update | 3.0, 1.0, 92% | 1.8 |
| | 60% | 500 MB/s | Q update | 4.7, 0.92, 68% | |
| | | | P update | 3.0, 0.92, 84% | |
| | 85% | 2 GB/s | Q update | 2.8, 0.66, 91% | |
| | | | P update | 2.4, 0.41, 92% | |

- In memory every pass is limited by computation (95% of 4 cores), and
  processes the file at 0.4-0.7 GB/s: a pass that reads the whole file
  once keeps up with a disk faster than that (on 4 cores; twice as fast on
  8), and waits for a slower one.
- With 8,192 individuals or more (one Q segment), the Q-side passes of dev
  read the file 40 times (`nbi / n_cores` = 157 / 4), at the speed of the
  disk (477 MB/s on the 500 MB/s disk, 550 MB/s on the faster one, whose
  page faults do not reach the 2 GB/s of a sequential read), with 24-29% of
  the CPU: 35-45 times their time in memory, whatever the fraction cached
  below 100%. The P-side passes read the file once, at 80-83% of the CPU,
  1.5 times their time in memory. The Q side is the whole problem.
- With fewer individuals (segments), dev re-reads each segment from the
  cache, and every pass reads the file at most once: 1.0-1.7 times the time
  in memory, still limited by computation on the faster disk (91-95%), and
  partly by the disk at 500 MB/s (68% for the Q update).
- The page cache is not LRU: dev's forward passes re-read only 41-66% of
  the file with 85% of it cacheable (a pure LRU would make them read all of
  it). Linux keeps pages used twice on its active list, which a forward
  scan does not evict.
- Streaming, measured the same way (one call per pass, so one view per
  pass: these runs do not alternate directions, as each view starts
  forwards), takes 2.5-4.8 s per pass at 82-94% of the CPU, reading
  0.84-1.0 of the file per pass with 10,000 individuals. With 1,000
  individuals and 85% cacheable, it read the whole file per pass, where dev
  read 41-66% of it, and was 15-35% slower than dev: without alternating
  directions, the chunks and the prefetch defeat the kernel's reuse of the
  cache. Within a fit (one view, alternating), streaming read less than dev
  there (7.7 against 10.6 GB, below).

## The design

A view in streaming mode (`data$stream`) makes every pass walk the loci in
chunks of about `chunk_bytes` bytes of the file (whole tiles of 64 loci),
one parallel region per chunk, all the threads working on the same chunk:

- **Q side, one segment**: each block of individuals is accumulated over
  a chunk, then over the next, from where it stopped. Its accumulators (and
  its log-likelihood) are kept between chunks (I (K + K (K + 1) / 2)
  doubles), copied into thread-local buffers for each chunk (accumulating
  in place, in one large array, made the passes 15-40% slower). Chunks have
  an even number of loci, so the loci are paired as in one run: the sums
  are bit-identical.
- **Q side, segments**: a chunk is a group of whole segments; the tasks are
  those of the usual layout. A segment can be much larger than
  `chunk_bytes` (187 MB against 32 MB chunks for 1,000 x 1,500,000 with
  500 MB of memory), so only the first `chunk_bytes` of the next chunk are
  prefetched (see "The prefetch of the Q side", below).
- **P side, log-likelihood, deviance**: blocks of loci are independent (or
  keep their own partial sums), so the chunks can be taken in any order.
- **Prefetch** (step 2 of the plan): before computing a chunk, the main
  thread starts a helper thread (`Prefetcher`, a `std::thread`) that calls
  `madvise(MADV_WILLNEED)` on the columns of the next chunk, merged into
  runs of adjacent columns, in pieces of 128 KB. Linux reads only the
  start of a larger range: at most the larger of the disk's read-ahead and
  its largest request (4 MB in the container: a 256 MB call left 1,024 of
  65,536 pages resident; 128 KB on the desktop's NVMe, where the first
  pieces of 2 MB read 1/16 of each chunk; see "Desktop"), whereas each
  call starts its reading and returns at once (in the container, 256 MB
  advised in 17 ms, all resident 0.5 s later). On cached pages it
  costs nothing (0.04 ms per 256 MB). A source of genotype likelihoods can
  give a prefetch function (optional element `prefetch` of its list, with
  `locus_bytes`; tidypopgl's does).
- **Alternating directions** (step 3): the P side, the log-likelihood and
  the deviance take the chunks in reverse order after a forward pass, and in
  order after a backward one (`direction()`), so that each starts with the
  chunks the previous pass read last, still in the cache. The Q side always
  runs forwards (its sums are accumulated over the loci in order). In a fit,
  Q-side and P-side passes alternate (EM: Q then P; each iteration: P, Q, P
  and the fused Q of the next point), so every pass starts where the last
  one ended, except after a rejected quasi-Newton point (two Q passes in a
  row). Within a chunk the loci are read forwards: the kernel reads ahead
  of forward page faults only (backward reads within chunks doubled the
  page faults), and with the chunk order reversed it stops reading ahead
  even of forward faults, so the alternation only works with the prefetch
  (92% of the file cacheable: 929,000 waiting faults without the prefetch,
  3,900 with it).

The results are bit-identical to those of the usual layout
(`tests/testthat/test_radmix_stream.R`), whatever the size of the chunks
and the number of threads, so whether a fit streams does not matter for
reproducibility.

**When**: `radmix_stream_mode()` streams when the bytes of the file read by
a pass (whole columns of the FBM, or `locus_bytes` per locus of a source)
exceed 0.8 times the memory available (`radmix_mem_available()`: on Linux,
the smaller of `MemAvailable` and the room under the limit of the cgroup,
v1 or v2; `Inf` elsewhere, so it never streams on its own on macOS or
Windows). The option `tidygenclust.radmix_stream` forces it on or off.
Chunks are a sixteenth of the memory available, between 16 and 256 MB
(`tidygenclust.radmix_chunk_mb`): the chunk being computed and the one
being prefetched take memory from the cache, and 256 MB chunks with 860 MB
of cache read 7.5 GB where 32-64 MB chunks read 5.0-5.3 GB, in the same
time.

## Benchmark

Container with 4 cores, 15 GB. Simulated genotypes (`make_data.R`: 3
populations, 1% missing), 10,000 x 150,000 (one Q segment) and 1,000 x
1,500,000 (8 segments), 1.5 GB each. One fit (`bench_fit.R`: K = 4, one EM
step, 2 iterations, about 10 passes), from a cold page cache, under a cgroup
v1 memory cap of the R process plus a fraction of the file (the caps were
the peak memory of the fit plus 100 MB plus 60% or 85% of the file, which
left about 67% and 92% of the file to the page cache)
(`bench_stream.sh`, `run_ooc.sh`), on the disk as it is (about 2 GB/s,
backed by the host's cache) or throttled with the blkio cgroup (500 and
150 MB/s). "dev" is the build of `dev` before streaming; the other rows
are streaming with everything on, without the prefetch, and without
alternating directions. Every run gave the same estimates (hash of Q, P and
the log-likelihoods). "Faults" are major page faults (that waited for the
disk); "read" is `read_bytes` of `/proc/self/io`. In memory (no cap), the
same fit takes 25 s (10,000 x 150,000) and 30 s (1,000 x 1,500,000), and
reads the file once (1.4 GB). Chunks of 64 MB.

10,000 x 150,000:

| cache | disk | version | time (s) | read (GB) | faults |
|---|---|---|---|---|---|
| 67% | 2 GB/s | dev | 518 | 294 | 4,024,000 |
| 67% | 2 GB/s | streaming | 27.3 | 8.2 | 12,000 |
| 67% | 2 GB/s | no prefetch | 36.5 | 8.1 | 147,000 |
| 67% | 2 GB/s | no alternation | 31.4 | 15.1 | 14,000 |
| 67% | 500 MB/s | dev | 618 | 294 | 4,039,000 |
| 67% | 500 MB/s | streaming | 29.4 | 8.0 | 7,000 |
| 67% | 500 MB/s | no prefetch | 39.1 | 8.2 | 185,000 |
| 67% | 500 MB/s | no alternation | 34.4 | 14.1 | 13,000 |
| 67% | 150 MB/s | streaming | 67.7 | 8.3 | 7,000 |
| 67% | 150 MB/s | no prefetch | 65.0 | 7.2 | 140,000 |
| 92% | 2 GB/s | dev | 513 | 295 | 3,997,000 |
| 92% | 2 GB/s | streaming | 27.7 | 5.1 | 3,900 |
| 92% | 2 GB/s | no prefetch | 57.5 | 6.5 | 929,000 |
| 92% | 2 GB/s | no alternation | 30.4 | 15.9 | 14,000 |
| 92% | 500 MB/s | dev | 618 | 294 | 4,031,000 |
| 92% | 500 MB/s | streaming | 29.3 | 5.4 | 4,100 |
| 92% | 500 MB/s | no prefetch | 45.4 | 4.2 | 590,000 |
| 92% | 500 MB/s | no alternation | 31.4 | 11.1 | 21,000 |
| 92% | 150 MB/s | streaming | 50.2 | 4.8 | 6,200 |
| 92% | 150 MB/s | no prefetch | 56.0 | 4.6 | 542,000 |

1,000 x 1,500,000 (8 segments of 187 MB, each a chunk of the Q side):

| cache | disk | version | time (s) | read (GB) | faults |
|---|---|---|---|---|---|
| 67% | 2 GB/s | dev | 31.9 | 13.8 | 18,000 |
| 67% | 2 GB/s | streaming | 29.9 | 8.1 | 1,600 |
| 67% | 2 GB/s | no prefetch | 31.3 | 8.0 | 11,000 |
| 67% | 2 GB/s | no alternation | 33.3 | 15.5 | 6,000 |
| 67% | 500 MB/s | dev | 38.7 | 13.3 | 19,000 |
| 67% | 500 MB/s | streaming | 34.9 | 8.1 | 1,100 |
| 67% | 500 MB/s | no prefetch | 35.6 | 6.6 | 10,000 |
| 67% | 500 MB/s | no alternation | 39.5 | 14.7 | 3,800 |
| 67% | 150 MB/s | streaming | 74.2 | 8.5 | 1,600 |
| 67% | 150 MB/s | no prefetch | 68.8 | 7.1 | 11,000 |
| 92% | 2 GB/s | dev | 30.7 | 10.6 | 13,000 |
| 92% | 2 GB/s | streaming | 29.0 | 7.7 | 3,600 |
| 92% | 2 GB/s | no prefetch | 29.1 | 6.0 | 9,700 |
| 92% | 2 GB/s | no alternation | 29.2 | 7.7 | 1,300 |
| 92% | 500 MB/s | dev | 36.3 | 10.7 | 14,000 |
| 92% | 500 MB/s | streaming | 31.8 | 7.4 | 3,500 |
| 92% | 500 MB/s | no prefetch | 33.5 | 6.0 | 9,900 |
| 92% | 500 MB/s | no alternation | 37.3 | 12.3 | 6,400 |
| 92% | 150 MB/s | streaming | 66.7 | 7.2 | 700 |
| 92% | 150 MB/s | no prefetch | 68.3 | 6.8 | 9,200 |

Data about 3 times the memory (a third of the file cacheable: the peak
memory of the fit plus 477 MB; chunks of 32 MB, the default for that
memory; `run_ooc_third.sh`, `results/third_of_memory*.csv`). Two of the 22
single runs were disturbed (1.5-2 times slower than their repeats: the
machine is a shared VM), and are replaced by the medians of their repeats
(*):

| data | disk | dev | streaming | no prefetch | no alternation |
|---|---|---|---|---|---|
| 10,000 x 150,000 (in memory: 25 s) | 2 GB/s | 523 s, 295 GB | 33.1 s*, 12.6 GB | 36.7 s, 12.0 GB | 32.5 s*, 15.7 GB |
| | 500 MB/s | 619 s, 295 GB | 35.2 s, 13.2 GB | 39.7 s, 11.9 GB | 34.7 s, 15.5 GB |
| | 150 MB/s | (about 35 min) | 89.7 s, 11.8 GB | 89.9 s, 12.0 GB | 103.4 s, 14.8 GB |
| 1,000 x 1,500,000 (in memory: 30 s) | 2 GB/s | 34.2 s, 15.1 GB | 38.2 s*, 16.5 GB | 34.1 s, 11.6 GB | 37.9 s, 19.6 GB |
| | 500 MB/s | 43.7 s, 15.6 GB | 49.6 s, 16.1 GB | 42.5 s, 11.7 GB | 54.5 s, 19.7 GB |
| | 150 MB/s | 116.1 s, 15.1 GB | 122.0 s, 15.3 GB | 96.3 s, 11.6 GB | 145.4 s, 19.5 GB |

- 10,000 individuals: streaming takes 1.3-1.4 times the time in memory on
  the faster disks, dev 21-25 times; on the slowest disk, every version is
  bound by its reads (12 GB at 150 MB/s is 80 s).
- 1,000 individuals: dev stays within 1.1-1.5 times the time in memory on
  the faster disks. Streaming without the prefetch is as fast, and faster
  on the slowest disk (96 against 116 s, 23% fewer bytes thanks to the
  alternation); with the prefetch it reads 4 GB more, which costs 8 s at
  500 MB/s and 26 s at 150 MB/s: the chunks of its Q side are whole
  segments (187 MB here, against 32 MB chunks elsewhere), and the segment
  computed plus the one prefetched (374 MB) do not leave the 477 MB cache
  room to keep the first while all its blocks of individuals are processed.
  The prefetch of the Q side is now capped at `chunk_bytes` (next
  section).
- The alternation saves 15-40% of the bytes, which pays on the slowest
  disk (10,000 individuals: 90 against 103 s; 1,000: 96 against 145 s).
- System CPU time was 8-9 s with the prefetch (with or without the
  alternation) and 19 s without it (ten times as many page faults that
  waited for the disk), out of 120 s of CPU in all.

(The dev rows at 67% come from the first run of the matrix, with the same
build of dev; the streaming rows of that first run, before the prefetch was
split into 2 MB pieces and the chunks read forwards, are not shown.)

Reading the results:

- Many individuals (one Q segment): chunks are most of the gain. Each pass
  reads the file at most once, so streaming fits take 1.1-1.2 times the
  time in memory (25 s) where dev took 20-25 times as long (513-618 s,
  294 GB read for a 1.5 GB file).
- Few individuals (segments): dev reads each segment from the cache while
  its blocks are processed, so it is not pathological (1.7 times the bytes
  of streaming); streaming is 5-12% faster on the faster disks.
- The prefetch keeps the threads from waiting on page faults (with one Q
  segment, 147,000-929,000 waiting faults without it, 4,000-12,000 with
  it): 25-50% of the time on the faster disks. With the chunk order
  reversed, the kernel stops reading ahead, so the alternation needs the
  prefetch. It reads 5-20% more bytes (the chunk prefetched while one is
  computed evicts some of what the next pass would reuse), which costs
  that much when the disk is the bottleneck (150 MB/s, 67% cacheable:
  67.7 s against 65.0 s, and 74.2 against 68.8 s); at 92% it still won
  (50.2 against 56.0 s).
- Alternating directions halves the bytes read (15 to 8 GB at 67%, 16 to
  5 GB at 92%, with one segment), 5-15% of the time on these disks; more
  on a slower one, where bytes are time.
- On the slowest disk, every version is bound by the bytes read (8 GB at
  150 MB/s is 55 s): only reading less helps (smaller chunks; see "Left to
  do").

### The prefetch of the Q side

With segments, the Q side used to prefetch the whole next chunk (one or
more segments), which with a small cache evicted the segment being
computed before all its blocks of individuals were done. The prefetch now
reads at most `chunk_bytes` of the next chunk (`max_prefetch` in
`q_side_stream()`); the rest is read by the kernel's read-ahead of the
threads' forward page faults. 1,000 x 1,500,000, full fit as above:

| cache | disk | capped | uncapped | no prefetch | dev |
|---|---|---|---|---|---|
| 33% | 2 GB/s | 33.1, 34.9 s; 14.0 GB | 34.7, 38.2 s; 15.8-16.8 GB | 34.1-36.6 s; 11.6-12.2 GB | 34.2-35.8 s; 15.1-16.2 GB |
| 33% | 500 MB/s | 45.2 s; 12.8 GB | 49.6 s; 16.1 GB | 42.5 s; 11.7 GB | 43.7 s; 15.6 GB |
| 33% | 150 MB/s | 103.5 s; 12.6 GB | 122.0 s; 15.3 GB | 96.3 s; 11.6 GB | 116.1 s; 15.1 GB |
| 67% | 2 GB/s | 34.8, 34.7 s; 8.5 GB | 35.0, 34.1 s; 9.3 GB | 31.3 s; 8.0 GB | 31.9 s; 13.8 GB |

(The runs at 2 GB/s alternated the capped and uncapped builds; the others
are single runs, compared with those of the matrix above.) The cap removes
most of the extra reads, which saves 9-15% of the time on the slower disks
and costs nothing on the fast one. Streaming still reads 1 GB more than
without the prefetch here, and on the slowest disk without the prefetch is
7% faster; with 10,000 individuals (one segment, unchanged by the cap) the
prefetch pays (see above). Whether few individuals should prefetch at all
is for the desktop.

## Conventional fits

Data that fit in memory never stream (0.8 times the memory available), and
the usual layout is unchanged except for three arguments of
`accum_q_block()` (`init`, `out`, `ll0`: true, true, 0). Streaming off
against dev, with the file in the page cache, 4 threads, K = 4, processes
alternated dev, new, new, dev, dev, new, identical estimates:

| data | dev (s) | new (s) |
|---|---|---|
| fit, 10,000 x 150,000 (3 iterations, 2 EM steps) | 40.1, 38.6, 37.8 | 41.6, 40.1, 37.1 |
| fit, 1,000 x 1,500,000 | 42.6, 37.0, 34.5 | 38.2, 35.8, 35.5 |

and each pass (`bench_kernels.R`, the median of 3 runs in each process;
the medians of the 3 processes):

| pass | 10,000 x 150,000: dev, new (s) | 1,000 x 1,500,000: dev, new (s) |
|---|---|---|
| `update_q` | 2.85, 3.04 | 2.70, 2.61 |
| `update_q_loglik` | 3.69, 3.76 | 3.32, 3.12 |
| `update_p` | 2.16, 2.19 | 2.20, 2.20 |
| `loglik` | 2.10, 1.94 | 2.19, 2.00 |
| `em` | 4.06, 3.70 | 3.71, 3.46 |

The differences go both ways and are within the variation between runs of
the same build on this machine (a shared VM: up to 15-20%), so no slowdown
could be measured; confirm on the desktop before merging.

Streaming forced on in memory: on 10,000 x 150,000 the Q side is 15-20%
faster (all the threads read the same columns, which stay in the cache),
and on 1,000 x 1,500,000 the passes are within 5% (9 runs each). The
variation between runs of the same build on this machine (a shared VM) was
up to 15-20%.

## Desktop (manica01, 24-25 November 2026)

i9-13900, 8 P-cores (`taskset -c 0-7`, `OMP_PROC_BIND=close`), 125 GB,
NVMe Samsung 970 EVO Plus (`read_ahead_kb` and `max_sectors_kb` 128),
Linux 6.8 with cgroup v2. Builds at -O2 with the kit of the tag
`archive/newspeed` (`build_variants.sh`, AVX2 clones checked): base =
51b4e9f (`ada196_article` before the merge of streaming), new = 001667d
(the merge, with the blocked Hessians), new2 = b07f087 (the prefetch in
pieces of 128 KB, below). The same data and fit as in the container
(`make_data.R`; `bench_fit.R`: K = 4, one EM step, 2 iterations), on 8
cores. The hard disks of the desktop are not mounted (and would be outside
the folders exempt from the antivirus), so the slow disks are the NVMe
throttled to 500 and 150 MB/s (`IOReadBandwidthMax` of a system scope,
with sudo); the memory caps are scopes of the user's systemd. Scripts:
`data-raw/radmix/stream/desktop/` (the driver and its phases); results:
`data-raw/radmix/stream/results/desktop_*`. Every run gave the same
estimates as in memory (one hash per data set, 72 runs).

### In memory: no slowdown, except small data on one thread

The kit's `run_bench.sh` (full fits to convergence, 3 balanced rounds) and
`MODE=kernels`, base against new, on eur (379 x 54,051) and hgdp
(1,043 x 130,943); median seconds:

| data | K | 1 thread: base, new | 8 threads: base, new |
|---|---|---|---|
| eur | 3 | 4.86, 5.06 (+4.1%) | 0.66, 0.67 |
| eur | 4 | 6.07, 6.26 (+3.2%) | 0.81, 0.82 |
| hgdp | 5 | 41.61, 40.95 (-1.6%) | 5.09, 5.14 |
| hgdp | 10 | 132.75, 133.36 (+0.5%) | 17.75, 17.93 |

Identical log-likelihood histories, and 1 and 8 threads identical. The
eur fits on one thread are 3-4% slower in every round (ranges 4.86-4.88
against 5.05-5.07 s, and 6.02-6.14 against 6.26-6.33 s); the passes put
it in `update_p` and `loglik` (+5-7% on eur, 1 thread; not on hgdp). The
cause is not known yet (open, below).

The fit of the simulated data in memory (8 cores, warm page cache, 3 runs
each, alternated; `desktop_in_memory.csv`):

| data | base | new, not streaming | new, streaming forced |
|---|---|---|---|
| 10,000 x 150,000 | 6.59 s (6.51-6.62) | 6.56 s (6.50-6.63) | 5.61 s (5.54-5.88) |
| 1,000 x 1,500,000 | 5.98 s (5.94-6.03) | 5.97 s (5.97-6.11) | 6.17 s (6.12-6.25) |

### The prefetch read only 1/16 of each chunk

On the desktop, the prefetch first made no difference: as many major page
faults with it as without (156,000-226,000 per fit with 10,000
individuals). Each `madvise(MADV_WILLNEED)` call makes Linux read at most
the larger of the disk's read-ahead and its largest request
(`force_page_cache_ra()`: `max(io_pages, ra_pages)`), 4 MB in the
container but 128 KB here (the kernel's default read-ahead), so the pieces
of 2 MB read 128 KB each: 4,096 of the 65,536 pages of a 256 MB range,
against all of them within 50 ms in pieces of 128 KB or 64 KB (`mincore()`
on manica01). b07f087 calls it in pieces of 128 KB, which every Linux
block device honours (128 KB is the default read-ahead; NVMe, SATA,
virtio and USB disks all allow requests at least that large). Alternated
runs, unthrottled (`desktop_prefetch_fix.csv`): 10,000 x 150,000 at 92%
cacheable, 9.7-10.1 s and 156,000-173,000 faults with pieces of 2 MB,
7.4-7.5 s and 13,000-20,000 with pieces of 128 KB, 10.1-10.2 s without
the prefetch; 1,000 x 1,500,000 at 67%, 8.6-8.9 s for all three.

### Out of core

Seconds per fit (bytes read; GB here are 2^30 bytes), new2 for streaming (`desktop_third_pf128.csv`,
`desktop_ooc_pf128.csv`), base from the first run (`desktop_third.csv`,
`desktop_ooc.csv`); chunks of 32 MB at 33% cacheable, 64 MB otherwise.
Base was not run at 150 MB/s with 10,000 individuals (about 35 minutes
each), nor at 67-92% at 150 MB/s.

10,000 x 150,000 (in memory: 6.6 s):

| cache | disk | base | streaming | no prefetch | no alternation |
|---|---|---|---|---|---|
| 33% | NVMe | 105.5 (148 GB) | 9.1 (12.2 GB) | 11.4 (11.0 GB) | 10.5 (16.9 GB) |
| 33% | 500 MB/s | 318.0 (148 GB) | 25.2 (11.3 GB) | 24.9 (11.0 GB) | 33.5 (15.6 GB) |
| 33% | 150 MB/s | | 81.4 (11.3 GB) | 79.7 (11.0 GB) | 110.2 (15.4 GB) |
| 67% | NVMe | 104.2 (148 GB) | 7.9 (9.5 GB) | 10.4 (8.5 GB) | 9.6 (15.5 GB) |
| 67% | 500 MB/s | 317.9 (148 GB) | 20.3 (8.6 GB) | 20.1 (8.5 GB) | 32.5 (15.1 GB) |
| 67% | 150 MB/s | | 64.1 (8.7 GB) | 62.1 (8.5 GB) | |
| 92% | NVMe | 23.9 (16.8 GB) | 7.2 (6.7 GB) | 9.8 (5.4 GB) | 8.8 (15.1 GB) |
| 92% | 500 MB/s | 37.5 (16.8 GB) | 15.6 (6.3 GB) | 17.2 (7.0 GB) | 32.1 (15.0 GB) |
| 92% | 150 MB/s | | 46.0 (6.2 GB) | 53.4 (7.2 GB) | |

1,000 x 1,500,000 (in memory: 6.0 s):

| cache | disk | base | streaming | no prefetch | no alternation |
|---|---|---|---|---|---|
| 33% | NVMe | 9.4 (15.4 GB) | 9.7 (11.9 GB) | 9.0 (10.6 GB) | 10.7 (16.6 GB) |
| 33% | 500 MB/s | 33.3 (15.4 GB) | 26.5 (11.6 GB) | 24.4 (10.6 GB) | 35.6 (16.4 GB) |
| 33% | 150 MB/s | 110.3 (15.4 GB) | 83.3 (11.5 GB) | 77.7 (10.6 GB) | 116.8 (16.3 GB) |
| 67% | NVMe | 9.3 (14.1 GB) | 8.6 (9.0 GB) | 8.4 (8.5 GB) | 10.1 (15.5 GB) |
| 67% | 500 MB/s | 30.8 (14.0 GB) | 21.3 (9.0 GB) | 20.6 (8.7 GB) | 33.9 (15.7 GB) |
| 67% | 150 MB/s | | 63.1 (8.6 GB) | 62.5 (8.5 GB) | |
| 92% | NVMe | 8.5 (5.6 GB) | 7.5 (4.7 GB) | 8.2 (5.9 GB) | 8.2 (7.6 GB) |
| 92% | 500 MB/s | 15.1 (5.3 GB) | 14.4 (5.1 GB) | 14.0 (4.6 GB) | 18.4 (7.6 GB) |
| 92% | 150 MB/s | | 37.4 (4.8 GB) | 40.7 (5.3 GB) | |

- On the NVMe, streaming takes 1.1-1.4 times its time in memory with
  10,000 individuals (base 3.6-16 times: it reads the 1.5 GB file about
  100 times at 33-67%), and 1.25-1.6 times with 1,000 individuals, where
  base is as fast (9.4 against 9.7 s at 33%) and slower only on the
  throttled disks.
- On the throttled disks, every fit is bound by its reads (11.3 GB at
  150 MB/s is 81 s): the time follows the bytes, which fall as more of
  the file is cached (81, 64, 46 s at 33, 67, 92%, 10,000 individuals).
- The prefetch (fixed) saves 20-27% on the NVMe with 10,000 individuals,
  where the threads otherwise wait on page faults, 8-14% at 92% cacheable
  on the slowest disk, and 8-9% at 92% on the NVMe (1,000 individuals) and
  at 500 MB/s (10,000); it costs 7-9% at 33% with 1,000 individuals, on
  every disk (it reads 10-12% more), and is otherwise within 3.5%.
- The alternation saves 28-56% of the bytes, and 22-51% of the time on the
  throttled disks.
- The per-pass picture (`desktop_passes.csv`, 33% cacheable, the build
  with pieces of 2 MB): with 10,000 individuals, base's Q side reads the
  file 20 times per pass (28.6 GB, 19.8 s, against 0.70 s in memory);
  streaming's, once (1.5 s).

## The article's runs (manica01, 25-26 November 2026)

The benchmark article's section "Data larger than the memory": the
article's own scripts (`vignettes/articles/code/admixture_benchmark_large*`),
`gt_radmix()` with the default options, K = 5, 8 P-cores, tidygenclust
0.1.5.9007; caps by `systemd-run`, throttles as above. Results in
`vignettes/articles/data/admixture_benchmark/{mid,large}_*.csv`.

| data | cap | disk | time (s) | times in memory | read (GB) |
|---|---|---|---|---|---|
| 50,000 x 100,000 (5 GB) | none | | 214 (was 255) | 1.0 | 0 |
| | 4 GiB | NVMe | 277 | 1.29 (was 1.08) | 238 |
| | 1 GiB | NVMe | 332 | 1.55 (was 1.30) | 476 |
| | 4 GiB | 500 MB/s | 547 | 2.6 | 225 |
| | 1 GiB | 500 MB/s | 917 | 4.3 | 444 |
| | 4 GiB | 150 MB/s | 1,648 | 7.7 | 237 |
| | 1 GiB | 150 MB/s | 2,972 | 13.9 | 441 |
| 250,000 x 200,000 (50 GB) | none | | 1,711 (was 2,173) | 1.0 | 0 |
| | 8 GiB | NVMe | 2,290 | 1.34 (was 1.05) | 3,579 |
| | 8 GiB | 500 MB/s | 7,330 | 4.3 | 3,532 |

The fits in memory were rerun on 27 November with 79e21c4 (chunks for
the Q side in memory, see "Chunks in memory"; `phases/inmem_cq.sh`, the
article's script with `LARGE_READ_MBPS=" "`), the "was" values are those
of 25-26 November without them. The streamed runs always took the Q side
in chunks, so the first comparison understated the cost of the disk: the
NVMe runs gained the chunks' 15-20% and lost it again to the reads. A
second 5 GB fit in memory took 221.0 s (214.4 s in the article): 3%
between two runs.

Identical log-likelihoods within each data set (22 and 17 iterations);
the capped runs peaked at their cap (R and the page cache), R itself at
0.21 and 0.40 GB. Throttled, the time is that of the reads (3,532 GB at
500 MB/s is 7,064 s). On the NVMe, the 50 GB fit read 1.56 GB/s on
average, against the 2.1 GB/s at which the same passes go through the
data in memory, hence 34%; with 5 GB under 1 GiB, 1.43 GB/s against
2.2 GB/s, hence 55%.

The one-thread regression of the gate (eur, 3-4%) did not reproduce in
5 rounds (base 4.89, merge 4.88, fix 4.88 s at K = 3). The merge had
stopped GCC from inlining `loglik_slice()` into `sum_slices()` (a second
copy of its loop for streaming); 524b347 keeps one copy, and the
log-likelihood pass is 7-11% faster than before the merge (eur K = 3:
0.054, 0.058, 0.050 s; hgdp K = 5: 0.419, 0.414, 0.373 s). The 8-thread
HGDP fits of the merge and of the article's build are equal to base
(5.07-5.14 against 5.09 s); on one thread they were 1.6% faster and 3%
slower in two sessions: code placement, within noise.

## Chunks in memory (manica01, 26-27 November 2026)

Andrea, 26 Nov: keep the defaults (stream above 0.8 x the memory
available, prefetch in pieces of 128 KB, alternating directions), and
take the Q side in chunks for data in memory too. c521e36: `chunk_q`, set
by `radmix_stream_mode()` when the data are not streamed (unless the
option `tidygenclust.radmix_stream` is `FALSE`, and not for sources of
genotype likelihoods), sends a Q side with one segment (8,192 individuals
or more) through `q_side_stream()`, without the prefetch. The results are
identical (tests in `test_radmix_stream.R`, chunks of one byte).

The chunk size, `phases/chunkq_bench.sh` (build c521e36, `bench_fit.R`,
K = 5, one EM step and 2 iterations, 1 for the 50 GB data; 8 P-cores,
warm page cache, alternated rounds; 23:53-00:21, nothing else running).
Medians (min-max) in seconds, `chunkq.csv` in the work folder:

| data | usual | 8 MB | 16 MB | 32 MB | 64 MB |
|---|---|---|---|---|---|
| 10,000 x 150,000 (3 rounds) | 7.39 (7.27-7.40) | 6.27 (6.10-6.31) | 6.26 (6.17-6.35) | 6.40 (6.34-6.49) | 6.37 (6.36-6.46) |
| 50,000 x 100,000 (3) | 25.42 (25.03-25.98) | 21.25 (21.04-21.80) | 21.62 (20.95-21.65) | 21.06 (21.05-21.17) | 21.13 (21.10-21.54) |
| 250,000 x 200,000 (2) | 163.1 (162.4-163.8) | | 135.3 (135.3-135.4) | 133.0 (132.7-133.3) | 132.6 (132.5-132.7) |

Chunks are 15-19% faster than the usual layout, with identical hashes
and log-likelihoods; the sizes are within about 2% of each other. 64 MB
(`radmix_chunk_q_mb()`) is best on the largest data, where 16 MB is a
single tile of 64 loci, and 1.7% behind 16 MB on the smallest. Not
tested: larger K, where a chunk's P (loci x K doubles, up to 1.25 MB at
K = 20 with 8,192 individuals and 64 MB) may outgrow the L2 cache. The
streamed runs of the article are unchanged (chunk_q is off when
streaming); its in-memory baselines are rerun with this build.

## Left to do

- The SVD start of genotype likelihoods (`init = "svd"`) has the same
  pattern in its product X b (`radmix_gl_dosage_prod_cpp()`, one task per
  slice of individuals over all the loci): it would read the file many
  times per product out of core. The same chunked loop, with the
  accumulators kept between chunks, would fix it. The SVD start of
  genotypes is `tidypopgen::gt_pca_partialSVD()`, whose memory and reads
  are bigstatsr's.
- On Windows `prefetch()` does nothing (`PrefetchVirtualMemory()` would
  do), and `radmix_mem_available()` is `Inf` outside Linux, so streaming is
  never chosen on its own there.
- Memory that still grows with the data: P (J x K doubles, about 9 copies
  in the fit) and, in streaming mode with 8,192 individuals or more, the
  accumulators kept between chunks (I (K + K (K + 1) / 2) doubles).
- The chunk size of `chunk_q` at large K (10, 20): 16 against 64 MB on
  10,000 x 150,000 (see "Chunks in memory").
- Timing on a real hard disk and on a network file system (NFS, Lustre):
  the desktop's NVMe and the container's disk are both fast, and the
  throttle caps the throughput, not the latency. Pieces of 128 KB may make
  many small requests on a network file system; if they cost there, read
  the limit of the disk from `/sys` and use larger pieces where it allows.
  A 10-second check that the prefetch covers its range on a machine:
  `mincore()` after `madvise()` (the Python script of 25 November, to be
  added to `data-raw/radmix/stream/`).
- Few individuals, file a little larger than the cache: on the desktop,
  streaming is within 3% of dev on the NVMe at 33% cacheable (9.7 against
  9.4 s) and faster everywhere else, so streaming need not be limited to
  one Q segment. The prefetch costs 7-9% there (it reads 10-12% more);
  whether to prefetch less in that case is open.

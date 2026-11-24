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
  runs of adjacent columns, in pieces of 2 MB. Linux reads only the first
  4 MB of a larger range (a 256 MB call left 1,024 of 65,536 pages
  resident), whereas each 2 MB call starts its reading and returns at once
  (256 MB advised in 17 ms, all resident 0.5 s later). On cached pages it
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
- Streaming in memory made the Q side of many individuals faster: chunking
  the Q side of one segment could become the default for data that fit in
  memory too, after timing it on the desktop.
- Timing on a real disk (the container's disk is backed by the host's
  cache, and the throttle only caps its throughput, not its latency).
- Few individuals, file a little larger than the cache: the chunks and the
  prefetch can make streaming read more than dev (step 1), and dev is
  already close to its time in memory there. Whether streaming should
  apply only to one Q segment (the case that needs it), or prefetch less
  (smaller chunks, or not past what the cache can hold), is for the
  desktop.

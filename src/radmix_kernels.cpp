// The passes over the genotypes of gt_radmix(): the log-likelihood, the
// deviance of held-out genotypes, the EM step, and the gradients, Hessians
// and Newton updates of the two blocks of parameters. See radmix_kernels.h
// for the model, the notation and the layout of q and p, and
// dev_notes/admixture_implementation.md for the algorithm as a whole.
//
// Origin of the code: the first version of these passes was ported from
// src/algorithms_inner.jl and src/loops.jl of OpenADMIXTURE.jl (MIT
// licence, see inst/COPYRIGHTS). They have since been rewritten for
// tidygenclust: the derivatives (OpenADMIXTURE.jl differentiates a
// homogeneous form of the log-likelihood with respect to q), the data
// access, the layout of the work, the accumulation, the log-likelihood,
// the deviance and genotype likelihoods below are new. The quadratic
// programs and the projections that they call are still translations of
// OpenADMIXTURE.jl (radmix_qp_sweep.cpp, radmix_projections.cpp).
//
// Per-genotype quantities
// -----------------------
// For individual i and locus j, write g = g_ij (copies of the alternate
// allele), h = n_i - g (copies of the reference allele, n_i = 2 for diploids,
// 1 for pseudohaploids), and qp = qp_ij = sum_k q_ki p_kj. A missing (or
// hidden) genotype has g = h = 0, so it drops out of every sum below without
// a test. The genotype contributes l(qp) = g log(qp) + h log(1 - qp) to the
// log-likelihood, and every pass is built from the following coefficients:
//   a  = g / qp,  b = h / (1 - qp)            (EM step)
//   d1 = l'(qp)  = a - b                      (gradients)
//   d2 = -l''(qp) = g / qp^2 + h / (1 - qp)^2 (Hessians; d2 >= 0)
// With genotype likelihoods (from another package), g and h are the
// expected counts of allele2 and allele1 given the likelihoods and qp, and
// the log-likelihood is that of the likelihoods (see slice_gl() and
// loglik_gl_slice()); every pass is otherwise the same.
//
// Gradients and Hessians, by the chain rule
// -----------------------------------------
// qp_ij is linear in q_.i (with coefficients p_.j) and in p_.j (with
// coefficients q_.i), so
//   dL / dq_ki  = sum_j d1 p_kj,   -d2L / dq_ki dq_k2i = sum_j d2 p_kj p_k2j
//   dL / dp_kj  = sum_i d1 q_ki,   -d2L / dp_kj dp_k2j = sum_i d2 q_ki q_k2i
// The two blocks have the same structure: Q-side sums run over the loci of
// one individual, weighted by the columns of p; P-side sums over the
// individuals of one locus, weighted by the columns of q. The negative
// Hessians are sums of d2 v v' with d2 >= 0, hence positive semi-definite
// (L is concave in each column). Given p, L is a sum of one term per
// individual, and given q one term per locus, so each column of q or p has
// its own K x K Hessian, and is updated on its own.
//
// Layout of the work
// ------------------
// The genotypes are stored individuals x loci, column-major, so the
// genotypes of consecutive individuals at one locus are adjacent bytes. The
// unit of work is a slice: one locus and up to TILE_I consecutive
// individuals. For each slice, qp, the decoded genotypes and the
// coefficients are computed into small thread-local buffers (Slice), with
// loops over the individuals of the slice that have no branches and
// vectorise (except the table lookups that decode the genotypes, see
// slice_geno()). The sums are then accumulated:
// - Q side (one sum per individual): a block of TILE_I individuals is
//   streamed over the loci, with one accumulator vector of length TILE_I
//   per gradient entry and per entry of the Hessians, each updated by a
//   unit-stride loop over the individuals of the slice;
// - P side (one sum per locus): for a block of TILE_J loci, each sum is a
//   dot product over the individuals of the slice, and the TILE_I x K corner
//   of q stays in cache while the loci of the block are processed.
// Loci are processed two at a time: the coefficients of both loci are
// computed first, and the accumulators are then read and written once for
// both, which halves their memory traffic (the accumulation is limited by
// loads and stores rather than arithmetic). A range with an odd number of
// loci is padded with a locus whose coefficients are zero.
// Only the K (K + 1) / 2 entries of each symmetric Hessian with k2 >= k are
// accumulated; the others are filled in at the end (mirror_hess()).
//
// Parallelism and determinism
// ---------------------------
// Q-side work is parallelised over blocks of individuals and P-side work
// over blocks of loci, so each thread writes its own output columns: there
// are no races and no reductions across threads. With few individuals the
// Q side also splits the loci into segments (q_segments()), whose partial
// sums are added up in a fixed order. The scalar sums (log-likelihood,
// deviance) keep one partial sum per block of loci, added up in block order.
// Block and segment boundaries depend only on the dimensions of the data,
// never on the number of threads, and each block is processed start to
// finish by one thread, so the results are bit-identical for any n_cores.
// (`omp simd reduction` regroups a sum, but in a way fixed at compile time.)
//
// Thread safety: inside the parallel regions, the threads only read
// immutable data (the genotypes, q, p and the transposed copy of q) and
// write thread-local buffers or their own output columns. They never call
// R: R is not thread safe, and an R error would jump out of the parallel
// region. All the checks happen before, in the R wrappers
// (radmix_exports.cpp, radmix_fit.cpp) and the GenoView constructor.
//
// Streaming
// ---------
// With a view in streaming mode (data larger than the memory, see
// radmix_geno_access.h), each pass walks the loci in chunks (of about
// chunk_bytes bytes of the backing file, in whole tiles), one parallel
// region per chunk, so that all the threads read the same part of the file,
// and a helper thread (Prefetcher) asks the kernel to read the next chunk
// while they compute. The blocks and segments, and the order in which each
// sum is accumulated, are those of the usual layout, so the results are
// bit-identical:
// - Q side, one segment: each block of individuals is accumulated over the
//   loci of a chunk, then of the next, from where it stopped: its
//   accumulators (and its log-likelihood) are kept between chunks, in
//   I (K + K (K + 1) / 2) doubles (accum_q_block() with init and out
//   false), and chunks have an even number of loci, so the loci are paired
//   as in one run;
// - Q side, several segments: a chunk is a group of whole segments, whose
//   tasks are those of the usual layout (it can be larger than chunk_bytes,
//   so only its first chunk_bytes are prefetched);
// - P side, log-likelihood, deviance: each block of loci is independent (or
//   keeps its own partial sums), so the chunks can be processed in any
//   order: in order or in reverse order, as GenoView::direction() says
//   (the loci of a chunk are always read forwards).
// Streaming is never used for data that fit in memory, so the usual layout
// keeps its speed.

#include <algorithm>
#include <cmath>
#include <thread>
#include <utility>
#include <vector>
#include "radmix_geno_access.h"
#include "radmix_kernels.h"
#include "radmix_projections.h"
#include "radmix_qp_sweep.h"

namespace radmix {

// Function multiversioning. R compiles packages for the baseline x86-64
// instruction set (SSE2: 2 doubles per instruction, no fused multiply-add),
// which leaves most of the speed of the `omp simd` loops unused. On x86-64
// Linux with GCC >= 12, the two accumulation functions, which hold the hot
// loops, are also compiled for the x86-64-v3 level (AVX2: 4 doubles per
// instruction, and FMA), and the version that suits the processor is chosen
// once, when the package is loaded (an ifunc). This makes the Newton updates
// and EM steps 15-30% faster. The log-likelihood and the deviance passes
// are left alone (their AVX2 versions were slower, with log() called for
// some genotypes inside their loops); the log-likelihood computed by the
// fused pass (update_q_loglik()) is inlined into the clone of
// accum_q_block(). The two versions can differ by rounding (FMA and wider
// vectors change the order of operations), so results are
// reproducible on a given machine but not bit-identical across machines;
// defining RADMIX_NO_MULTIVERSION when compiling (e.g. in ~/.R/Makevars)
// builds only the baseline version, whose results are the same on every
// x86-64 machine.
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 12 && \
    defined(__x86_64__) && defined(__linux__) && \
    !defined(RADMIX_NO_MULTIVERSION)
#define RADMIX_MULTIVERSION \
  __attribute__((target_clones("arch=x86-64-v3", "default")))
#else
#define RADMIX_MULTIVERSION
#endif

// A function called from a multiversioned one is compiled for AVX2 only if
// it is inlined into it, and GCC does not inline the larger ones on its own:
// slice_geno(), which computes qp and decodes the genotypes of every slice,
// remained a call to its baseline version from every AVX2 clone (and so did
// slice_coef(), which calls it, once slice_geno() was inlined into it).
// Such functions are therefore forced inline.
#if defined(__GNUC__)
#define RADMIX_INLINE inline __attribute__((always_inline))
#else
#define RADMIX_INLINE inline
#endif

// Tile sizes. A vector of TILE_I accumulators takes 512 bytes, so the
// Q-side accumulators (K + K (K + 1) / 2 vectors) and the TILE_I x K corner
// of q stay in the L1 or L2 cache for moderate K. loglik_slice() needs
// TILE_I == LOG_BATCH^2, and the padding of odd ranges of loci on the P side
// needs TILE_J to be even (see accum_p_block()).
//
// TILE_J, MIN_Q_TASKS (below) and the OpenMP schedule of the passes can be
// changed at compile time, for benchmarks only (e.g. -DRADMIX_TILE_J=32
// -DRADMIX_MIN_Q_TASKS=64 -DRADMIX_OMP_SCHEDULE=guided in the CXXFLAGS of a
// private Makevars). On the desktop (i9-13900, 8 threads, HGDP K = 5 and
// 10, 15 Nov 2026), MIN_Q_TASKS = 128 made fits 3-6% faster than 32, and
// 2-3% faster on 1 thread; TILE_J = 32 changed nothing, and a guided
// schedule was 5-7% slower than the dynamic one. The schedule never
// changes the results (each block is processed start to finish by one
// thread). TILE_J changes the partial sums
// of loglik() and deviance(), and MIN_Q_TASKS the segments of the Q side, so
// with other values results differ by rounding (and test_radmix_kernels.R,
// which checks the segments of the default, fails), but they are still the
// same for any number of threads.
#ifndef RADMIX_TILE_J
#define RADMIX_TILE_J 64
#endif
#ifndef RADMIX_MIN_Q_TASKS
#define RADMIX_MIN_Q_TASKS 128
#endif
#ifndef RADMIX_OMP_SCHEDULE
#define RADMIX_OMP_SCHEDULE dynamic
#endif
#ifndef RADMIX_HESS_BLOCK_MIN_K
#define RADMIX_HESS_BLOCK_MIN_K 7
#endif
static const size_t TILE_I = 64;
static const size_t TILE_J = RADMIX_TILE_J;

// Register blocking of the Hessians, from K = HESS_BLOCK_MIN_K. On the Q
// side, the
// sums hT[i, (k, k2)] += d2_ij p_kj p_k2j over the loci (K (K + 1) / 2 per
// individual, against K for the gradients: 54-56% of a pass at K = 10)
// keep their accumulators in memory, so each step of the loop below loads
// three values and stores one for two multiply-adds. Blocked, the d2 of
// HB_L loci (HB_L x TILE_I) and the products of their frequencies (HB_L x
// n_pair, padded to a multiple of HB_P pairs) are buffered first; then, for
// each group of HB_P pairs and of 8 individuals, HB_P x 8 accumulators
// (eight AVX2 registers) are summed over the HB_L loci before being added
// to hT, so each multiply-add loads about a quarter of a value. It is the
// same loops in another order (unroll and jam), in portable C++ (no
// intrinsics: the compiler vectorises the innermost loops for the baseline
// and the AVX2 clones alike), and it adds the same products in another
// order, so results change by rounding only (still identical for any
// number of threads). The P side is blocked in the same way, with the roles
// of the individuals and the loci swapped (see accum_p_block()). On the
// desktop (8 threads, HGDP), the Q side took 0.86-0.94 of the time per
// iteration at K = 8-12, and the P side then 0.88-0.93 of the rest, with
// identical iterations and estimates within 3e-12; K = 7 gained 6%, K = 6
// and below nothing, so they keep the plain loops. (In single precision,
// the blocked sums were no faster per iteration than in double precision,
// and changed the paths of the fits; blocks of 6 pairs, 12 accumulators,
// were 1.3-1.4 times slower than blocks of 4; adding the gradients to the
// Q-side blocks gained nothing.)
static const int HESS_BLOCK_MIN_K = RADMIX_HESS_BLOCK_MIN_K;
static const int HB_L = 16;
static const int HB_P = 4;
static const int HB_V = 8;
// the number of pairs (k, k2) of the buffers, padded to a multiple of HB_P
static inline size_t hess_pairs_pad(size_t n_pair) {
  return (n_pair + HB_P - 1) / HB_P * HB_P;
}

// The Hessian sums of nl buffered loci (see above): d2b (nl x TILE_I, zero
// beyond the first ni individuals) and ppb (nl x npair_pad, zero beyond
// n_pair), added into hT (n_pair x TILE_I).
static RADMIX_INLINE void hess_block(const double* __restrict d2b,
                                     const double* __restrict ppb, int nl,
                                     size_t n_pair, size_t npair_pad,
                                     size_t ni, double* __restrict hT) {
  const size_t ni_pad = (ni + HB_V - 1) / HB_V * HB_V;
  for (size_t pr = 0; pr < n_pair; pr += HB_P) {
    const size_t npr = std::min((size_t) HB_P, n_pair - pr);
    for (size_t ii0 = 0; ii0 < ni_pad; ii0 += HB_V) {
      double acc[HB_P][HB_V];
#pragma GCC unroll 4
      for (int q = 0; q < HB_P; q++)
#pragma omp simd
        for (int v = 0; v < HB_V; v++) acc[q][v] = 0;
      for (int l = 0; l < nl; l++) {
        const double* x = d2b + (size_t) l * TILE_I + ii0;
        const double* c = ppb + (size_t) l * npair_pad + pr;
#pragma GCC unroll 4
        for (int q = 0; q < HB_P; q++) {
          const double cq = c[q];
#pragma omp simd
          for (int v = 0; v < HB_V; v++) acc[q][v] += x[v] * cq;
        }
      }
      for (size_t q = 0; q < npr; q++) {
        double* h = hT + TILE_I * (pr + q) + ii0;
#pragma omp simd
        for (int v = 0; v < HB_V; v++) h[v] += acc[q][v];
      }
    }
  }
}

// Buffer the d2 of one locus (ni values, then zeros) and the products of
// its frequencies (n_pair values, then zeros) as row l of d2b and ppb.
static RADMIX_INLINE void hess_buffer(const double* d2, const double* pj,
                                      int K, size_t ni, int l, double* d2b,
                                      double* ppb, size_t npair_pad) {
  double* x = d2b + (size_t) l * TILE_I;
  for (size_t ii = 0; ii < ni; ii++) x[ii] = d2[ii];
  for (size_t ii = ni; ii < TILE_I; ii++) x[ii] = 0;
  double* c = ppb + (size_t) l * npair_pad;
  size_t pr = 0;
  for (int k = 0; k < K; k++)
    for (int k2 = k; k2 < K; k2++) c[pr++] = pj[k] * pj[k2];
  for (; pr < npair_pad; pr++) c[pr] = 0;
}
static_assert(TILE_J % 2 == 0, "TILE_J must be even");
// (the blocked P-side Hessians store the d2 of a block of loci in rows of
// TILE_I values)
static_assert(TILE_J <= TILE_I, "TILE_J must not exceed TILE_I");
// (the bit mask of a training view gives one 64-bit word per slice, see
// slice_geno())
static_assert(TILE_I == 64, "slice_geno() requires TILE_I == 64");

// the number of blocks of size `tile` needed to cover n (the last one may be
// shorter; 0 for n = 0)
static inline size_t n_blocks(size_t n, size_t tile) {
  return (n + tile - 1) / tile;
}

// The helper thread of a streaming pass (see "Streaming" at the top):
// start(j0, j1) waits for the previous request, then asks the view to
// prefetch loci [j0, j1) on a new thread, and returns at once. Prefetching
// is only a hint, so a thread that cannot be created is ignored. Used on
// the main thread only, between the parallel regions of the chunks.
class Prefetcher {
public:
  explicit Prefetcher(const GenoView& geno) : _geno(geno) {}
  ~Prefetcher() { wait(); }
  void start(size_t j0, size_t j1) {
    wait();
    if (j0 >= j1) return;
    try {
      _thread = std::thread([this, j0, j1] { _geno.prefetch(j0, j1); });
    } catch (...) {
    }
  }
  void wait() {
    if (_thread.joinable()) _thread.join();
  }
private:
  const GenoView& _geno;
  std::thread _thread;
};

// The number of blocks of TILE_J loci (defined below) per chunk of a
// streaming pass
static size_t chunk_blocks(const GenoView& geno);

// Fill in a symmetric K x K matrix (column-major) whose entries H[k2 + K k]
// with k2 >= k (the lower triangle) have been computed, by copying each of
// them into H[k + K k2] (the upper triangle).
static inline void mirror_hess(double* H, int K) {
  for (int k = 0; k < K; k++)
    for (int k2 = k + 1; k2 < K; k2++) H[k + K * k2] = H[k2 + K * k];
}

// The transpose of q (K x I): an I x K column-major matrix, in which the
// values q_k. of consecutive individuals are contiguous, so that the loops
// over the individuals of a slice read q with unit stride (and SIMD loads).
// Every pass builds it once, before its parallel region (it takes K I
// doubles, the size of q, and costs much less than the pass itself).
static std::vector<double> transpose(const double* q, int K, size_t I) {
  std::vector<double> qT((size_t) K * I);
  for (size_t i = 0; i < I; i++)
    for (int k = 0; k < K; k++) qT[i + I * k] = q[k + K * i];
  return qT;
}

// The transpose of the columns idx[0], ..., idx[n - 1] of q: an n x K
// column-major matrix, whose row f is individual idx[f] (the "compact"
// transpose of the free individuals of the Q side, see q_side()).
static std::vector<double> transpose_subset(const double* q, int K,
                                            const size_t* idx, size_t n) {
  std::vector<double> qT((size_t) K * n);
  for (size_t f = 0; f < n; f++)
    for (int k = 0; k < K; k++) qT[f + n * k] = q[k + K * idx[f]];
  return qT;
}

// Thread-local buffers for one slice (one locus x up to TILE_I individuals):
// qp, the decoded genotypes g and h, and the coefficients defined at the top
// of the file, a and b (EM) or d1 and d2 (Newton). loglik_slice() uses a and
// b as scratch, and accum_p_block() uses tmp. Each thread allocates its
// slices once, at the start of a parallel region, and reuses them for every
// slice it processes; every buffer is recomputed for each slice. Entry ii of
// each buffer refers to individual i0 + ii of the slice; only the first
// ni = i1 - i0 entries are meaningful. With genotype likelihoods, e1, e2 and
// obs hold the decoded likelihoods of the slice, and rf = obs / f (see
// slice_gl()).
struct Slice {
  std::vector<double> qp, g, h, a, b, d1, d2, tmp, e1, e2, obs, rf;
  Slice() : qp(TILE_I), g(TILE_I), h(TILE_I), a(TILE_I), b(TILE_I),
    d1(TILE_I), d2(TILE_I), tmp(TILE_I), e1(TILE_I), e2(TILE_I),
    obs(TILE_I), rf(TILE_I) {}
  // zero coefficients for the first ni individuals: the padding locus of a
  // range with an odd number of loci, which adds nothing to any sum
  void zero_coef(size_t ni) {
    std::fill(a.begin(), a.begin() + ni, 0.0);
    std::fill(b.begin(), b.begin() + ni, 0.0);
    std::fill(d1.begin(), d1.begin() + ni, 0.0);
    std::fill(d2.begin(), d2.begin() + ni, 0.0);
  }
};

// Genotype likelihoods: the expected allele counts of the slice of locus j
// and individuals [i0, i0 + ni), given their likelihoods and s.qp (already
// computed by slice_geno()). With likelihoods (1, e1, e2) of 0, 1 and 2
// copies of allele2 (from the source, see radmix_geno_access.h; allele2 is
// often the minor allele) and x = qp, the frequency of allele2,
// the probabilities of the three genotypes given the data and the
// estimates are proportional to
//   w0 = (1 - x)^2,  w1 = 2 e1 x (1 - x),  w2 = e2 x^2,
// and the expected counts of allele2 and allele1 are
//   g = (w1 + 2 w2) / f,  h = (2 w0 + w1) / f,  f = w0 + w1 + w2
// (g + h = 2; h is computed directly rather than as 2 - g, which would
// lose precision when g is close to 2). The kernels then use g and h as
// they use genotypes, which gives:
// - the EM step of NGSadmix (Skotte et al. 2013): the E step over the
//   unobserved genotypes and over the ancestry of each allele copy;
// - d1 = g / x - h / (1 - x) = f'(x) / f(x), the derivative of the
//   log-likelihood log f(x) of the genotype likelihoods, so the gradients
//   are exact;
// - d2 = g / x^2 + h / (1 - x)^2 >= 0, the negative second derivative of the
//   expected log-likelihood of the genotypes (the function that the EM step
//   maximises), rather than of log f(x), which is not concave in x (its
//   Hessian can be indefinite, and the quadratic programs of the Newton
//   updates need a positive semi-definite one). Each Newton update is then
//   a Newton step on that concave function, which has the same gradient as
//   the log-likelihood at the current estimates; the quasi-Newton step is
//   accepted on the log-likelihood itself (loglik_gl_slice()).
//   This d2 exceeds the curvature of log f(x) by the missing information
//   (the uncertainty of the genotype), so its steps are short where the
//   genotypes are uncertain. With gl_curvature = "observed" (see
//   slice_coef()), the Newton updates use the curvature of log f(x)
//   itself where it is large enough. For it, rf = obs / f is kept in the
//   slice.
// f > 0 always (the source gives positive likelihoods), and x of exactly 0
// or 1 gives finite counts (g = 0 or h = 0): with likelihoods no genotype is
// impossible. An individual without data gets g = h = 0, so that it drops
// out of every sum, as a missing genotype.
// The kind of data is a template parameter (GL) of the functions that
// process slices (slice_geno(), slice_coef(), the accumulation functions and
// sum_slices()), chosen once per block of work by the drivers (q_side(),
// p_side(), loglik()), so that the passes over genotypes are compiled as if
// genotype likelihoods did not exist: testing the kind in every slice made
// them about 3% slower (and more when the decoding of the likelihoods was
// inlined into them, which the call through the source now rules out).
//
// With idx (the free individuals of the Q side of a supervised fit, see
// slice_geno()), entry ii of the slice is individual idx[ii]: the source
// fills a run of consecutive individuals in one call, as without idx, and
// scattered ones one at a time.
//
// With a cross-validation mask (a training or held-out view), an individual
// whose likelihoods at this locus are hidden gets obs = 0, as one without
// data: the folds are those of genotypes (cv_visible(), or one word of the
// bit mask per slice), so the likelihoods of hard calls give the same fits
// and deviances as the genotypes.
static RADMIX_INLINE void slice_gl(const GenoView& geno, size_t j,
                                   size_t i0, size_t ni, Slice& s,
                                   const size_t* idx = nullptr) {
  if (idx == nullptr || idx[ni - 1] - idx[0] == ni - 1) {
    geno.gl_slice(idx == nullptr ? i0 : idx[0], ni, j, s.e1.data(),
                  s.e2.data(), s.obs.data());
  } else {
    for (size_t ii = 0; ii < ni; ii++)
      geno.gl_slice(idx[ii], 1, j, s.e1.data() + ii, s.e2.data() + ii,
                    s.obs.data() + ii);
  }
  if (geno.cv_masked()) {
    double* __restrict ob = s.obs.data();
    if (idx == nullptr && geno.cv_bits()) {
      // (i0 is a multiple of TILE_I = 64 without idx, see slice_geno())
      const uint64_t w = geno.cv_word(i0, j);
      for (size_t ii = 0; ii < ni; ii++)
        ob[ii] *= (double) ((w >> ii) & 1);
    } else {
      for (size_t ii = 0; ii < ni; ii++) {
        const size_t i = idx == nullptr ? i0 + ii : idx[ii];
        ob[ii] *= geno.cv_visible(i, j) ? 1.0 : 0.0;
      }
    }
  }
  const double* __restrict qp = s.qp.data();
  const double* __restrict e1 = s.e1.data();
  const double* __restrict e2 = s.e2.data();
  const double* __restrict obs = s.obs.data();
  double* __restrict g = s.g.data();
  double* __restrict h = s.h.data();
  double* __restrict rfs = s.rf.data();
#pragma omp simd
  for (size_t ii = 0; ii < ni; ii++) {
    const double x = qp[ii], y = 1.0 - x;
    const double w0 = y * y, w1 = 2.0 * e1[ii] * x * y, w2 = e2[ii] * x * x;
    // (obs zeroes the counts of an individual without data, without a
    // branch)
    const double rf = obs[ii] / (w0 + w1 + w2);
    g[ii] = (w1 + 2.0 * w2) * rf;
    h[ii] = (2.0 * w0 + w1) * rf;
    rfs[ii] = rf;
  }
}

// The floor of the observed curvature of genotype likelihoods, as a
// fraction of the complete-data curvature (see slice_coef()).
static const double GL_CURV_FLOOR = 0.25;

// qp and the decoded genotypes of the slice of locus j and individuals
// [i0, i1): s.qp, s.g and s.h for ii = i - i0 in [0, i1 - i0). qT is the
// transposed q (I x K, see transpose()), pj the K frequencies of locus j
// (p + K j); i1 - i0 must be at most TILE_I. A missing (or hidden) genotype
// gets g = h = 0, so that it drops out of every sum without a test. With
// genotype likelihoods (GL true, for a view with is_gl()), g and h are the
// expected allele counts given qp (see slice_gl()).
//
// With idx (the free individuals of the Q side, see q_side()), [i0, i1) are
// rows of a compact transpose qT with ldq rows (transpose_subset()), and
// entry ii of the slice is individual idx[ii] of the data (idx points to
// the entry of row i0). Without it (idx = NULL, the default), ldq is I and
// the individuals are i0 to i1 - 1.
template <bool GL>
static RADMIX_INLINE void slice_geno(const GenoView& geno,
                                     const double* qT, const double* pj,
                                     int K, size_t j, size_t i0, size_t i1,
                                     Slice& s, const size_t* idx = nullptr,
                                     size_t ldq = 0) {
  const size_t ni = i1 - i0;
  const size_t I = idx == nullptr ? geno.nrow() : ldq;
  // (__restrict tells the compiler that these buffers do not overlap, which
  // it needs to vectorise the loops below)
  double* __restrict qp = s.qp.data();
  double* __restrict g = s.g.data();
  double* __restrict h = s.h.data();
  // qp = sum_k qT[, k] p_kj, one population at a time: each pass is an axpy
  // over contiguous memory, which vectorises (a loop over k innermost would
  // be a short reduction with stride I)
  std::fill(qp, qp + ni, 0.0);
  for (int k = 0; k < K; k++) {
    const double pk = pj[k];
    const double* __restrict qk = qT + I * k + i0;
#pragma omp simd
    for (size_t ii = 0; ii < ni; ii++) qp[ii] += qk[ii] * pk;
  }
  // genotype likelihoods (a constant: see slice_gl())
  if (GL) {
    slice_gl(geno, j, i0, ni, s, idx);
    return;
  }
  // Decoding the genotypes needs table lookups, which do not vectorise, so
  // it has loops of its own. Missing and hidden genotypes are zeroed by
  // multiplying by 0 rather than by a branch, which would be mispredicted
  // when they are scattered at random (as the folds of cross-validation
  // are).
  // Fast path: when the individuals of the slice are adjacent rows of the
  // backing file and are decoded alike (no pseudohaploids), their bytes are
  // read directly and decoded with two tables giving g and h (both 0 for a
  // missing genotype). This gives the same values as the general loop below.
  // (with idx, the individuals of the slice are a run of consecutive
  // individuals of the genotypes when their indices are, and can then take
  // the paths below from the first of them, gi0)
  const size_t gi0 = idx == nullptr ? i0 : idx[0];
  const bool run = idx == nullptr || idx[ni - 1] - idx[0] == ni - 1;
  if (run) {
    const double* gtab;
    const double* htab;
    const unsigned char* __restrict bytes =
      geno.slice_bytes(gi0, ni, j, gtab, htab);
    if (bytes != nullptr) {
      if (!geno.cv_masked()) {
        for (size_t ii = 0; ii < ni; ii++) {
          const unsigned char b = bytes[ii];
          g[ii] = gtab[b];
          h[ii] = htab[b];
        }
      } else if (idx == nullptr && geno.cv_bits()) {
        // a training view with its precomputed bit mask: one word gives the
        // visibility of the whole slice (i0 is a multiple of TILE_I = 64;
        // with idx, gi0 need not be, and the bits are read one by one below)
        const uint64_t w = geno.cv_word(i0, j);
        for (size_t ii = 0; ii < ni; ii++) {
          const unsigned char b = bytes[ii];
          const double vis = (double) ((w >> ii) & 1);
          g[ii] = vis * gtab[b];
          h[ii] = vis * htab[b];
        }
      } else {
        for (size_t ii = 0; ii < ni; ii++) {
          const unsigned char b = bytes[ii];
          const double vis = geno.cv_visible(gi0 + ii, j) ? 1.0 : 0.0;
          g[ii] = vis * gtab[b];
          h[ii] = vis * htab[b];
        }
      }
      return;
    }
    // consecutive individuals whose rows are not adjacent in the backing
    // file (after filtering or reordering the individuals): when every
    // individual is decoded alike, their bytes are gathered from the column
    // of the locus through their rows and decoded with the same tables (the
    // same values as the general loop below, at a fraction of its cost)
    const unsigned char* __restrict col = geno.locus_bytes(j, gtab, htab);
    if (col != nullptr) {
      if (!geno.cv_masked()) {
        for (size_t ii = 0; ii < ni; ii++) {
          const unsigned char b = col[geno.row(gi0 + ii)];
          g[ii] = gtab[b];
          h[ii] = htab[b];
        }
      } else if (idx == nullptr && geno.cv_bits()) {
        // the bit mask, as in the fast path above
        const uint64_t w = geno.cv_word(i0, j);
        for (size_t ii = 0; ii < ni; ii++) {
          const unsigned char b = col[geno.row(gi0 + ii)];
          const double vis = (double) ((w >> ii) & 1);
          g[ii] = vis * gtab[b];
          h[ii] = vis * htab[b];
        }
      } else {
        for (size_t ii = 0; ii < ni; ii++) {
          const unsigned char b = col[geno.row(gi0 + ii)];
          const double vis = geno.cv_visible(gi0 + ii, j) ? 1.0 : 0.0;
          g[ii] = vis * gtab[b];
          h[ii] = vis * htab[b];
        }
      }
      return;
    }
    // general path: pseudohaploids (individuals decoded differently)
    double gij;
    for (size_t ii = 0; ii < ni; ii++) {
      const double obs = geno.get(gi0 + ii, j, gij) ? 1.0 : 0.0;
      g[ii] = obs * gij;
      h[ii] = obs * (geno.n_trials(gi0 + ii) - gij);
    }
    return;
  }
  // free individuals that are not consecutive: each genotype is read
  // through the index of its individual (the folds of cross-validation
  // follow these indices, so they are the same as in the other passes).
  // When every individual is decoded alike, their bytes are read from the
  // column of the locus through their rows, and decoded with the tables of
  // the fast path (the same values as the general loop below).
  const double* gtab;
  const double* htab;
  const unsigned char* __restrict col = geno.locus_bytes(j, gtab, htab);
  if (col != nullptr) {
    if (!geno.cv_masked()) {
      for (size_t ii = 0; ii < ni; ii++) {
        const unsigned char b = col[geno.row(idx[ii])];
        g[ii] = gtab[b];
        h[ii] = htab[b];
      }
    } else {
      for (size_t ii = 0; ii < ni; ii++) {
        const unsigned char b = col[geno.row(idx[ii])];
        const double vis = geno.cv_visible(idx[ii], j) ? 1.0 : 0.0;
        g[ii] = vis * gtab[b];
        h[ii] = vis * htab[b];
      }
    }
    return;
  }
  double gij;
  for (size_t ii = 0; ii < ni; ii++) {
    const size_t i = idx[ii];
    const double obs = geno.get(i, j, gij) ? 1.0 : 0.0;
    g[ii] = obs * gij;
    h[ii] = obs * (geno.n_trials(i) - gij);
  }
}

// qp, genotypes and coefficients of a slice (see slice_geno()): a and b if
// NEWTON is false (EM step), d1 and d2 if it is true (gradients and
// Hessians). A term whose genotype count is 0 contributes 0, even when qp is
// exactly 0 or 1: the count is divided by 1 instead of qp (or 1 - qp), so
// that 0 / 0 never gives NaN. A count above 0 with qp of exactly 0 (or 1)
// does give an infinite coefficient: that genotype is impossible under the
// current estimates (the Newton update then discards the step of that
// column, see newton_column()).
template <bool NEWTON, bool GL>
static RADMIX_INLINE void slice_coef(const GenoView& geno,
                                     const double* qT, const double* pj,
                                     int K, size_t j, size_t i0, size_t i1,
                                     Slice& s, const size_t* idx = nullptr,
                                     size_t ldq = 0) {
  slice_geno<GL>(geno, qT, pj, K, j, i0, i1, s, idx, ldq);
  const size_t ni = i1 - i0;
  const double* __restrict qp = s.qp.data();
  const double* __restrict g = s.g.data();
  const double* __restrict h = s.h.data();
  double* __restrict a = s.a.data();
  double* __restrict b = s.b.data();
  double* __restrict d1 = s.d1.data();
  double* __restrict d2 = s.d2.data();
#pragma omp simd
  for (size_t ii = 0; ii < ni; ii++) {
    // x = (g > 0) ? qp : 1 and y = (h > 0) ? 1 - qp : 1, written with
    // arithmetic masks so that the loop has no branches
    const double mg = (double) (g[ii] > 0), mh = (double) (h[ii] > 0);
    const double x = mg * qp[ii] + (1.0 - mg);
    const double y = mh * (1.0 - qp[ii]) + (1.0 - mh);
    if (NEWTON) {
      // two divisions instead of four: g / x = g * (1 / x), and
      // g / x^2 = (g / x) * (1 / x)
      const double rx = 1.0 / x, ry = 1.0 / y;
      const double ax = g[ii] * rx, by = h[ii] * ry;
      d1[ii] = ax - by;
      d2[ii] = ax * rx + by * ry;
    } else {
      a[ii] = g[ii] / x;
      b[ii] = h[ii] / y;
    }
  }
  // Genotype likelihoods with gl_curvature = "observed": the negative
  // second derivative of log f(x) itself, d1^2 - f''(x) / f(x) with
  // f''(x) = 2 (1 - 2 e1 + e2) (see slice_gl()), instead of the
  // complete-data d2 above. The observed curvature can be small, or
  // negative (log f(x) is not concave), where the genotypes are uncertain,
  // so it is floored at GL_CURV_FLOOR times the complete-data d2, which
  // keeps the Hessians positive semi-definite and the steps bounded. An
  // individual without data keeps d2 = 0 (d1 = rf = 0). On simulated data
  // at depths 1-2 and K = 3-10 (desktop, Nov 2026), fits took about 0.7
  // times as long (fewer iterations, at the same cost per iteration) and
  // reached maxima at least as good; without a floor, they were a little
  // faster but found the best maxima a little less often.
  if (NEWTON && GL && geno.observed_curvature()) {
    const double* __restrict e1 = s.e1.data();
    const double* __restrict e2 = s.e2.data();
    const double* __restrict rf = s.rf.data();
#pragma omp simd
    for (size_t ii = 0; ii < ni; ii++) {
      const double obs_curv =
        d1[ii] * d1[ii] - 2.0 * (1.0 - 2.0 * e1[ii] + e2[ii]) * rf[ii];
      d2[ii] = std::max(obs_curv, GL_CURV_FLOOR * d2[ii]);
    }
  }
}

// x^e for e in {0, 1, 2}, without branches, where x is already 1 when
// e = 0 (as in loglik_slice()): x when e < 2, x x when e = 2 (the factor
// in brackets is x or exactly 1, as m2 = e (e - 1) / 2 is 0, 0 or 1). For
// other values of e the result is meaningless, but loglik_slice() never
// uses it then. (Only arithmetic: a mask written as (double) (e > 1) keeps
// the loop of loglik_slice() from being vectorised in the baseline x86-64
// code, as do min() and max().)
static inline double pow012(double x, double e) {
  const double m2 = 0.5 * e * (e - 1.0);
  return x * (m2 * x + (1.0 - m2));
}

// The log-likelihood of a slice takes LOG_BATCH logarithms instead of one or
// two per genotype: it multiplies the likelihoods of LOG_BATCH genotypes
// before taking a logarithm. Each factor is at least LOG_MIN_Z (smaller ones
// are handled separately) and at most 1, so a product of 8 is at least
// 1e-240, far from underflow (the smallest normal double is about
// 2.2e-308). The slice of TILE_I individuals is split into LOG_BATCH
// products of LOG_BATCH factors each.
static const int LOG_BATCH = 8;
static const double LOG_MIN_Z = 1e-30;
static_assert(LOG_BATCH * LOG_BATCH == TILE_I,
              "loglik_slice() requires LOG_BATCH^2 == TILE_I");

// The log-likelihood of the ni genotypes of a slice (after slice_geno()),
// sum_ii g log(qp) + h log(1 - qp) (a term with a zero count is 0). For a
// hard call, g and h are 0, 1 or 2, and the contribution is log(z) with
// z = qp^g (1 - qp)^h, which needs no logarithm to compute; the z are
// multiplied in batches before taking logarithms (see LOG_BATCH). Other
// genotypes (dosages, from imputed or rescaled genotypes) and the rare z
// below LOG_MIN_Z use the direct formula. A missing genotype has
// g = h = 0, hence z = 1. Uses s.a and s.b as scratch for z and for a flag
// (all TILE_I entries of s.a, as pass 2 pads them). ni is at most TILE_I.
static RADMIX_INLINE double loglik_slice_impl(size_t ni, Slice& s) {
  const double* __restrict qp = s.qp.data();
  const double* __restrict g = s.g.data();
  const double* __restrict h = s.h.data();
  double* __restrict z = s.a.data();
  double* __restrict c = s.b.data();
  // Pass 1 (vectorised): z = qp^g (1 - qp)^h, assuming that g and h are 0,
  // 1 or 2, and a flag c that is 0 if and only if they are (c_sum = 0 if
  // all of them are). As in slice_coef(), x = qp if g > 0 and 1 otherwise,
  // and y = 1 - qp if h > 0 and 1 otherwise, so that pow012() only has to
  // choose between x and x^2 (z is exactly qp, qp^2, 1 - qp, ... as a
  // product of the same factors). For counts of 0, 1 or 2, the mask g > 0
  // is g (3 - g) / 2 (0, 1, 1; other counts are flagged by c), written with
  // arithmetic only so that the loop vectorises (see pow012()).
  double c_sum = 0.0;
#pragma omp simd reduction(+:c_sum)
  for (size_t ii = 0; ii < ni; ii++) {
    const double mg = 0.5 * g[ii] * (3.0 - g[ii]);
    const double mh = 0.5 * h[ii] * (3.0 - h[ii]);
    const double x = mg * qp[ii] + (1.0 - mg);
    const double y = mh * (1.0 - qp[ii]) + (1.0 - mh);
    z[ii] = pow012(x, g[ii]) * pow012(y, h[ii]);
    c[ii] = std::fabs(g[ii] * (g[ii] - 1.0) * (g[ii] - 2.0)) +
      std::fabs(h[ii] * (h[ii] - 1.0) * (h[ii] - 2.0));
    c_sum += c[ii];
  }
  // Pass 2: LOG_BATCH interleaved products (product b collects z[b],
  // z[b + 8], z[b + 16], ...), so that the inner loop multiplies 8
  // independent accumulators by 8 contiguous values and vectorises. A short
  // slice is padded with ones (log(1) = 0) to keep the loop bounds fixed.
  std::fill(z + ni, z + TILE_I, 1.0);
  double prod[LOG_BATCH];
  for (int b = 0; b < LOG_BATCH; b++) prod[b] = 1.0;
  for (int m = 0; m < LOG_BATCH; m++)
#pragma omp simd
    for (int b = 0; b < LOG_BATCH; b++) prod[b] *= z[LOG_BATCH * m + b];
  // Pass 3 (scalar, rarely taken): the genotypes that cannot go into the
  // products, dosages (c > 0) and z below LOG_MIN_Z, get explicit logarithms
  // and a neutral z = 1, and the products are formed again. This includes
  // z = 0 (an impossible genotype, qp of exactly 0 with g > 0, or 1 with
  // h > 0), which correctly gives -Inf. The factors of hard calls are at
  // most 1, so a z below LOG_MIN_Z makes its product smaller than LOG_MIN_Z:
  // the pass is needed only if c_sum > 0 or a product is below LOG_MIN_Z,
  // which never happens for hard calls with estimates inside
  // [1e-5, 1 - 1e-5] (z >= 1e-10). Otherwise the result is the same as if
  // the pass had run (it would not have changed anything). (| rather than ||
  // avoids a second branch.)
  double ll = 0.0;
  double p_min = prod[0];
  for (int b = 1; b < LOG_BATCH; b++) p_min = std::min(p_min, prod[b]);
  if ((c_sum != 0) | !(p_min >= LOG_MIN_Z)) {
    for (size_t ii = 0; ii < ni; ii++) {
      if ((z[ii] < LOG_MIN_Z) | (c[ii] != 0)) {
        z[ii] = 1.0;
        if (g[ii] > 0) ll += g[ii] * std::log(qp[ii]);
        if (h[ii] > 0) ll += h[ii] * std::log(1.0 - qp[ii]);
      }
    }
    for (int b = 0; b < LOG_BATCH; b++) prod[b] = 1.0;
    for (int m = 0; m < LOG_BATCH; m++)
#pragma omp simd
      for (int b = 0; b < LOG_BATCH; b++) prod[b] *= z[LOG_BATCH * m + b];
  }
  for (int b = 0; b < LOG_BATCH; b++) ll += std::log(prod[b]);
  return ll;
}

// The log-likelihood of the ni individuals of a slice of genotype
// likelihoods (after slice_geno()): the sum of log(f / (1 + e1 + e2)) over
// the individuals with data, where f = (1 - qp)^2 + 2 e1 qp (1 - qp) +
// e2 qp^2 (see slice_gl()). Dividing by 1 + e1 + e2 normalises the
// likelihoods of each individual to sum to 1, as in the beagle files read
// by NGSadmix, so that the log-likelihood does not depend on the constant
// factor of the likelihoods of the source; it is the log of the
// probability of the data under the estimates, up to the probability of the
// reads given each genotype. Individuals without data contribute 0. Each
// factor z = f / (1 + e1 + e2) is a weighted mean of the three genotype
// probabilities, at most 1, and is far from LOG_MIN_Z in practice (as long
// as the likelihoods are not many orders of magnitude apart, e.g. within
// about 32 nats of the best, as the tests clamp them), so the factors are
// multiplied in batches before taking logarithms, as in
// loglik_slice(); smaller ones take the direct formula. Uses s.a as scratch
// (all TILE_I entries). ni is at most TILE_I.
static RADMIX_INLINE double loglik_gl_slice(size_t ni, Slice& s) {
  const double* __restrict qp = s.qp.data();
  const double* __restrict e1 = s.e1.data();
  const double* __restrict e2 = s.e2.data();
  const double* __restrict obs = s.obs.data();
  double* __restrict z = s.a.data();
  // z = f / (1 + e1 + e2), and 1 for an individual without data
#pragma omp simd
  for (size_t ii = 0; ii < ni; ii++) {
    const double x = qp[ii], y = 1.0 - x;
    const double f = y * y + 2.0 * e1[ii] * x * y + e2[ii] * x * x;
    z[ii] = obs[ii] * (f / (1.0 + e1[ii] + e2[ii])) + (1.0 - obs[ii]);
  }
  // the rare factors that cannot go into the products (see loglik_slice())
  double ll = 0.0;
  for (size_t ii = 0; ii < ni; ii++) {
    if (z[ii] < LOG_MIN_Z) {
      ll += std::log(z[ii]);
      z[ii] = 1.0;
    }
  }
  // LOG_BATCH interleaved products, padded with ones (see loglik_slice())
  std::fill(z + ni, z + TILE_I, 1.0);
  double prod[LOG_BATCH];
  for (int b = 0; b < LOG_BATCH; b++) prod[b] = 1.0;
  for (int m = 0; m < LOG_BATCH; m++)
#pragma omp simd
    for (int b = 0; b < LOG_BATCH; b++) prod[b] *= z[LOG_BATCH * m + b];
  for (int b = 0; b < LOG_BATCH; b++) ll += std::log(prod[b]);
  return ll;
}

// loglik_slice() for loglik() and deviance() (sum_slices(), baseline code)
static double loglik_slice(size_t ni, Slice& s) {
  return loglik_slice_impl(ni, s);
}

// loglik_slice() for the Q-side passes that also return the log-likelihood
// (accum_q_block() with LOGLIK), inlined into accum_q_block(), so that it is
// also compiled for AVX2 in its clone. Its logarithms are calls to log() of
// the C library either way (8 per slice of 64 genotypes, or more for
// dosages), and its products and masks are exact, so each slice gets the
// same value as in loglik() (up to the rounding of dosages, whose terms can
// be contracted into FMAs). This made the fused pass about 20% faster than
// with a call to the baseline version (laptop, HGDP, K = 5); the AVX2
// clones of the whole log-likelihood pass that were measured slower before
// (see RADMIX_MULTIVERSION) also called log() for every genotype of a slice
// with a dosage or a tiny likelihood. Defining RADMIX_FUSED_LL_NOINLINE when
// compiling calls the baseline version instead, for benchmarks. With
// genotype likelihoods (GL), it calls loglik_gl_slice(), as loglik() does.
#if defined(__GNUC__) && defined(RADMIX_FUSED_LL_NOINLINE)
template <bool GL>
__attribute__((noinline)) static double loglik_slice_q(size_t ni, Slice& s) {
  return GL ? loglik_gl_slice(ni, s) : loglik_slice_impl(ni, s);
}
#else
template <bool GL>
static RADMIX_INLINE double loglik_slice_q(size_t ni, Slice& s) {
  return GL ? loglik_gl_slice(ni, s) : loglik_slice_impl(ni, s);
}
#endif

// Q side: the sums over loci [j0, j1) of each individual of the block
// [i0, i1) (all the loci, unless the work is split into segments, see
// q_side()), and, if LOGLIK, the log-likelihood of their genotypes, which
// is returned (0 without LOGLIK):
// - NEWTON: grad_k = sum_j d1 p_kj (grad, K x ni) and
//   hess_kk2 = sum_j d2 p_kj p_k2j (hess, K x K x ni), the gradient and the
//   negative Hessian of the log-likelihood with respect to q_.i;
// - EM: grad_k = sum_j a p_kj + b (1 - p_kj) (grad, K x ni; hess unused),
//   the expected number of allele copies of individual i from population
//   k, divided by q_ki (see em()).
// - LOGLIK without NEWTON: the log-likelihood alone, with no sums (grad,
//   hess, gT and hT are not used): the slices only need qp and the
//   genotypes (slice_geno()). Used for the reference individuals of a
//   supervised fit, whose columns of q are not updated (update_q_impl()).
// Outputs are overwritten. qT is the transposed q (I x K), p is K x J, and
// i1 - i0 must be at most TILE_I. s holds two slices (one per locus of a
// step); gT (K vectors of TILE_I values) and hT (K (K + 1) / 2 vectors of
// TILE_I values, only used if NEWTON) are scratch for the accumulators.
// With a list of free individuals (free_idx, see q_side()), [i0, i1) are
// positions in that list, and qT is their compact transpose, with ldq rows.
// init, out and ll0 are true, true and 0 except in a streaming pass, which
// accumulates a block over several chunks of loci: init = false resumes the
// accumulators gT and hT as they are, out = false leaves the sums in them
// (grad and hess are then not written), and the log-likelihood returned
// adds that of the block to ll0.
//
// The accumulators are kept "transposed": one vector of length TILE_I per
// gradient entry k, and one per Hessian entry (k, k2) with k2 >= k, so that
// each update is a unit-stride loop over the individuals of the slice,
// vectorised across individuals. Vector k (or pair idx, in the order (0,0),
// (0,1), ..., (0,K-1), (1,1), ..., (K-1,K-1)) starts at gT + TILE_I k (or
// hT + TILE_I idx), whatever ni is; a short last block simply leaves the
// tail of each vector unused. They are copied into the output layout once,
// at the end. Each step of the loop over the loci reads and writes every
// accumulator once for two loci (see "Layout of the work" at the top).
// BLOCK (with NEWTON): the Hessian sums are register-blocked (see
// HESS_BLOCK_MIN_K), with the buffers in hB (HB_L x (TILE_I + npair_pad)
// doubles).
template <bool NEWTON, bool GL, bool LOGLIK, bool BLOCK>
RADMIX_MULTIVERSION
static double accum_q_block(const GenoView& geno, const double* qT,
                            const double* p, int K, size_t i0, size_t i1,
                            size_t j0, size_t j1, double* grad,
                            double* hess, Slice* s, double* gT, double* hT,
                            double* hB, const size_t* free_idx, size_t ldq,
                            bool init, bool out, double ll0) {
  const size_t ni = i1 - i0, KK = (size_t) K * K;
  // the log-likelihood alone (see above): qp and the genotypes of each
  // slice, and its log-likelihood, in the order of the loci as below (EM
  // steps, which use s.a and s.b, cannot also compute it, as loglik_slice()
  // overwrites them)
  if (!NEWTON && LOGLIK) {
    const size_t* fi = free_idx == nullptr ? nullptr : free_idx + i0;
    double ll = ll0;
    for (size_t j = j0; j < j1; j++) {
      slice_geno<GL>(geno, qT, p + K * j, K, j, i0, i1, s[0], fi, ldq);
      ll += loglik_slice_q<GL>(ni, s[0]);
    }
    return ll;
  }
  // the number of entries (k, k2) with k2 >= k of a symmetric K x K matrix
  const size_t n_pair = (size_t) K * (K + 1) / 2;
  // (the whole vectors are zeroed, including the tail beyond ni, which is
  // never read; a streaming pass resumes the accumulators of the block from
  // the previous chunk instead, see "Streaming" at the top)
  if (init) {
    std::fill(gT, gT + K * TILE_I, 0.0);
    if (NEWTON) std::fill(hT, hT + n_pair * TILE_I, 0.0);
  }
  // the buffers of the blocked Hessians, and the number of loci in them
  // (added into hT at the end of every call, so nothing is carried over
  // between the chunks of a streaming pass; its chunks are whole blocks of
  // TILE_J loci, a multiple of HB_L, so the loci are buffered in the same
  // groups as in one run over all of them, and the sums are the same)
  const size_t npair_pad = hess_pairs_pad(n_pair);
  double* d2b = hB;
  double* ppb = BLOCK ? hB + HB_L * TILE_I : nullptr;
  int nl = 0;
  // the coefficients of the two loci of each step (s[0] for locus j, s[1]
  // for locus j + 1); a and b are used by EM, d1 and d2 by Newton
  const double* __restrict a1 = s[0].a.data();
  const double* __restrict b1 = s[0].b.data();
  const double* __restrict a2 = s[1].a.data();
  const double* __restrict b2 = s[1].b.data();
  const double* __restrict d11 = s[0].d1.data();
  const double* __restrict d21 = s[0].d2.data();
  const double* __restrict d12 = s[1].d1.data();
  const double* __restrict d22 = s[1].d2.data();
  // the log-likelihood of the slices of the block, in the order of the loci
  // (continuing from ll0, that of the previous chunks of a streaming pass)
  double ll = ll0;
  // with a list of free individuals, those of the block (see q_side())
  const size_t* fi = free_idx == nullptr ? nullptr : free_idx + i0;
  for (size_t j = j0; j < j1; j += 2) {
    // the frequencies of the two loci; the padding locus of an odd range
    // (j + 1 == j1) has zero coefficients, so any frequencies will do
    const bool two = j + 1 < j1;
    const double* pj = p + K * j;
    const double* pj2 = two ? pj + K : pj;
    slice_coef<NEWTON, GL>(geno, qT, pj, K, j, i0, i1, s[0], fi, ldq);
    if (two)
      slice_coef<NEWTON, GL>(geno, qT, pj2, K, j + 1, i0, i1, s[1], fi, ldq);
    else s[1].zero_coef(ni);
    // The log-likelihood comes from the same qp and genotypes as the
    // coefficients (slice_coef() leaves s.qp, s.g and s.h in place, and,
    // for genotype likelihoods, s.e1, s.e2 and s.obs), so it costs no extra
    // reading of the data. loglik_slice() overwrites s.a and s.b, which
    // Newton steps do not use (the padding locus has no genotypes and is
    // skipped).
    if (LOGLIK) {
      ll += loglik_slice_q<GL>(ni, s[0]);
      if (two) ll += loglik_slice_q<GL>(ni, s[1]);
    }
    if (NEWTON) {
      // gradient: gT[, k] += d1 p_kj, for both loci
      for (int k = 0; k < K; k++) {
        const double pk = pj[k], pk2 = pj2[k];
        double* __restrict gk = gT + TILE_I * k;
#pragma omp simd
        for (size_t ii = 0; ii < ni; ii++)
          gk[ii] += d11[ii] * pk + d12[ii] * pk2;
      }
      // Hessian: hT[, (k, k2)] += d2 p_kj p_k2j for k2 >= k
      if (BLOCK) {
        // buffer the two loci (the padding locus has d2 = 0), and add up
        // the buffered loci every HB_L loci and at the end
        hess_buffer(d21, pj, K, ni, nl, d2b, ppb, npair_pad);
        hess_buffer(d22, pj2, K, ni, nl + 1, d2b, ppb, npair_pad);
        nl += 2;
        if (nl == HB_L || j + 2 >= j1) {
          hess_block(d2b, ppb, nl, n_pair, npair_pad, ni, hT);
          nl = 0;
        }
      } else {
        // the products of the frequencies do not depend on the
        // individual, so they are formed outside the loop over individuals
        double* __restrict hk = hT;
        for (int k = 0; k < K; k++) {
          for (int k2 = k; k2 < K; k2++, hk += TILE_I) {
            const double c1 = pj[k] * pj[k2], c2 = pj2[k] * pj2[k2];
#pragma omp simd
            for (size_t ii = 0; ii < ni; ii++)
              hk[ii] += d21[ii] * c1 + d22[ii] * c2;
          }
        }
      }
    } else {
      // EM: gT[, k] += a p_kj + b (1 - p_kj), for both loci (pk, ok: p_kj
      // and 1 - p_kj of locus j; pk2, ok2: those of locus j + 1)
      for (int k = 0; k < K; k++) {
        const double pk = pj[k], ok = 1.0 - pj[k];
        const double pk2 = pj2[k], ok2 = 1.0 - pj2[k];
        double* __restrict gk = gT + TILE_I * k;
#pragma omp simd
        for (size_t ii = 0; ii < ni; ii++)
          gk[ii] += (a1[ii] * pk + b1[ii] * ok) +
            (a2[ii] * pk2 + b2[ii] * ok2);
      }
    }
  }
  // copy the accumulators into the output layout: grad K x ni and hess
  // K x K x ni, one column (or K x K matrix) per individual; the Hessian
  // entries are read in the order in which they were accumulated (idx), into
  // the lower triangle, and then mirrored (not before the last chunk of a
  // streaming pass)
  if (!out) return ll;
  for (size_t ii = 0; ii < ni; ii++) {
    for (int k = 0; k < K; k++) grad[k + K * ii] = gT[ii + TILE_I * k];
    if (NEWTON) {
      double* H = hess + KK * ii;
      size_t idx = 0;
      for (int k = 0; k < K; k++)
        for (int k2 = k; k2 < K; k2++, idx++)
          H[k2 + K * k] = hT[ii + TILE_I * idx];
      mirror_hess(H, K);
    }
  }
  return ll;
}

// Two-dimensional blocking of the Q side. Q-side sums are split into blocks
// of TILE_I individuals, each owned by one thread; but with few individuals
// there are few blocks (379 individuals give 6, 50 individuals a single
// one), which cannot keep many threads busy. Population genetic data often
// have tens to hundreds of individuals and many loci, so when there are
// fewer than MIN_Q_TASKS blocks the loci are also split into segments, and
// each (block, segment) pair is a separate task (see q_side()).
//
// Memory: the partial sums of the segments take n_seg I (K + K^2) doubles.
// Segments are only used with fewer than MIN_Q_TASKS blocks, i.e.
// I < MIN_Q_TASKS TILE_I = 8192, and n_seg = ceil(MIN_Q_TASKS / n_blocks), so
// n_seg I stays below 2 MIN_Q_TASKS TILE_I = 16384: under 15 MB for K = 10.

// the target number of Q-side tasks, enough to balance many threads: with
// 32, the HGDP (1043 individuals, 17 blocks) had 34 tasks, and fits on 8
// threads were 3-6% slower than with 128 (136 tasks; see above)
static const size_t MIN_Q_TASKS = RADMIX_MIN_Q_TASKS;
// the minimum number of loci of a segment, so that each task has enough
// work to amortise its overheads
static const size_t MIN_SEG_LOCI = 256;

// The number of segments of loci of the Q side (1 = no segments). It depends
// only on the dimensions of the data, so that results do not depend on the
// number of threads.
size_t q_segments(size_t I, size_t J) {
  const size_t nbi = n_blocks(I, TILE_I);
  // enough blocks already (or no individuals: nothing to split)
  if (nbi == 0 || nbi >= MIN_Q_TASKS) return 1;
  // enough segments for about MIN_Q_TASKS tasks (e.g. 6 blocks need 6
  // segments, giving 36 tasks), but none shorter than MIN_SEG_LOCI loci
  size_t n_seg = (MIN_Q_TASKS + nbi - 1) / nbi;
  n_seg = std::min(n_seg, std::max((size_t) 1, J / MIN_SEG_LOCI));
  return n_seg;
}

// The number of blocks of TILE_J loci per chunk of a streaming pass: about
// chunk_bytes bytes of the backing file (see GenoView::chunk_loci()), at
// least one block.
static size_t chunk_blocks(const GenoView& geno) {
  const size_t b = geno.chunk_loci() / TILE_J;
  return b < 1 ? 1 : b;
}

// The Q-side driver: the sums of accum_q_block() over all the loci for
// every individual, then finish(i0, i1, grad, hess) for each block of
// TILE_I individuals [i0, i1), with the complete sums of the block (grad
// K x ni and, if NEWTON, hess K x K x ni, ni = i1 - i0; without NEWTON,
// hess is not meaningful). finish() runs inside the parallel region, on the
// thread that owns the block, so it can write the output columns of its
// individuals without synchronisation; it must not call R. Its arguments
// point to thread-local buffers, valid only during the call. qT is the
// transposed q (I x K), p is K x J.
//
// With a list of free individuals (free_idx, n_free values, in increasing
// order; see "Fixed individuals" in radmix_kernels.h), only those are
// processed: the blocks are blocks of TILE_I consecutive entries of the
// list, qT is the compact transpose of their columns (n_free x K, see
// transpose_subset()), and [i0, i1) passed to finish() are positions in
// the list (entry f is individual free_idx[f]). The sums of the fixed
// individuals are never computed. The free individuals of a slice are read
// through their indices, unless they are consecutive (slice_geno()); the
// arithmetic of each sum is the same as without the list.
//
// With one segment, each block is accumulated over all the loci by one
// thread, into thread-local buffers. With several, each (segment, block)
// task writes its partial sums into its own slot of grad_part and
// hess_part, and the partial sums of each individual are then added up in
// segment order before finish(). Both stages run in a single parallel
// region, separated by the barrier at the end of the first loop.
//
// With LOGLIK (and NEWTON), *ll receives the log-likelihood L(q, p), from
// the same slices: each task keeps the log-likelihood of its own genotypes
// (a block of individuals over a segment of loci), and these partial sums
// are added up at the end in task order, which depends only on the
// dimensions of the data, so *ll does not depend on the number of threads.
// It equals loglik() up to rounding (the terms are the same, added up in
// another order). Without LOGLIK, ll is not used. With a list of
// individuals, LOGLIK gives the log-likelihood of their genotypes alone,
// which is only part of L(q, p) (update_q_impl() adds up those of the free
// and of the fixed individuals).
//
// LOGLIK without NEWTON computes the log-likelihood alone (see
// accum_q_block()), with the same tasks: there are no sums, and finish() is
// not called.
template <bool NEWTON, bool LOGLIK, class Finish>
static void q_side_stream(const GenoView& geno, const double* qT,
                          const double* p, int K, int n_cores, double* ll,
                          Finish finish, const size_t* free_idx,
                          size_t n_free);

template <bool NEWTON, bool LOGLIK, class Finish>
static void q_side(const GenoView& geno, const double* qT, const double* p,
                   int K, int n_cores, double* ll, Finish finish,
                   const size_t* free_idx = nullptr, size_t n_free = 0) {
  if (geno.stream()) {
    q_side_stream<NEWTON, LOGLIK>(geno, qT, p, K, n_cores, ll, finish,
                                  free_idx, n_free);
    return;
  }
  // whether the pass computes sums (anything but the log-likelihood alone)
  const bool SUMS = NEWTON || !LOGLIK;
  // I: the number of individuals processed (all of them, or those of the
  // list), which is also the number of rows of qT
  const size_t I = free_idx == nullptr ? geno.nrow() : n_free;
  const size_t J = geno.ncol(), KK = (size_t) K * K;
  const size_t n_pair = (size_t) K * (K + 1) / 2;
  const size_t nbi = n_blocks(I, TILE_I), n_seg = q_segments(I, J);
  // the partial sums of individual i in segment seg start at
  // K (seg I + i) (and KK (seg I + i) for the Hessians)
  std::vector<double> grad_part(n_seg > 1 && SUMS ? n_seg * K * I : 0);
  std::vector<double> hess_part(n_seg > 1 && NEWTON ? n_seg * KK * I : 0);
  // with LOGLIK, the log-likelihood of each task, at index seg * nbi + bl
  // (one segment: bl)
  std::vector<double> ll_part(LOGLIK ? n_seg * nbi : 0);
  // the accumulation function for the kind of data (see slice_gl())
  // (with the blocked Hessians from HESS_BLOCK_MIN_K)
  const bool block = NEWTON && K >= HESS_BLOCK_MIN_K;
  const auto accum = geno.is_gl() ?
    (block ? accum_q_block<NEWTON, true, LOGLIK, true> :
     accum_q_block<NEWTON, true, LOGLIK, false>) :
    (block ? accum_q_block<NEWTON, false, LOGLIK, true> :
     accum_q_block<NEWTON, false, LOGLIK, false>);
#pragma omp parallel num_threads(n_cores)
{
  // thread-local buffers: two slices, the accumulators, and the sums of one
  // block of individuals
  Slice s[2];
  std::vector<double> gT(SUMS ? K * TILE_I : 0);
  std::vector<double> hT(NEWTON ? n_pair * TILE_I : 0);
  std::vector<double> hB(
    block ? HB_L * (TILE_I + hess_pairs_pad(n_pair)) : 0);
  std::vector<double> grad(SUMS ? K * TILE_I : 0);
  std::vector<double> hess(NEWTON ? KK * TILE_I : 0);
  if (n_seg > 1) {
    // tasks t = seg * nbi + bl: blocks vary fastest, so threads that start
    // together read the same columns of the backing file
    const ptrdiff_t n_tasks = (ptrdiff_t) (nbi * n_seg);
#pragma omp for schedule(RADMIX_OMP_SCHEDULE)
    for (ptrdiff_t t = 0; t < n_tasks; t++) {
      const size_t seg = t / nbi, bl = t % nbi;
      const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
      // n_seg contiguous, nearly equal ranges of loci
      const size_t j0 = J * seg / n_seg, j1 = J * (seg + 1) / n_seg;
      const double ll_task = accum(
        geno, qT, p, K, i0, i1, j0, j1,
        SUMS ? &grad_part[K * (seg * I + i0)] : NULL,
        NEWTON ? &hess_part[KK * (seg * I + i0)] : NULL, s, gT.data(),
        NEWTON ? hT.data() : NULL, block ? hB.data() : NULL, free_idx, I,
        true, true, 0.0);
      if (LOGLIK) ll_part[t] = ll_task;
    }
    // (the implicit barrier at the end of the loop above guarantees that
    // all the partial sums are written before they are added up below)
  }
  // one task per block of individuals: accumulate it over all the loci, or
  // add up its partial sums, then hand the complete sums to finish()
#pragma omp for schedule(RADMIX_OMP_SCHEDULE)
  for (ptrdiff_t bl = 0; bl < (ptrdiff_t) nbi; bl++) {
    const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
    if (n_seg == 1) {
      const double ll_task = accum(
        geno, qT, p, K, i0, i1, 0, J, grad.data(),
        NEWTON ? hess.data() : NULL, s, gT.data(),
        NEWTON ? hT.data() : NULL, block ? hB.data() : NULL, free_idx, I,
        true, true, 0.0);
      if (LOGLIK) ll_part[bl] = ll_task;
    } else if (SUMS) {
      // add up the partial sums of the segments, always in segment order
      for (size_t i = i0; i < i1; i++) {
        double* g = &grad[K * (i - i0)];
        for (int k = 0; k < K; k++) {
          double sum = 0.0;
          for (size_t seg = 0; seg < n_seg; seg++)
            sum += grad_part[K * (seg * I + i) + k];
          g[k] = sum;
        }
        if (!NEWTON) continue;
        double* H = &hess[KK * (i - i0)];
        for (size_t kk = 0; kk < KK; kk++) {
          double sum = 0.0;
          for (size_t seg = 0; seg < n_seg; seg++)
            sum += hess_part[KK * (seg * I + i) + kk];
          H[kk] = sum;
        }
      }
    }
    if (SUMS) finish(i0, i1, grad.data(), hess.data());
  }
}
  // the log-likelihood, adding up the tasks in order (outside the parallel
  // region)
  if (LOGLIK) {
    double total = 0.0;
    for (size_t t = 0; t < ll_part.size(); t++) total += ll_part[t];
    *ll = total;
  }
}

// q_side() for a view in streaming mode (see "Streaming" at the top): the
// same tasks, sums and order of accumulation, walked in chunks of loci,
// with the next chunk prefetched. A Q-side pass always runs forwards.
template <bool NEWTON, bool LOGLIK, class Finish>
static void q_side_stream(const GenoView& geno, const double* qT,
                          const double* p, int K, int n_cores, double* ll,
                          Finish finish, const size_t* free_idx,
                          size_t n_free) {
  const bool SUMS = NEWTON || !LOGLIK;
  const size_t I = free_idx == nullptr ? geno.nrow() : n_free;
  const size_t J = geno.ncol(), KK = (size_t) K * K;
  const size_t n_pair = (size_t) K * (K + 1) / 2;
  const size_t nbi = n_blocks(I, TILE_I), n_seg = q_segments(I, J);
  geno.direction(true);
  // (with the blocked Hessians from HESS_BLOCK_MIN_K, as in q_side())
  const bool block = NEWTON && K >= HESS_BLOCK_MIN_K;
  const auto accum = geno.is_gl() ?
    (block ? accum_q_block<NEWTON, true, LOGLIK, true> :
     accum_q_block<NEWTON, true, LOGLIK, false>) :
    (block ? accum_q_block<NEWTON, false, LOGLIK, true> :
     accum_q_block<NEWTON, false, LOGLIK, false>);
  const size_t hb_size = block ? HB_L * (TILE_I + hess_pairs_pad(n_pair)) : 0;
  // the chunks, as ranges of loci, and (several segments) of segments
  std::vector<size_t> cut_j(1, 0), cut_seg(1, 0);
  if (n_seg == 1) {
    // an even number of loci per chunk (TILE_J is even), so that the loci
    // are paired as in one run over all of them
    const size_t C = chunk_blocks(geno) * TILE_J;
    for (size_t j = C; j < J; j += C) cut_j.push_back(j);
    cut_j.push_back(J);
  } else {
    // groups of whole segments, of about the size of a chunk
    const double seg_loci = (double) J / n_seg;
    size_t per = (size_t) std::floor(chunk_blocks(geno) * TILE_J / seg_loci);
    if (per < 1) per = 1;
    for (size_t sg = per; sg < n_seg; sg += per) {
      cut_seg.push_back(sg);
      cut_j.push_back(J * sg / n_seg);
    }
    cut_seg.push_back(n_seg);
    cut_j.push_back(J);
  }
  const size_t n_chunks = cut_j.size() - 1;
  // one segment: the accumulators and the log-likelihood of each block,
  // kept from one chunk to the next
  const bool keep = n_seg == 1 && n_chunks > 1;
  std::vector<double> acc_g(keep && SUMS ? nbi * K * TILE_I : 0);
  std::vector<double> acc_h(keep && NEWTON ? nbi * n_pair * TILE_I : 0);
  // several segments: the partial sums of each (segment, block), as in
  // q_side()
  std::vector<double> grad_part(n_seg > 1 && SUMS ? n_seg * K * I : 0);
  std::vector<double> hess_part(n_seg > 1 && NEWTON ? n_seg * KK * I : 0);
  std::vector<double> ll_part(LOGLIK ? n_seg * nbi : 0, 0.0);
  // The prefetch of the next chunk reads at most chunk_blocks() blocks of
  // loci, the size of a chunk elsewhere: a chunk of whole segments can be
  // much larger (a segment of 187 MB against 32 MB chunks), and the chunk
  // computed must stay in the cache until all its blocks of individuals are
  // done, so prefetching all of the next one could evict it (with a third
  // of the file cacheable, that read 4 GB more in a fit than no prefetch).
  // The rest of a large chunk is read by the kernel's read-ahead of the
  // threads' forward page faults.
  const size_t max_prefetch = chunk_blocks(geno) * TILE_J;
  Prefetcher prefetcher(geno);
  for (size_t c = 0; c < n_chunks; c++) {
    if (c + 1 < n_chunks)
      prefetcher.start(cut_j[c + 1],
                       std::min(cut_j[c + 2], cut_j[c + 1] + max_prefetch));
    const bool first = c == 0, last = c + 1 == n_chunks;
    if (n_seg == 1) {
      const size_t j0 = cut_j[c], j1 = cut_j[c + 1];
#pragma omp parallel num_threads(n_cores)
{
      // the accumulators are thread-local, as in q_side(), and copied from
      // and to the block's slots between chunks (accumulating straight into
      // the slots, in one large array, made the pass 15-40% slower)
      Slice s[2];
      std::vector<double> gT(SUMS ? K * TILE_I : 0);
      std::vector<double> hT(NEWTON ? n_pair * TILE_I : 0);
      std::vector<double> hB(hb_size);
      std::vector<double> grad(last && SUMS ? K * TILE_I : 0);
      std::vector<double> hess(last && NEWTON ? KK * TILE_I : 0);
#pragma omp for schedule(RADMIX_OMP_SCHEDULE)
      for (ptrdiff_t bl = 0; bl < (ptrdiff_t) nbi; bl++) {
        const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
        double* slot_g = SUMS ? &acc_g[bl * K * TILE_I] : NULL;
        double* slot_h = NEWTON ? &acc_h[bl * n_pair * TILE_I] : NULL;
        if (!first) {
          if (SUMS) std::copy(slot_g, slot_g + K * TILE_I, gT.begin());
          if (NEWTON)
            std::copy(slot_h, slot_h + n_pair * TILE_I, hT.begin());
        }
        const double l = accum(
          geno, qT, p, K, i0, i1, j0, j1,
          last && SUMS ? grad.data() : NULL,
          last && NEWTON ? hess.data() : NULL, s,
          SUMS ? gT.data() : NULL, NEWTON ? hT.data() : NULL,
          block ? hB.data() : NULL, free_idx, I, first, last,
          LOGLIK ? ll_part[bl] : 0.0);
        if (LOGLIK) ll_part[bl] = l;
        if (last && SUMS) finish(i0, i1, grad.data(), hess.data());
        if (!last) {
          if (SUMS) std::copy(gT.begin(), gT.end(), slot_g);
          if (NEWTON) std::copy(hT.begin(), hT.end(), slot_h);
        }
      }
}
    } else {
      const size_t s0 = cut_seg[c], s1 = cut_seg[c + 1];
      const ptrdiff_t t0 = (ptrdiff_t) (s0 * nbi),
        t1 = (ptrdiff_t) (s1 * nbi);
#pragma omp parallel num_threads(n_cores)
{
      Slice s[2];
      std::vector<double> gT(SUMS ? K * TILE_I : 0);
      std::vector<double> hT(NEWTON ? n_pair * TILE_I : 0);
      std::vector<double> hB(hb_size);
#pragma omp for schedule(RADMIX_OMP_SCHEDULE)
      for (ptrdiff_t t = t0; t < t1; t++) {
        const size_t seg = t / nbi, bl = t % nbi;
        const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
        const size_t j0 = J * seg / n_seg, j1 = J * (seg + 1) / n_seg;
        const double l = accum(
          geno, qT, p, K, i0, i1, j0, j1,
          SUMS ? &grad_part[K * (seg * I + i0)] : NULL,
          NEWTON ? &hess_part[KK * (seg * I + i0)] : NULL, s, gT.data(),
          NEWTON ? hT.data() : NULL, block ? hB.data() : NULL, free_idx, I,
          true, true, 0.0);
        if (LOGLIK) ll_part[t] = l;
      }
}
    }
  }
  prefetcher.wait();
  // several segments: add up the partial sums of each block in segment
  // order, as q_side() does
  if (n_seg > 1 && SUMS) {
#pragma omp parallel num_threads(n_cores)
{
    std::vector<double> grad(K * TILE_I);
    std::vector<double> hess(NEWTON ? KK * TILE_I : 0);
#pragma omp for schedule(RADMIX_OMP_SCHEDULE)
    for (ptrdiff_t bl = 0; bl < (ptrdiff_t) nbi; bl++) {
      const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
      for (size_t i = i0; i < i1; i++) {
        double* g = &grad[K * (i - i0)];
        for (int k = 0; k < K; k++) {
          double sum = 0.0;
          for (size_t seg = 0; seg < n_seg; seg++)
            sum += grad_part[K * (seg * I + i) + k];
          g[k] = sum;
        }
        if (!NEWTON) continue;
        double* H = &hess[KK * (i - i0)];
        for (size_t kk = 0; kk < KK; kk++) {
          double sum = 0.0;
          for (size_t seg = 0; seg < n_seg; seg++)
            sum += hess_part[KK * (seg * I + i) + kk];
          H[kk] = sum;
        }
      }
      finish(i0, i1, grad.data(), hess.data());
    }
}
  }
  if (LOGLIK) {
    double total = 0.0;
    for (size_t t = 0; t < ll_part.size(); t++) total += ll_part[t];
    *ll = total;
  }
}

// The free individuals of a mask of fixed individuals (see "Fixed
// individuals" in radmix_kernels.h), in increasing order, for q_side(); the
// columns of the fixed individuals of q are copied into q_new.
static std::vector<size_t> free_individuals(const unsigned char* fixed,
                                            const double* q, int K,
                                            size_t I, double* q_new) {
  std::vector<size_t> free_idx;
  for (size_t i = 0; i < I; i++) {
    if (fixed[i]) std::copy(q + K * i, q + K * (i + 1), q_new + K * i);
    else free_idx.push_back(i);
  }
  return free_idx;
}

// P side: the sums over all the individuals of each locus of the block
// [j0, j1):
// - NEWTON: out1 = grad (K x nj), grad_k = sum_i d1 q_ki, and out2 = hess
//   (K x K x nj), hess_kk2 = sum_i d2 q_ki q_k2i, the gradient and the
//   negative Hessian of the log-likelihood with respect to p_.j;
// - EM: out1 = sum_i a q_ki and out2 = sum_i b q_ki (both K x nj), which
//   em() turns into expected allele counts.
// Outputs are overwritten. They must have room for TILE_J loci: the padding
// locus of a block with an odd number of loci adds zeros to the slot after
// its last locus (such a block is the last one, as TILE_J is even, so it has
// fewer than TILE_J loci and the spare slot exists). qT is the transposed q
// (I x K), p is K x J, j1 - j0 is at most TILE_J, and s holds two slices.
//
// Loop order: tiles of TILE_I individuals (outer), then the loci of the
// block, two at a time (inner). The TILE_I x K corner of qT used by a tile is
// reused for all the loci of the block, so it stays in cache, and each sum
// is a dot product over the individuals of the tile; for each locus, the
// tiles are always added in the same order (i0 = 0, TILE_I, ...).
//
// BLOCK (with NEWTON): the Hessian sums are register-blocked as on the Q
// side (see HESS_BLOCK_MIN_K), with the roles of the individuals and the
// loci swapped: H[(k, k2), j] = sum_i d2_ij q_ki q_k2i, with the products
// q_ki q_k2i of the tile formed once (qq, TILE_I x npair_pad), the d2 of the
// tile and the loci of the block stored by individual (d2t, TILE_I x TILE_I,
// one row per individual, one column per locus), and the sums of the block
// in hp (n_pair x TILE_I, one column per locus), all in pB. The tiles are
// still added in order, so the results do not depend on the number of
// threads.
template <bool NEWTON, bool GL, bool BLOCK>
RADMIX_MULTIVERSION
static void accum_p_block(const GenoView& geno, const double* qT,
                          const double* p, int K, size_t j0, size_t j1,
                          double* out1, double* out2, Slice* s,
                          double* pB) {
  const size_t I = geno.nrow(), nj = j1 - j0, KK = (size_t) K * K;
  // out2 holds a K x K matrix per locus with NEWTON, K values otherwise
  const size_t n2 = NEWTON ? KK : K;
  // the number of output slots used, including that of a padding locus
  const size_t n_slots = nj + nj % 2;
  // the sums are accumulated over the tiles of individuals, so they start
  // at zero
  std::fill(out1, out1 + K * n_slots, 0.0);
  std::fill(out2, out2 + n2 * n_slots, 0.0);
  // the coefficients of the two loci of each step (s[0] for locus j, s[1]
  // for locus j + 1); a and b are used by EM, d1 and d2 by Newton
  const double* __restrict a1 = s[0].a.data();
  const double* __restrict b1 = s[0].b.data();
  const double* __restrict a2 = s[1].a.data();
  const double* __restrict b2 = s[1].b.data();
  const double* __restrict d11 = s[0].d1.data();
  const double* __restrict d21 = s[0].d2.data();
  const double* __restrict d12 = s[1].d1.data();
  const double* __restrict d22 = s[1].d2.data();
  // scratch for the products d2 q_ki of the two loci
  double* __restrict wq1 = s[0].tmp.data();
  double* __restrict wq2 = s[1].tmp.data();
  // the buffers of the blocked Hessians (see above); the columns of d2t
  // beyond the last locus stay at zero
  const size_t n_pair = (size_t) K * (K + 1) / 2;
  const size_t npair_pad = hess_pairs_pad(n_pair);
  double* qq = pB;
  double* d2t = BLOCK ? qq + TILE_I * npair_pad : nullptr;
  double* hp = BLOCK ? d2t + TILE_I * TILE_I : nullptr;
  if (BLOCK) {
    std::fill(d2t, d2t + TILE_I * TILE_I, 0.0);
    std::fill(hp, hp + n_pair * TILE_I, 0.0);
  }
  for (size_t i0 = 0; i0 < I; i0 += TILE_I) {
    const size_t i1 = std::min(I, i0 + TILE_I), ni = i1 - i0;
    if (BLOCK) {
      // the products of the ancestry proportions of each individual of the
      // tile, reused by every locus of the block
      for (size_t ii = 0; ii < ni; ii++) {
        double* c = qq + ii * npair_pad;
        size_t pr = 0;
        for (int k = 0; k < K; k++) {
          const double qk = qT[I * k + i0 + ii];
          for (int k2 = k; k2 < K; k2++) c[pr++] = qk * qT[I * k2 + i0 + ii];
        }
        for (; pr < npair_pad; pr++) c[pr] = 0;
      }
    }
    // loci j0 + jj and j0 + jj + 1 (the second is the padding locus, with
    // zero coefficients, when jj + 1 == nj)
    for (size_t jj = 0; jj < nj; jj += 2) {
      const size_t j = j0 + jj;
      slice_coef<NEWTON, GL>(geno, qT, p + K * j, K, j, i0, i1, s[0]);
      if (jj + 1 < nj)
        slice_coef<NEWTON, GL>(geno, qT, p + K * (j + 1), K, j + 1, i0, i1,
                               s[1]);
      else
        s[1].zero_coef(ni);
      // the output columns of the two loci (those of the second follow
      // those of the first)
      double* o1 = out1 + K * jj;
      double* o2 = out2 + n2 * jj;
      for (int k = 0; k < K; k++) {
        const double* __restrict qk = qT + I * k + i0;
        if (NEWTON) {
          // gradient: sum_i d1 q_ki
          double sg1 = 0.0, sg2 = 0.0;
#pragma omp simd reduction(+:sg1, sg2)
          for (size_t ii = 0; ii < ni; ii++) {
            sg1 += d11[ii] * qk[ii];
            sg2 += d12[ii] * qk[ii];
          }
          o1[k] += sg1;
          o1[K + k] += sg2;
          // (blocked: the Hessians are summed after the loci of the block)
          if (BLOCK) continue;
          // Hessian column k, rows k2 >= k: sum_i (d2 q_ki) q_k2i; the
          // product in brackets depends only on k, so it is formed once,
          // and each entry is then a single dot product
#pragma omp simd
          for (size_t ii = 0; ii < ni; ii++) {
            wq1[ii] = d21[ii] * qk[ii];
            wq2[ii] = d22[ii] * qk[ii];
          }
          double* H1 = o2 + K * k;
          double* H2 = o2 + KK + K * k;
          for (int k2 = k; k2 < K; k2++) {
            const double* __restrict qk2 = qT + I * k2 + i0;
            double sh1 = 0.0, sh2 = 0.0;
#pragma omp simd reduction(+:sh1, sh2)
            for (size_t ii = 0; ii < ni; ii++) {
              sh1 += wq1[ii] * qk2[ii];
              sh2 += wq2[ii] * qk2[ii];
            }
            H1[k2] += sh1;
            H2[k2] += sh2;
          }
        } else {
          // EM: sum_i a q_ki and sum_i b q_ki
          double sa1 = 0.0, sb1 = 0.0, sa2 = 0.0, sb2 = 0.0;
#pragma omp simd reduction(+:sa1, sb1, sa2, sb2)
          for (size_t ii = 0; ii < ni; ii++) {
            sa1 += a1[ii] * qk[ii];
            sb1 += b1[ii] * qk[ii];
            sa2 += a2[ii] * qk[ii];
            sb2 += b2[ii] * qk[ii];
          }
          o1[k] += sa1;
          o2[k] += sb1;
          o1[K + k] += sa2;
          o2[K + k] += sb2;
        }
      }
      // (blocked: store the d2 of the two loci, one column each)
      if (NEWTON && BLOCK) {
        for (size_t ii = 0; ii < ni; ii++) d2t[ii * TILE_I + jj] = d21[ii];
        if (jj + 1 < nj)
          for (size_t ii = 0; ii < ni; ii++)
            d2t[ii * TILE_I + jj + 1] = d22[ii];
      }
    }
    // (blocked: the Hessian sums of the tile, over its individuals, for
    // every locus of the block)
    if (NEWTON && BLOCK) hess_block(d2t, qq, (int) ni, n_pair, npair_pad,
                                    nj, hp);
  }
  // (blocked: the sums into the K x K matrices of the loci)
  if (NEWTON && BLOCK) {
    for (size_t jj = 0; jj < nj; jj++) {
      double* H = out2 + KK * jj;
      size_t idx = 0;
      for (int k = 0; k < K; k++)
        for (int k2 = k; k2 < K; k2++, idx++)
          H[k2 + K * k] = hp[TILE_I * idx + jj];
    }
  }
  // fill in the other half of each Hessian, once all the tiles of
  // individuals have been added
  if (NEWTON)
    for (size_t jj = 0; jj < nj; jj++) mirror_hess(out2 + KK * jj, K);
}

// The P-side driver: the sums of accum_p_block() for every locus,
// parallelised over blocks of TILE_J loci, then finish(j0, j1, out1, out2)
// for each block of loci [j0, j1), with the sums of the block (out1 K x nj,
// out2 K x K x nj with NEWTON and K x nj without, nj = j1 - j0). As in
// q_side(), finish() runs on the thread that owns the block, must not call
// R, and receives thread-local buffers valid only during the call. qT is
// the transposed q (I x K), p is K x J. Each block covers all the
// individuals, so no partial sums need combining across threads.
template <bool NEWTON, class Finish>
static void p_side(const GenoView& geno, const double* qT, const double* p,
                   int K, int n_cores, Finish finish) {
  const size_t J = geno.ncol(), KK = (size_t) K * K;
  const ptrdiff_t nb = (ptrdiff_t) n_blocks(J, TILE_J);
  // the accumulation function for the kind of data (see slice_gl())
  // (with the blocked Hessians from HESS_BLOCK_MIN_K)
  const bool block = NEWTON && K >= HESS_BLOCK_MIN_K;
  const auto accum = geno.is_gl() ?
    (block ? accum_p_block<NEWTON, true, true> :
     accum_p_block<NEWTON, true, false>) :
    (block ? accum_p_block<NEWTON, false, true> :
     accum_p_block<NEWTON, false, false>);
  const size_t n_pair = (size_t) K * (K + 1) / 2;
  const size_t p_scratch = block ?
    TILE_I * (hess_pairs_pad(n_pair) + TILE_I) + n_pair * TILE_I : 0;
  if (geno.stream()) {
    // streaming (see "Streaming" at the top): the blocks chunk by chunk,
    // with the chunks in order or in reverse order; each block is
    // independent, so the order does not change the results
    const bool backward = geno.direction(false);
    const size_t cb = chunk_blocks(geno), n_chunks = n_blocks(nb, cb);
    Prefetcher prefetcher(geno);
    for (size_t cc = 0; cc < n_chunks; cc++) {
      const size_t c = backward ? n_chunks - 1 - cc : cc;
      const ptrdiff_t b0 = c * cb, b1 = std::min((size_t) nb, (c + 1) * cb);
      if (cc + 1 < n_chunks) {
        const size_t cn = backward ? c - 1 : c + 1;
        prefetcher.start(cn * cb * TILE_J,
                         std::min(J, (cn + 1) * cb * TILE_J));
      }
#pragma omp parallel num_threads(n_cores)
{
      Slice s[2];
      std::vector<double> out1(K * TILE_J), out2((NEWTON ? KK : K) * TILE_J);
      std::vector<double> pB(p_scratch);
#pragma omp for schedule(RADMIX_OMP_SCHEDULE)
      for (ptrdiff_t x = 0; x < b1 - b0; x++) {
        // (forwards within the chunk, see GenoView::direction())
        const ptrdiff_t bl = b0 + x;
        const size_t j0 = bl * TILE_J, j1 = std::min(J, j0 + TILE_J);
        accum(geno, qT, p, K, j0, j1, out1.data(), out2.data(), s,
              pB.data());
        finish(j0, j1, out1.data(), out2.data());
      }
}
    }
    return;
  }
#pragma omp parallel num_threads(n_cores)
{
  // thread-local buffers: two slices, and the sums of one block of loci
  // (TILE_J slots, so that the spare slot of a padding locus always exists)
  Slice s[2];
  std::vector<double> out1(K * TILE_J), out2((NEWTON ? KK : K) * TILE_J);
  std::vector<double> pB(p_scratch);
#pragma omp for schedule(RADMIX_OMP_SCHEDULE)
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    const size_t j0 = bl * TILE_J, j1 = std::min(J, j0 + TILE_J);
    accum(geno, qT, p, K, j0, j1, out1.data(), out2.data(), s, pB.data());
    finish(j0, j1, out1.data(), out2.data());
  }
}
}

// Sums over all the genotypes, for loglik() and deviance(): NS running sums,
// to which slice_sum(ni, s, acc) adds the contribution of each slice (after
// slice_geno(), which fills s.qp, s.g and s.h). Parallelised over blocks of
// TILE_J loci; each block keeps its own sums, which are added up in block
// order at the end, so the totals do not depend on the number of threads.
// q (K x I) and p (K x J) are in the usual layout; total receives the NS
// sums. slice_sum runs inside the parallel region, so it must not call R.
// Within a block, each locus is processed for all the tiles of individuals
// in turn, so the bytes of a column of the backing file are read in order.
template <int NS, bool GL, class SliceSum>
static void sum_slices(const GenoView& geno, const double* q, const double* p,
                       int K, int n_cores, double* total,
                       SliceSum slice_sum) {
  const size_t I = geno.nrow(), J = geno.ncol();
  std::vector<double> qT = transpose(q, K, I);
  const ptrdiff_t nb = (ptrdiff_t) n_blocks(J, TILE_J);
  std::vector<double> partial(nb * NS, 0.0);
  // the sums of block bl into partial
  auto block = [&](ptrdiff_t bl, Slice& s) {
    const size_t j1 = std::min(J, (bl + 1) * TILE_J);
    double acc[NS] = {};
    for (size_t j = bl * TILE_J; j < j1; j++) {
      for (size_t i0 = 0; i0 < I; i0 += TILE_I) {
        const size_t i1 = std::min(I, i0 + TILE_I);
        slice_geno<GL>(geno, qT.data(), p + K * j, K, j, i0, i1, s);
        slice_sum(i1 - i0, s, acc);
      }
    }
    for (int c = 0; c < NS; c++) partial[NS * bl + c] = acc[c];
  };
  if (geno.stream()) {
    // streaming (see "Streaming" at the top): the blocks chunk by chunk,
    // with the chunks in order or in reverse order; each keeps its own sums,
    // so the order does not change the results
    const bool backward = geno.direction(false);
    const size_t cb = chunk_blocks(geno), n_chunks = n_blocks(nb, cb);
    Prefetcher prefetcher(geno);
    for (size_t cc = 0; cc < n_chunks; cc++) {
      const size_t c = backward ? n_chunks - 1 - cc : cc;
      const ptrdiff_t b0 = c * cb, b1 = std::min((size_t) nb, (c + 1) * cb);
      if (cc + 1 < n_chunks) {
        const size_t cn = backward ? c - 1 : c + 1;
        prefetcher.start(cn * cb * TILE_J,
                         std::min(J, (cn + 1) * cb * TILE_J));
      }
#pragma omp parallel num_threads(n_cores)
{
      Slice s;
#pragma omp for schedule(RADMIX_OMP_SCHEDULE)
      for (ptrdiff_t x = 0; x < b1 - b0; x++)
        block(b0 + x, s);
}
    }
  } else {
#pragma omp parallel num_threads(n_cores)
{
  Slice s;
  // (dynamic scheduling balances the load if some threads are slower, e.g.
  // while the backing file is being paged in)
#pragma omp for schedule(RADMIX_OMP_SCHEDULE)
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    const size_t j1 = std::min(J, (bl + 1) * TILE_J);
    double acc[NS] = {};
    for (size_t j = bl * TILE_J; j < j1; j++) {
      for (size_t i0 = 0; i0 < I; i0 += TILE_I) {
        const size_t i1 = std::min(I, i0 + TILE_I);
        slice_geno<GL>(geno, qT.data(), p + K * j, K, j, i0, i1, s);
        slice_sum(i1 - i0, s, acc);
      }
    }
    for (int c = 0; c < NS; c++) partial[NS * bl + c] = acc[c];
  }
}
  }
  // the totals, adding up the blocks in order (outside the parallel region)
  for (int c = 0; c < NS; c++) {
    total[c] = 0.0;
    for (ptrdiff_t bl = 0; bl < nb; bl++) total[c] += partial[NS * bl + c];
  }
}

// The log-likelihood (see radmix_kernels.h): the sum of loglik_slice() (or
// of loglik_gl_slice(), for genotype likelihoods) over all the slices. (The
// kind of data is tested once, outside the pass: see slice_gl().)
double loglik(const GenoView& geno, const double* q, const double* p, int K,
              int n_cores) {
  double ll;
  if (geno.is_gl()) {
    sum_slices<1, true>(geno, q, p, K, n_cores, &ll,
                        [](size_t ni, Slice& s, double* acc) {
                          acc[0] += loglik_gl_slice(ni, s);
                        });
  } else {
    sum_slices<1, false>(geno, q, p, K, n_cores, &ll,
                         [](size_t ni, Slice& s, double* acc) {
                           acc[0] += loglik_slice(ni, s);
                         });
  }
  return ll;
}

// The binomial deviance of the visible genotypes, the prediction error of
// cross-validation (Alexander and Lange 2011). A held-out genotype g of an
// individual with n allele copies is predicted by mu = n qp, the expected
// number of alternate alleles under the estimates fitted without it, and
// the error is
//   d(g, mu) = g log(g / mu) + (n - g) log((n - g) / (n - mu))
//            = g log(g / (n qp)) + h log(h / (n (1 - qp))),
// with 0 log 0 = 0: 0 for a perfect prediction, growing as mu moves away
// from g. (ADMIXTURE gives this formula for diploids, n = 2.) With a
// held-out view (radmix_geno_access.h), only the genotypes of the scored
// fold are visible. *dev receives the sum of the deviances and *count the
// number of visible genotypes, so that the error can be averaged over all
// the held-out genotypes of all the folds.
//
// With genotype likelihoods, the genotypes are not observed, and the error
// of an individual with data at a locus is
//   d = log max_x f(x) - log f(qp),
//   f(x) = (1 - x)^2 + 2 e1 x (1 - x) + e2 x^2,
// the log-likelihood ratio of the best frequency of allele2 for that
// individual and locus alone (the "saturated" model, as g / n for a
// genotype) against the prediction qp (see deviance_gl()). For the
// likelihoods of a hard call, this is the deviance of the genotype above:
// f is (1 - x)^2, 2 x (1 - x) or x^2 (up to a constant factor), whose
// maxima are at x = 0, 1/2 and 1. *count is then the number of visible
// individuals and loci with data.
static void deviance_gl(const GenoView& geno, const double* q,
                        const double* p, int K, int n_cores, double* total);

void deviance(const GenoView& geno, const double* q, const double* p, int K,
              int n_cores, double* dev, double* count) {
  // total[0]: the sum of the deviances, total[1]: the number of genotypes
  double total[2];
  if (geno.is_gl()) {
    deviance_gl(geno, q, p, K, n_cores, total);
    *dev = total[0];
    *count = total[1];
    return;
  }
  sum_slices<2, false>(geno, q, p, K, n_cores, total,
                       [](size_t ni, Slice& s, double* acc) {
    // (a plain scalar loop: cross-validation computes the deviance once per
    // fold, a single pass that costs little next to the fit of the fold)
    for (size_t ii = 0; ii < ni; ii++) {
      const double g = s.g[ii], h = s.h[ii], qp = s.qp[ii];
      // a visible genotype has g + h = n_i > 0, a hidden or missing one 0
      if (g + h <= 0) continue;
      const double n = g + h;
      // (the estimates of a fit are projected into [1e-5, 1 - 1e-5], so qp
      // is not exactly 0 or 1 here; the tests on g and h implement
      // 0 log 0 = 0)
      if (g > 0) acc[0] += g * std::log(g / (n * qp));
      if (h > 0) acc[0] += h * std::log(h / (n * (1.0 - qp)));
      acc[1] += 1.0;
    }
  });
  *dev = total[0];
  *count = total[1];
}

// The deviance of genotype likelihoods (see deviance()). f is a quadratic
// in x, a x^2 + b x + 1 with a = 1 - 2 e1 + e2 and b = 2 e1 - 2, so its
// maximum on [0, 1] is at an end (f(0) = 1, f(1) = e2) or, if f is concave
// (a < 0) with its vertex inside (0, 1), at the vertex x* = -b / (2 a),
// where f(x*) = 1 - b^2 / (4 a) (a sum of positive terms, as a < 0). f > 0
// on [0, 1], as the likelihoods of a source are positive.
static void deviance_gl(const GenoView& geno, const double* q,
                        const double* p, int K, int n_cores, double* total) {
  sum_slices<2, true>(geno, q, p, K, n_cores, total,
                      [](size_t ni, Slice& s, double* acc) {
    // (a plain scalar loop, as for genotypes)
    for (size_t ii = 0; ii < ni; ii++) {
      // an individual without data, or hidden by the mask
      if (!(s.obs[ii] > 0)) continue;
      const double e1 = s.e1[ii], e2 = s.e2[ii], x = s.qp[ii];
      const double y = 1.0 - x;
      const double a = 1.0 - 2.0 * e1 + e2, b = 2.0 * e1 - 2.0;
      double f_max = std::max(1.0, e2);
      if (a < 0) {
        const double xs = -b / (2.0 * a);
        if (xs > 0 && xs < 1) f_max = std::max(f_max, 1.0 - b * b / (4.0 * a));
      }
      const double f_x = y * y + 2.0 * e1 * x * y + e2 * x * x;
      acc[0] += std::log(f_max / f_x);
      acc[1] += 1.0;
    }
  });
}

// One EM step (the algorithm of FRAPPE, Tang et al. 2005), updating both
// blocks from the same (q, p).
//
// Each allele copy of individual i at locus j comes from one of the K
// ancestral populations. Given an alternate allele, the probability that it
// came from population k is q_ki p_kj / qp; given a reference allele, it is
// q_ki (1 - p_kj) / (1 - qp). The E step therefore assigns to population k
//   g q_ki p_kj / qp + h q_ki (1 - p_kj) / (1 - qp)
//     = q_ki [a p_kj + b (1 - p_kj)]
// expected allele copies, of which q_ki p_kj a are alternate alleles. The M
// step sets q_ki to the share of the allele copies of individual i assigned
// to k, and p_kj to the share of alternate alleles among the copies at locus
// j assigned to k. Summed over k, the copies assigned at one locus are
// a qp + b (1 - qp) = g + h = n_i, so q_ki is normalised by n_i J_i, where
// J_i is the number of loci observed in individual i; normalising by the sum
// over k gives the same result, and keeps q on the simplex exactly.
//
// The Q block takes a pass over blocks of individuals, the P block a pass
// over blocks of loci, both from the same (q, p) (and the same transposed
// q), so one EM step costs two passes over the genotypes.
void em(const GenoView& geno, const double* q, const double* p, int K,
        double* q_new, double* p_new, int n_cores,
        const unsigned char* fixed) {
  const size_t I = geno.nrow();
  std::vector<double> qT = transpose(q, K, I);
  // With fixed individuals, their columns are copied into q_new, and the Q
  // side only processes the free ones (free_idx, with the compact transpose
  // of their columns, qT_free); the P side below uses every individual.
  std::vector<size_t> free_idx;
  std::vector<double> qT_free;
  if (fixed != nullptr) {
    free_idx = free_individuals(fixed, q, K, I, q_new);
    qT_free = transpose_subset(q, K, free_idx.data(), free_idx.size());
  }
  // Q: q_new_ki proportional to q_ki sum_j [a p_kj + b (1 - p_kj)]
  // (f indexes the processed individuals, i the individuals of the data;
  // they are the same without fixed individuals)
  auto q_finish = [&](size_t i0, size_t i1, const double* sums,
                      const double*) {
    for (size_t f = i0; f < i1; f++) {
      const size_t i = fixed == nullptr ? f : free_idx[f];
      const double* si = sums + K * (f - i0);
      const double* qi = q + K * i;
      double* qn = q_new + K * i;
      double denom = 0.0;
      for (int k = 0; k < K; k++) {
        qn[k] = si[k] * qi[k];
        denom += qn[k];
      }
      // denom is 0 only for an individual without any observed genotype:
      // there is no information, so q is left unchanged
      if (denom > 0)
        for (int k = 0; k < K; k++) qn[k] /= denom;
      else
        for (int k = 0; k < K; k++) qn[k] = qi[k];
    }
  };
  if (fixed == nullptr) {
    q_side<false, false>(geno, qT.data(), p, K, n_cores, NULL, q_finish);
  } else if (!free_idx.empty()) {
    q_side<false, false>(geno, qT_free.data(), p, K, n_cores, NULL, q_finish,
                         free_idx.data(), free_idx.size());
  }
  // P: p_new = A / (A + B), where A = p sum_i a q and B = (1 - p) sum_i b q
  // are the expected numbers of alternate and reference allele copies at
  // locus j from population k
  p_side<false>(geno, qT.data(), p, K, n_cores,
                [&](size_t j0, size_t j1, const double* sa,
                    const double* sb) {
    for (size_t j = j0; j < j1; j++) {
      for (int k = 0; k < K; k++) {
        // x indexes the K x J matrices, y the K x nj sums of the block
        const size_t x = k + K * j, y = k + K * (j - j0);
        const double A = p[x] * sa[y], B = (1.0 - p[x]) * sb[y];
        // A + B = 0 means that no allele copy is assigned to k at this locus
        // (e.g. all its genotypes are missing): the frequency is kept. The
        // result is not projected, so a locus where all the observed
        // alleles are reference alleles gets p = 0 exactly (and one where
        // they are all alternate alleles p = 1); the other passes handle qp
        // of exactly 0 or 1 (slice_coef(), loglik_slice()).
        p_new[x] = (A + B > 0) ? A / (A + B) : p[x];
      }
    }
  });
}

// The gradients and negative Hessians of the columns of q (see
// radmix_kernels.h): the sums of q_side(), copied block by block into grad
// (K x I) and hess (K x K x I). Only used by the tests (radmix_grad_q_cpp());
// the fit computes them inside update_q(), without storing them.
void grad_q(const GenoView& geno, const double* q, const double* p, int K,
            double* grad, double* hess, int n_cores) {
  const size_t KK = (size_t) K * K;
  std::vector<double> qT = transpose(q, K, geno.nrow());
  q_side<true, false>(geno, qT.data(), p, K, n_cores, NULL,
                      [&](size_t i0, size_t i1, const double* g,
                          const double* H) {
    std::copy(g, g + K * (i1 - i0), grad + K * i0);
    std::copy(H, H + KK * (i1 - i0), hess + KK * i0);
  });
}

// The same for the columns of p: the sums of p_side(), copied into grad
// (K x J) and hess (K x K x J). Only used by the tests (radmix_grad_p_cpp()).
void grad_p(const GenoView& geno, const double* q, const double* p, int K,
            double* grad, double* hess, int n_cores) {
  const size_t KK = (size_t) K * K;
  std::vector<double> qT = transpose(q, K, geno.nrow());
  p_side<true>(geno, qT.data(), p, K, n_cores,
               [&](size_t j0, size_t j1, const double* g, const double* H) {
    std::copy(g, g + K * (j1 - j0), grad + K * j0);
    std::copy(H, H + KK * (j1 - j0), hess + KK * j0);
  });
}

// The constrained Newton update of one column x (K values) of q (simplex =
// true) or p (simplex = false), from the gradient grad and the negative
// Hessian H of the log-likelihood with respect to that column: x_new is the
// projection of x + delta onto the feasible set, where the step delta
// maximises the quadratic approximation grad' delta - delta' H delta / 2
// subject to 0 <= x + delta <= 1 (and sum_k delta_k = 1 - sum_k x_k for q),
// solved exactly by qp_solve_column() (radmix_qp_sweep.cpp). The projection
// (radmix_projections.cpp) removes rounding errors and keeps every value
// at least 1e-5 away from 0 and 1, so that qp stays strictly inside (0, 1).
//
// A step that is not finite is discarded (the column keeps its value for
// this update, up to the projection). This happens when the gradient or the
// Hessian is not finite, i.e. when a genotype is impossible under the
// current estimates: qp of exactly 0 with g > 0, or of exactly 1 with h > 0
// (see slice_coef()). The projections keep every estimate that has been
// through a Newton or quasi-Newton step away from 0 and 1, so this is a
// safeguard for values that have not (the unprojected output of em(), or
// any q and p passed to the exported functions): the Inf or NaN would
// otherwise spread to the whole column, and through the next updates to the
// log-likelihood. The column gets another chance at the next update, once
// the other block has moved. (The return value of the solver is not used:
// delta is the best step it found in any case. A Hessian that is not
// numerically positive definite is shifted towards the identity until it
// is, see qp_solve_column_ridge(), rather than giving no step: at a large
// K, a nearly empty cluster otherwise froze whole columns, and the fit
// could stall and report convergence.)
//
// H (K x K), grad and x (K values) are read only; x_new (K values) must not
// overlap them. work provides the scratch of the solver and the buffer of
// delta, so each thread needs its own.
static inline void newton_column(const double* H, const double* grad,
                                 const double* x, int K, bool simplex,
                                 double* x_new, QPWork& work) {
  double* delta = work.delta.data();
  qp_solve_column_ridge(H, grad, x, K, simplex, delta, work);
  // discard a step with any non-finite component (see above)
  for (int k = 0; k < K; k++) {
    if (!std::isfinite(delta[k])) {
      std::fill(delta, delta + K, 0.0);
      break;
    }
  }
  // take the step, then project onto {sum = 1, q >= 1e-5} or
  // [1e-5, 1 - 1e-5]
  for (int k = 0; k < K; k++) x_new[k] = x[k] + delta[k];
  if (simplex) project_q(x_new, K);
  else project_p(x_new, K);
}

// The block update of q (see radmix_kernels.h). Each block of individuals
// solves the Newton steps of its columns as soon as their sums are
// complete, inside q_side()'s parallel region, so the gradients and
// Hessians are used while they are in cache, and those of all the
// individuals are never stored (except with segments, which are only used
// with few individuals). Every block gets its own QPWork, a few small
// vectors whose allocation costs nothing next to the accumulation. With
// LOGLIK, *ll also receives L(q, p) (see q_side()); the update is the same
// either way. With a mask of fixed individuals (a supervised analysis),
// their columns are copied, and neither their sums nor their quadratic
// programs are computed. With LOGLIK, a second Q-side pass over the fixed
// individuals alone then computes the log-likelihood of their genotypes
// (q_side() without NEWTON), which is added to that of the free ones: each
// genotype is still read once, and the update and the log-likelihood cost
// what a pass of update_q_loglik() without fixed individuals costs, less
// the sums and quadratic programs of the fixed individuals.
template <bool LOGLIK>
static void update_q_impl(const GenoView& geno, const double* q,
                          const double* p, int K, double* q_new,
                          int n_cores, double* ll,
                          const unsigned char* fixed) {
  const size_t KK = (size_t) K * K, I = geno.nrow();
  if (fixed == nullptr) {
    std::vector<double> qT = transpose(q, K, I);
    q_side<true, LOGLIK>(geno, qT.data(), p, K, n_cores, ll,
                         [&](size_t i0, size_t i1, const double* g,
                             const double* H) {
      QPWork work(K);
      for (size_t i = i0; i < i1; i++)
        newton_column(H + KK * (i - i0), g + K * (i - i0), q + K * i, K,
                      true, q_new + K * i, work);
    });
    return;
  }
  // with fixed individuals: copy their columns, and update the free ones
  // only (positions f in the list of free individuals, individuals
  // free_idx[f] of the data), from the compact transpose of their columns
  std::vector<size_t> free_idx = free_individuals(fixed, q, K, I, q_new);
  // (the log-likelihoods of the genotypes of the free and of the fixed
  // individuals, each independent of the number of threads)
  double ll_free = 0.0, ll_fixed = 0.0;
  if (!free_idx.empty()) {
    std::vector<double> qT =
      transpose_subset(q, K, free_idx.data(), free_idx.size());
    q_side<true, LOGLIK>(geno, qT.data(), p, K, n_cores, &ll_free,
                         [&](size_t i0, size_t i1, const double* g,
                             const double* H) {
      QPWork work(K);
      for (size_t f = i0; f < i1; f++) {
        const size_t i = free_idx[f];
        newton_column(H + KK * (f - i0), g + K * (f - i0), q + K * i, K,
                      true, q_new + K * i, work);
      }
    }, free_idx.data(), free_idx.size());
  }
  if (!LOGLIK) return;
  // the log-likelihood of the genotypes of the fixed individuals, from the
  // compact transpose of their columns
  std::vector<size_t> fixed_idx;
  for (size_t i = 0; i < I; i++)
    if (fixed[i]) fixed_idx.push_back(i);
  if (!fixed_idx.empty()) {
    std::vector<double> qT =
      transpose_subset(q, K, fixed_idx.data(), fixed_idx.size());
    q_side<false, true>(geno, qT.data(), p, K, n_cores, &ll_fixed,
                        [](size_t, size_t, const double*, const double*) {},
                        fixed_idx.data(), fixed_idx.size());
  }
  *ll = ll_free + ll_fixed;
}

void update_q(const GenoView& geno, const double* q, const double* p, int K,
              double* q_new, int n_cores, const unsigned char* fixed) {
  update_q_impl<false>(geno, q, p, K, q_new, n_cores, NULL, fixed);
}

// The block update of q and the log-likelihood of (q, p), in one pass (see
// radmix_kernels.h).
double update_q_loglik(const GenoView& geno, const double* q,
                       const double* p, int K, double* q_new, int n_cores,
                       const unsigned char* fixed) {
  double ll;
  update_q_impl<true>(geno, q, p, K, q_new, n_cores, &ll, fixed);
  return ll;
}

// The block update of p: the same for the loci, block by block, inside
// p_side()'s parallel region, so the Hessians of all the loci (K times the
// size of p) are never stored.
void update_p(const GenoView& geno, const double* q, const double* p, int K,
              double* p_new, int n_cores) {
  const size_t KK = (size_t) K * K;
  std::vector<double> qT = transpose(q, K, geno.nrow());
  p_side<true>(geno, qT.data(), p, K, n_cores,
               [&](size_t j0, size_t j1, const double* g, const double* H) {
    QPWork work(K);
    for (size_t j = j0; j < j1; j++)
      newton_column(H + KK * (j - j0), g + K * (j - j0), p + K * j, K, false,
                    p_new + K * j, work);
  });
}

}  // namespace radmix

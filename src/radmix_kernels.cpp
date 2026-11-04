// The passes over the genotypes of gt_radmix(): the log-likelihood, the
// deviance of held-out genotypes, the EM step, and the gradients, Hessians
// and Newton updates of the two blocks of parameters. See radmix_kernels.h
// for the model, the notation and the layout of q and p, and
// dev_notes/admixture_implementation.md for the algorithm as a whole.
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

#include <algorithm>
#include <cmath>
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
// and EM steps 15-30% faster. The log-likelihood and the deviance are left
// alone: their AVX2 versions are slower (they call log() from the C
// library inside their loops). The two versions can differ by rounding (FMA
// and wider vectors change the order of operations), so results are
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

// Tile sizes. A vector of TILE_I accumulators takes 512 bytes, so the
// Q-side accumulators (K + K (K + 1) / 2 vectors) and the TILE_I x K corner
// of q stay in the L1 or L2 cache for moderate K. loglik_slice() needs
// TILE_I == LOG_BATCH^2, and the padding of odd ranges of loci on the P side
// needs TILE_J to be even (see accum_p_block()).
static const size_t TILE_I = 64;
static const size_t TILE_J = 64;
static_assert(TILE_J % 2 == 0, "TILE_J must be even");

// the number of blocks of size `tile` needed to cover n (the last one may be
// shorter; 0 for n = 0)
static inline size_t n_blocks(size_t n, size_t tile) {
  return (n + tile - 1) / tile;
}

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

// Thread-local buffers for one slice (one locus x up to TILE_I individuals):
// qp, the decoded genotypes g and h, and the coefficients defined at the top
// of the file, a and b (EM) or d1 and d2 (Newton). loglik_slice() uses a and
// b as scratch, and accum_p_block() uses tmp. Each thread allocates its
// slices once, at the start of a parallel region, and reuses them for every
// slice it processes; every buffer is recomputed for each slice. Entry ii of
// each buffer refers to individual i0 + ii of the slice; only the first
// ni = i1 - i0 entries are meaningful.
struct Slice {
  std::vector<double> qp, g, h, a, b, d1, d2, tmp;
  Slice() : qp(TILE_I), g(TILE_I), h(TILE_I), a(TILE_I), b(TILE_I),
    d1(TILE_I), d2(TILE_I), tmp(TILE_I) {}
  // zero coefficients for the first ni individuals: the padding locus of a
  // range with an odd number of loci, which adds nothing to any sum
  void zero_coef(size_t ni) {
    std::fill(a.begin(), a.begin() + ni, 0.0);
    std::fill(b.begin(), b.begin() + ni, 0.0);
    std::fill(d1.begin(), d1.begin() + ni, 0.0);
    std::fill(d2.begin(), d2.begin() + ni, 0.0);
  }
};

// qp and the decoded genotypes of the slice of locus j and individuals
// [i0, i1): s.qp, s.g and s.h for ii = i - i0 in [0, i1 - i0). qT is the
// transposed q (I x K, see transpose()), pj the K frequencies of locus j
// (p + K j); i1 - i0 must be at most TILE_I. A missing (or hidden) genotype
// gets g = h = 0, so that it drops out of every sum without a test.
static inline void slice_geno(const GenoView& geno, const double* qT,
                              const double* pj, int K, size_t j, size_t i0,
                              size_t i1, Slice& s) {
  const size_t I = geno.nrow(), ni = i1 - i0;
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
  // Decoding the genotypes needs table lookups, which do not vectorise, so
  // it has loops of its own. Missing and hidden genotypes are zeroed by
  // multiplying by 0 rather than by a branch, which would be mispredicted
  // when they are scattered at random (as the folds of cross-validation
  // are).
  // Fast path: when the individuals of the slice are adjacent rows of the
  // backing file and are decoded alike (no pseudohaploids), their bytes are
  // read directly and decoded with two tables giving g and h (both 0 for a
  // missing genotype). This gives the same values as the general loop below.
  const double* gtab;
  const double* htab;
  const unsigned char* __restrict bytes =
    geno.slice_bytes(i0, ni, j, gtab, htab);
  if (bytes != nullptr) {
    if (!geno.cv_masked()) {
      for (size_t ii = 0; ii < ni; ii++) {
        const unsigned char b = bytes[ii];
        g[ii] = gtab[b];
        h[ii] = htab[b];
      }
    } else {
      for (size_t ii = 0; ii < ni; ii++) {
        const unsigned char b = bytes[ii];
        const double vis = geno.cv_visible(i0 + ii, j) ? 1.0 : 0.0;
        g[ii] = vis * gtab[b];
        h[ii] = vis * htab[b];
      }
    }
    return;
  }
  // general path: any order of the rows, pseudohaploids, masks
  double gij;
  for (size_t ii = 0; ii < ni; ii++) {
    const double obs = geno.get(i0 + ii, j, gij) ? 1.0 : 0.0;
    g[ii] = obs * gij;
    h[ii] = obs * (geno.n_trials(i0 + ii) - gij);
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
template <bool NEWTON>
static inline void slice_coef(const GenoView& geno, const double* qT,
                              const double* pj, int K, size_t j, size_t i0,
                              size_t i1, Slice& s) {
  slice_geno(geno, qT, pj, K, j, i0, i1, s);
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
}

// Q side: the sums over loci [j0, j1) of each individual of the block
// [i0, i1) (all the loci, unless the work is split into segments, see
// q_side()):
// - NEWTON: grad_k = sum_j d1 p_kj (grad, K x ni) and
//   hess_kk2 = sum_j d2 p_kj p_k2j (hess, K x K x ni), the gradient and the
//   negative Hessian of the log-likelihood with respect to q_.i;
// - EM: grad_k = sum_j a p_kj + b (1 - p_kj) (grad, K x ni; hess unused),
//   the expected number of allele copies of individual i from population
//   k, divided by q_ki (see em()).
// Outputs are overwritten. qT is the transposed q (I x K), p is K x J, and
// i1 - i0 must be at most TILE_I. s holds two slices (one per locus of a
// step); gT (K vectors of TILE_I values) and hT (K (K + 1) / 2 vectors of
// TILE_I values, only used if NEWTON) are scratch for the accumulators.
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
template <bool NEWTON>
RADMIX_MULTIVERSION
static void accum_q_block(const GenoView& geno, const double* qT,
                          const double* p, int K, size_t i0, size_t i1,
                          size_t j0, size_t j1, double* grad, double* hess,
                          Slice* s, double* gT, double* hT) {
  const size_t ni = i1 - i0, KK = (size_t) K * K;
  // the number of entries (k, k2) with k2 >= k of a symmetric K x K matrix
  const size_t n_pair = (size_t) K * (K + 1) / 2;
  // (the whole vectors are zeroed, including the tail beyond ni, which is
  // never read)
  std::fill(gT, gT + K * TILE_I, 0.0);
  if (NEWTON) std::fill(hT, hT + n_pair * TILE_I, 0.0);
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
  for (size_t j = j0; j < j1; j += 2) {
    // the frequencies of the two loci; the padding locus of an odd range
    // (j + 1 == j1) has zero coefficients, so any frequencies will do
    const bool two = j + 1 < j1;
    const double* pj = p + K * j;
    const double* pj2 = two ? pj + K : pj;
    slice_coef<NEWTON>(geno, qT, pj, K, j, i0, i1, s[0]);
    if (two) slice_coef<NEWTON>(geno, qT, pj2, K, j + 1, i0, i1, s[1]);
    else s[1].zero_coef(ni);
    if (NEWTON) {
      // gradient: gT[, k] += d1 p_kj, for both loci
      for (int k = 0; k < K; k++) {
        const double pk = pj[k], pk2 = pj2[k];
        double* __restrict gk = gT + TILE_I * k;
#pragma omp simd
        for (size_t ii = 0; ii < ni; ii++)
          gk[ii] += d11[ii] * pk + d12[ii] * pk2;
      }
      // Hessian: hT[, (k, k2)] += d2 p_kj p_k2j for k2 >= k; the products
      // of the frequencies do not depend on the individual, so they are
      // formed outside the loop over individuals
      double* __restrict hk = hT;
      for (int k = 0; k < K; k++) {
        for (int k2 = k; k2 < K; k2++, hk += TILE_I) {
          const double c1 = pj[k] * pj[k2], c2 = pj2[k] * pj2[k2];
#pragma omp simd
          for (size_t ii = 0; ii < ni; ii++)
            hk[ii] += d21[ii] * c1 + d22[ii] * c2;
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
  // the lower triangle, and then mirrored
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
// I < MIN_Q_TASKS TILE_I = 2048, and n_seg = ceil(MIN_Q_TASKS / n_blocks), so
// n_seg I stays below 2 MIN_Q_TASKS TILE_I = 4096: under 4 MB for K = 10.

// the target number of Q-side tasks, enough to balance a few dozen threads
static const size_t MIN_Q_TASKS = 32;
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
// With one segment, each block is accumulated over all the loci by one
// thread, into thread-local buffers. With several, each (segment, block)
// task writes its partial sums into its own slot of grad_part and
// hess_part, and the partial sums of each individual are then added up in
// segment order before finish(). Both stages run in a single parallel
// region, separated by the barrier at the end of the first loop.
template <bool NEWTON, class Finish>
static void q_side(const GenoView& geno, const double* qT, const double* p,
                   int K, int n_cores, Finish finish) {
  const size_t I = geno.nrow(), J = geno.ncol(), KK = (size_t) K * K;
  const size_t n_pair = (size_t) K * (K + 1) / 2;
  const size_t nbi = n_blocks(I, TILE_I), n_seg = q_segments(I, J);
  // the partial sums of individual i in segment seg start at
  // K (seg I + i) (and KK (seg I + i) for the Hessians)
  std::vector<double> grad_part(n_seg > 1 ? n_seg * K * I : 0);
  std::vector<double> hess_part(n_seg > 1 && NEWTON ? n_seg * KK * I : 0);
#pragma omp parallel num_threads(n_cores)
{
  // thread-local buffers: two slices, the accumulators, and the sums of one
  // block of individuals
  Slice s[2];
  std::vector<double> gT(K * TILE_I), hT(NEWTON ? n_pair * TILE_I : 0);
  std::vector<double> grad(K * TILE_I), hess(NEWTON ? KK * TILE_I : 0);
  if (n_seg > 1) {
    // tasks t = seg * nbi + bl: blocks vary fastest, so threads that start
    // together read the same columns of the backing file
    const ptrdiff_t n_tasks = (ptrdiff_t) (nbi * n_seg);
#pragma omp for schedule(dynamic)
    for (ptrdiff_t t = 0; t < n_tasks; t++) {
      const size_t seg = t / nbi, bl = t % nbi;
      const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
      // n_seg contiguous, nearly equal ranges of loci
      const size_t j0 = J * seg / n_seg, j1 = J * (seg + 1) / n_seg;
      accum_q_block<NEWTON>(geno, qT, p, K, i0, i1, j0, j1,
                            &grad_part[K * (seg * I + i0)],
                            NEWTON ? &hess_part[KK * (seg * I + i0)] : NULL,
                            s, gT.data(), NEWTON ? hT.data() : NULL);
    }
    // (the implicit barrier at the end of the loop above guarantees that
    // all the partial sums are written before they are added up below)
  }
  // one task per block of individuals: accumulate it over all the loci, or
  // add up its partial sums, then hand the complete sums to finish()
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < (ptrdiff_t) nbi; bl++) {
    const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
    if (n_seg == 1) {
      accum_q_block<NEWTON>(geno, qT, p, K, i0, i1, 0, J, grad.data(),
                            NEWTON ? hess.data() : NULL, s, gT.data(),
                            NEWTON ? hT.data() : NULL);
    } else {
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
    finish(i0, i1, grad.data(), hess.data());
  }
}
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
template <bool NEWTON>
RADMIX_MULTIVERSION
static void accum_p_block(const GenoView& geno, const double* qT,
                          const double* p, int K, size_t j0, size_t j1,
                          double* out1, double* out2, Slice* s) {
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
  for (size_t i0 = 0; i0 < I; i0 += TILE_I) {
    const size_t i1 = std::min(I, i0 + TILE_I), ni = i1 - i0;
    // loci j0 + jj and j0 + jj + 1 (the second is the padding locus, with
    // zero coefficients, when jj + 1 == nj)
    for (size_t jj = 0; jj < nj; jj += 2) {
      const size_t j = j0 + jj;
      slice_coef<NEWTON>(geno, qT, p + K * j, K, j, i0, i1, s[0]);
      if (jj + 1 < nj)
        slice_coef<NEWTON>(geno, qT, p + K * (j + 1), K, j + 1, i0, i1, s[1]);
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
#pragma omp parallel num_threads(n_cores)
{
  // thread-local buffers: two slices, and the sums of one block of loci
  // (TILE_J slots, so that the spare slot of a padding locus always exists)
  Slice s[2];
  std::vector<double> out1(K * TILE_J), out2((NEWTON ? KK : K) * TILE_J);
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    const size_t j0 = bl * TILE_J, j1 = std::min(J, j0 + TILE_J);
    accum_p_block<NEWTON>(geno, qT, p, K, j0, j1, out1.data(), out2.data(),
                          s);
    finish(j0, j1, out1.data(), out2.data());
  }
}
}

// x^e for e in {0, 1, 2}, without branches: the Lagrange interpolation
// polynomial in e through (0, 1), (1, x) and (2, x^2). It gives 1 for e = 0
// even if x = 0; for other values of e the result is meaningless, but
// loglik_slice() never uses it then.
static inline double pow012(double x, double e) {
  return 0.5 * (e - 1.0) * (e - 2.0) + x * e * (2.0 - e) +
    0.5 * x * x * e * (e - 1.0);
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
// (all TILE_I entries of s.a, as pass 3 pads them). ni is at most TILE_I.
static double loglik_slice(size_t ni, Slice& s) {
  const double* __restrict qp = s.qp.data();
  const double* __restrict g = s.g.data();
  const double* __restrict h = s.h.data();
  double* __restrict z = s.a.data();
  double* __restrict c = s.b.data();
  // Pass 1 (vectorised): z = qp^g (1 - qp)^h, assuming that g and h are 0,
  // 1 or 2, and a flag c that is 0 if and only if they are
#pragma omp simd
  for (size_t ii = 0; ii < ni; ii++) {
    z[ii] = pow012(qp[ii], g[ii]) * pow012(1.0 - qp[ii], h[ii]);
    c[ii] = std::fabs(g[ii] * (g[ii] - 1.0) * (g[ii] - 2.0)) +
      std::fabs(h[ii] * (h[ii] - 1.0) * (h[ii] - 2.0));
  }
  // Pass 2 (scalar, rarely taken): the genotypes that cannot go into the
  // products get explicit logarithms, and a neutral z = 1. This includes
  // z = 0 (an impossible genotype, qp of exactly 0 with g > 0, or 1 with
  // h > 0), which correctly gives -Inf. (| rather than || avoids a second
  // branch.)
  double ll = 0.0;
  for (size_t ii = 0; ii < ni; ii++) {
    if ((z[ii] < LOG_MIN_Z) | (c[ii] != 0)) {
      z[ii] = 1.0;
      if (g[ii] > 0) ll += g[ii] * std::log(qp[ii]);
      if (h[ii] > 0) ll += h[ii] * std::log(1.0 - qp[ii]);
    }
  }
  // Pass 3: LOG_BATCH interleaved products (product b collects z[b],
  // z[b + 8], z[b + 16], ...), so that the inner loop multiplies 8
  // independent accumulators by 8 contiguous values and vectorises. A short
  // slice is padded with ones (log(1) = 0) to keep the loop bounds fixed.
  std::fill(z + ni, z + TILE_I, 1.0);
  double prod[LOG_BATCH];
  for (int b = 0; b < LOG_BATCH; b++) prod[b] = 1.0;
  for (int m = 0; m < LOG_BATCH; m++)
#pragma omp simd
    for (int b = 0; b < LOG_BATCH; b++) prod[b] *= z[LOG_BATCH * m + b];
  for (int b = 0; b < LOG_BATCH; b++) ll += std::log(prod[b]);
  return ll;
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
template <int NS, class SliceSum>
static void sum_slices(const GenoView& geno, const double* q, const double* p,
                       int K, int n_cores, double* total,
                       SliceSum slice_sum) {
  const size_t I = geno.nrow(), J = geno.ncol();
  std::vector<double> qT = transpose(q, K, I);
  const ptrdiff_t nb = (ptrdiff_t) n_blocks(J, TILE_J);
  std::vector<double> partial(nb * NS, 0.0);
#pragma omp parallel num_threads(n_cores)
{
  Slice s;
  // (dynamic scheduling balances the load if some threads are slower, e.g.
  // while the backing file is being paged in)
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    const size_t j1 = std::min(J, (bl + 1) * TILE_J);
    double acc[NS] = {};
    for (size_t j = bl * TILE_J; j < j1; j++) {
      for (size_t i0 = 0; i0 < I; i0 += TILE_I) {
        const size_t i1 = std::min(I, i0 + TILE_I);
        slice_geno(geno, qT.data(), p + K * j, K, j, i0, i1, s);
        slice_sum(i1 - i0, s, acc);
      }
    }
    for (int c = 0; c < NS; c++) partial[NS * bl + c] = acc[c];
  }
}
  // the totals, adding up the blocks in order (outside the parallel region)
  for (int c = 0; c < NS; c++) {
    total[c] = 0.0;
    for (ptrdiff_t bl = 0; bl < nb; bl++) total[c] += partial[NS * bl + c];
  }
}

// The log-likelihood (see radmix_kernels.h): the sum of loglik_slice() over
// all the slices.
double loglik(const GenoView& geno, const double* q, const double* p, int K,
              int n_cores) {
  double ll;
  sum_slices<1>(geno, q, p, K, n_cores, &ll,
                [](size_t ni, Slice& s, double* acc) {
                  acc[0] += loglik_slice(ni, s);
                });
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
void deviance(const GenoView& geno, const double* q, const double* p, int K,
              int n_cores, double* dev, double* count) {
  // total[0]: the sum of the deviances, total[1]: the number of genotypes
  double total[2];
  sum_slices<2>(geno, q, p, K, n_cores, total,
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
        double* q_new, double* p_new, int n_cores) {
  std::vector<double> qT = transpose(q, K, geno.nrow());
  // Q: q_new_ki proportional to q_ki sum_j [a p_kj + b (1 - p_kj)]
  q_side<false>(geno, qT.data(), p, K, n_cores,
                [&](size_t i0, size_t i1, const double* sums,
                    const double*) {
    for (size_t i = i0; i < i1; i++) {
      const double* si = sums + K * (i - i0);
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
  });
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
  q_side<true>(geno, qT.data(), p, K, n_cores,
               [&](size_t i0, size_t i1, const double* g, const double* H) {
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
// delta is the best step it found in any case, 0 if it failed.)
//
// H (K x K), grad and x (K values) are read only; x_new (K values) must not
// overlap them. work provides the scratch of the solver and the buffer of
// delta, so each thread needs its own.
static inline void newton_column(const double* H, const double* grad,
                                 const double* x, int K, bool simplex,
                                 double* x_new, QPWork& work) {
  double* delta = work.delta.data();
  qp_solve_column(H, grad, x, K, simplex, delta, work);
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
// vectors whose allocation costs nothing next to the accumulation.
void update_q(const GenoView& geno, const double* q, const double* p, int K,
              double* q_new, int n_cores) {
  const size_t KK = (size_t) K * K;
  std::vector<double> qT = transpose(q, K, geno.nrow());
  q_side<true>(geno, qT.data(), p, K, n_cores,
               [&](size_t i0, size_t i1, const double* g, const double* H) {
    QPWork work(K);
    for (size_t i = i0; i < i1; i++)
      newton_column(H + KK * (i - i0), g + K * (i - i0), q + K * i, K, true,
                    q_new + K * i, work);
  });
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

// Port of the CPU kernels of OpenADMIXTURE.jl (src/loops.jl and the inner
// algorithms in src/algorithms_inner.jl). See radmix_kernels.h for
// conventions.
//
// The model
// ---------
// For individual i and locus j, let g_ij be the observed count of the alt
// allele and n_i the number of allele copies (2 for diploids, 1 for
// pseudohaploids). ADMIXTURE assumes
//   g_ij ~ Binomial(n_i, qp_ij),   qp_ij = sum_k q_ki p_kj,
// where q_ki is the fraction of the genome of individual i that comes from
// ancestral population k (each column of q lies on the simplex) and p_kj is
// the alt allele frequency of locus j in population k (each entry in [0, 1]).
// Up to a constant, the log-likelihood is
//   L = sum_{i,j observed} g_ij log(qp_ij) + (n_i - g_ij) log(1 - qp_ij).
// Missing genotypes simply do not contribute.
//
// Throughout this file we write, for one (i, j) pair,
//   g = g_ij, h = n_i - g_ij, qp = qp_ij,
//   a = g / qp,          b = h / (1 - qp),
//   wa = g / qp^2,       wb = h / (1 - qp)^2.
// a and b are the derivatives of g log(x) and h log(y) at x = qp and
// y = 1 - qp; wa and wb are minus the corresponding second derivatives.
// Every quantity below is a sum of these coefficients weighted by q or p.
//
// What each kernel computes (and its Julia counterpart)
// ----------------------------------------------------
// - loglik: L itself (loglikelihood() / loglikelihood_loop()).
// - em: one FRAPPE EM step (em_q!() / em_q_loop!() and
//   em_p!() / em_p_loop!()):
//     q_new_ki = q_ki sum_j [a p_kj + b (1 - p_kj)] / (n_i J_i)
//     p_new_kj = A / (A + B), A = p_kj sum_i a q_ki,
//                             B = (1 - p_kj) sum_i b q_ki.
// - grad_q: per individual, gradient and negative Hessian of L with
//   respect to q_.i (update_q_loop!()):
//     grad_k      = sum_j a p_kj + b (1 - p_kj)
//     hess_{k,k2} = sum_j wa p_kj p_k2j + wb (1 - p_kj)(1 - p_k2j)
// - grad_p: per locus, gradient and negative Hessian of L with respect
//   to p_.j (update_p_loop!()):
//     grad_k      = sum_i (a - b) q_ki
//     hess_{k,k2} = sum_i (wa + wb) q_ki q_k2i
// - update_q / update_p: the Newton step of update_q!() /
//   update_p!(): the gradient and Hessian above feed a small quadratic
//   program per column (radmix_qp_sweep.h), whose solution is added to q
//   or p and projected back onto the feasible set
//   (radmix_projections.h).
//
// Derivation of the Q-side formulas
// ---------------------------------
// The plain derivative of L with respect to q_ki is
//   sum_j g p_kj / qp - h p_kj / (1 - qp) = sum_j (a - b) p_kj.
// Julia instead uses the identity 1 - qp = sum_k q_k (1 - p_k), valid on the
// simplex, to write the likelihood in the homogeneous form
//   g log(sum_k q_k p_k) + h log(sum_k q_k (1 - p_k)),
// whose derivative is a p_k + b (1 - p_k). The two gradients differ by
// sum_j b, which is the same for every k. The Newton step is constrained to
// sum_k delta_k = 0 (it must stay on the simplex), and on that subspace a
// constant shift of the gradient adds nothing to grad' delta, so the step is
// unaffected. Likewise, for such delta, delta' (1 - p) = -delta' p, so the
// two forms of the Hessian agree on the subspace too. We follow Julia so
// that results match.
//
// Data layout
// -----------
// q is K x I and p is K x J, column-major, so the K values for one
// individual (or one locus) are contiguous. The FBM is individuals x loci,
// column-major, so the genotypes of consecutive individuals at one locus are
// adjacent bytes: loops that run over individuals innermost read memory
// sequentially. GenoView (radmix_geno_access.h) turns each raw byte into
// a dosage (or "missing") and handles subsets and reorderings of rows and
// columns.
//
// Parallel design
// ---------------
// Work is organised in tiles of TILE_I individuals x TILE_J loci. Q-side work
// (sums over loci for each individual) is parallelised over blocks of TILE_I
// individuals, P-side work (sums over individuals for each locus) over blocks
// of TILE_J loci. Each thread owns the output columns of its blocks, so there
// are no races and no cross-thread reductions, and every column is summed in
// the same order whatever the number of threads. With few individuals there
// are too few Q-side blocks to keep the threads busy, so the loci are also
// split into segments, whose partial sums are added up in a fixed order
// (two-dimensional blocking, see accum_q_segmented()).
//
// The innermost unit of work is a slice of one locus x up to TILE_I
// individuals (individuals are contiguous in the column-major FBM). For each
// slice, qp is computed once into a thread-local buffer with unit-stride loops
// over a transposed copy of q, the genotypes are decoded, and the likelihood
// coefficients a = g / qp, b = (n - g) / (1 - qp) (and wa = a / qp,
// wb = b / (1 - qp) for Hessians) are computed without branches, so that all
// these loops can be vectorised. The accumulation loops then run over the
// slice:
// - Q side: accumulators are stored as K (and K (K + 1) / 2 for the upper
//   triangle of the Hessian) vectors of length TILE_I, updated with
//   unit-stride vectorisable loops for every locus;
// - P side: each accumulator gets a (SIMD) dot product over the slice; the
//   TILE_I x K corner of q stays in L1 while the TILE_J loci are processed.
//
// Determinism: block boundaries depend only on I, J and the tile sizes, never
// on the number of threads, and each block is processed start to finish by a
// single thread in a fixed order. The only scalar reduction (the
// log-likelihood) is stored as one partial sum per block and the partial sums
// are added up serially. `#pragma omp simd reduction` may reassociate a sum,
// but in a way fixed at compile time (by the vector width), so the result is
// still bit-identical for any n_cores. schedule(dynamic) only changes which
// thread handles a block, not how the block is computed.
//
// Thread safety: everything that the parallel regions read (GenoView, q, p,
// the transposed copy of q) is immutable, everything they write is either
// thread-local scratch (Slice, gT/hT, QPWork) or a disjoint set of output
// columns. No R API (allocation, Rcpp::stop, R_CheckUserInterrupt, ...) is
// called inside a parallel region: R is not thread safe and an R error would
// longjmp out of the OpenMP region. All validation happens beforehand, in
// the Rcpp wrappers at the bottom of this file.
// (and in the GenoView constructor, which every wrapper, and the fitting
// loop in radmix_fit.cpp, runs before calling a kernel)
//
// Deviations from OpenADMIXTURE.jl
// --------------------------------
// The kernels follow the Julia code, except for the following, each
// documented where it is implemented:
// - the binomial number of trials is n_i (2 or 1, for pseudohaploids)
//   rather than the constant 2;
// - the EM update of q is normalised by the number of observed loci of
//   each individual rather than by J, which keeps q on the simplex with
//   missing data (em());
// - a Newton step that is not finite is discarded, keeping the current
//   column for that sweep (discard_nonfinite_step());
// - and, in the solver and projections called by the update kernels, mu is
//   floored at 0, vanishing pivots are never swept (radmix_qp_sweep.cpp)
//   and non-finite values are tolerated (radmix_projections.cpp).

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "radmix_geno_access.h"
#include "radmix_kernels.h"
#include "radmix_qp_sweep.h"
#include "radmix_projections.h"
using namespace Rcpp;

// internal code of gt_radmix() (see the note on the namespace in
// radmix_geno_access.h)
namespace radmix {


// Tile sizes. 64 individuals x 8 bytes = 512 bytes per K-vector of
// accumulators, which keeps the Q-side accumulators (K + K (K + 1) / 2
// vectors) and the TILE_I x K corner of q in L1/L2 for moderate K.
// loglik_slice() relies on TILE_I == LOG_BATCH^2 (see below).
static const size_t TILE_I = 64;
static const size_t TILE_J = 64;

// number of blocks of size `tile` needed to cover n
static inline size_t n_blocks(size_t n, size_t tile) {
  return (n + tile - 1) / tile;
}

// copy the upper triangle (k2 >= k, stored in column k) to the lower one
// The accumulators only compute the K (K + 1) / 2 entries with k2 >= k (the
// Hessians are symmetric, so the rest is redundant work). These are stored at
// H[k2 + K * k], i.e. row k2 >= column k; this loop fills in the symmetric
// entries H[k + K * k2]. (In row/column terms the computed half is the lower
// triangle of the column-major matrix; "upper" above refers to the index
// pairs (k, k2) with k2 >= k. Either way the result is the full matrix.)
static inline void mirror_hess(double* H, int K) {
  for (int k = 0; k < K; k++)
    for (int k2 = k + 1; k2 < K; k2++) H[k + K * k2] = H[k2 + K * k];
}

// 1 - p for every entry of p (length n). Computed once per kernel call and
// shared read-only by all threads, so that the Q-side inner loops do not
// recompute it for every individual block.
static std::vector<double> one_minus(const double* p, size_t n) {
  std::vector<double> omp(n);
  for (size_t x = 0; x < n; x++) omp[x] = 1.0 - p[x];
  return omp;
}

// q is K x I; returns its transpose (I x K, column-major)
// In qT, the values of q_k. for consecutive individuals are contiguous
// (qT + I * k + i0 is q_k,i0, q_k,i0+1, ...). The slice loops run over
// individuals innermost, so this gives them unit stride (and hence SIMD
// loads) instead of a stride of K. The copy also means the kernels never
// read R memory from inside a parallel region except through q and p, which
// are read only.
static std::vector<double> transpose(const double* q, int K, size_t I) {
  std::vector<double> qT((size_t) K * I);
  for (size_t i = 0; i < I; i++)
    for (int k = 0; k < K; k++) qT[i + I * k] = q[k + K * i];
  return qT;
}

// Thread-local buffers for one slice (one locus x up to TILE_I individuals)
// qp: sum_k q_ki p_kj; g: alt allele count; h: n_i - g; a, b, wa, wb: the
// likelihood coefficients defined at the top of the file; tmp: scratch.
// Each thread constructs one Slice at the start of a parallel region and
// reuses it for every slice it processes, so there is no allocation in the
// inner loops and no sharing between threads. Some functions reuse a and b
// as scratch (loglik_slice() stores z and c in them).
// accum_p_block() also overwrites a and wa with a - b and wa + wb, and
// uses tmp for the products (wa + wb) q_k. Every buffer is recomputed for
// each slice, so none of this outlives the slice.
struct Slice {
  std::vector<double> qp, g, h, a, b, wa, wb, tmp;
  Slice() : qp(TILE_I), g(TILE_I), h(TILE_I), a(TILE_I), b(TILE_I),
    wa(TILE_I), wb(TILE_I), tmp(TILE_I) {}
};

// qp and decoded genotypes for locus j and individuals [i0, i1). Missing
// genotypes get g = h = 0 (h = n - g), so they drop out of every sum.
// Inputs: qT (I x K, transposed q), pj (the K frequencies of locus j).
// Outputs: s.qp[ii], s.g[ii], s.h[ii] for ii = i - i0 in [0, i1 - i0).
// This plays the role of qp_block!() in loops.jl, but for a single locus.
static inline void slice_geno(const GenoView& geno, const double* qT,
                              const double* pj, int K, size_t j, size_t i0,
                              size_t i1, Slice& s) {
  const size_t I = geno.nrow(), ni = i1 - i0;
  // __restrict promises the compiler that these buffers do not alias, which
  // it needs in order to vectorise the loops below
  double* __restrict qp = s.qp.data();
  double* __restrict g = s.g.data();
  double* __restrict h = s.h.data();
  // qp = sum_k qT[, k] * p_kj, accumulated one population at a time: each
  // pass is an axpy over contiguous memory, which vectorises cleanly
  // (looping over k innermost would be a short reduction with stride I)
  std::fill(qp, qp + ni, 0.0);
  for (int k = 0; k < K; k++) {
    const double pk = pj[k];
    const double* __restrict qk = qT + I * k + i0;
#pragma omp simd
    for (size_t ii = 0; ii < ni; ii++) qp[ii] += qk[ii] * pk;
  }
  // Decode the genotypes. This loop has a table lookup, so it does not
  // vectorise; it is kept separate so that everything else can. Setting
  // both g and h to 0 for a missing (or hidden) genotype means it
  // contributes 0 to every coefficient (a, b, wa, wb) and to the
  // log-likelihood, with no further tests downstream. They are multiplied by
  // 1 or 0 rather than set in a branch, which would be mispredicted when
  // missing or hidden genotypes are scattered at random (as the folds of
  // cross-validation are); the values are the same (gij and n - gij are
  // finite and non-negative, so 0 times them is +0).
  // The bytes read here are those of individuals i0..i1-1 in column j of
  // the FBM: adjacent bytes when the gen_tibble keeps the rows of the
  // backing file in order, scattered within the column otherwise.
  double gij;
  for (size_t ii = 0; ii < ni; ii++) {
    const double obs = geno.get(i0 + ii, j, gij) ? 1.0 : 0.0;
    g[ii] = obs * gij;
    h[ii] = obs * (geno.n_trials(i0 + ii) - gij);
  }
}

// Likelihood coefficients for the slice. A term whose coefficient (g or
// n - g) is zero is set to 0 (dividing by 1 instead of qp or 1 - qp), so that
// qp of exactly 0 or 1 does not produce NaN.
// HESS selects whether the second-order weights wa and wb are also needed
// (Hessians) or only a and b (EM step, gradient-only uses). Being a template
// parameter, the test disappears at compile time.
// Outputs: s.a, s.b and, if HESS, s.wa, s.wb (plus everything slice_geno sets).
template <bool HESS>
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
  double* __restrict wa = s.wa.data();
  double* __restrict wb = s.wb.data();
#pragma omp simd
  for (size_t ii = 0; ii < ni; ii++) {
    // x = (g > 0) ? qp : 1 and y = (h > 0) ? 1 - qp : 1, written with
    // arithmetic masks so that the loop has no control flow
    // Why the guard matters: a homozygous reference genotype (g = 0) at a
    // locus where qp = 0 would otherwise give 0 / 0 = NaN, although its
    // contribution g / qp is genuinely 0 (likewise h = 0 with qp = 1, and
    // missing genotypes, where g = h = 0). If g > 0 and qp = 0 the
    // likelihood really is -Inf and a becomes Inf, as it should.
    const double mg = (double) (g[ii] > 0), mh = (double) (h[ii] > 0);
    const double x = mg * qp[ii] + (1.0 - mg);
    const double y = mh * (1.0 - qp[ii]) + (1.0 - mh);
    // first derivatives: d/dqp of g log(qp) and of h log(1 - qp) (up to sign)
    a[ii] = g[ii] / x;
    b[ii] = h[ii] / y;
    if (HESS) {
      // minus the second derivatives: g / qp^2 and h / (1 - qp)^2
      wa[ii] = a[ii] / x;
      wb[ii] = b[ii] / y;
    }
  }
}

// Q side: for individuals [i0, i1), sum over loci [j0, j1) (all loci, unless
// the work is split into segments of loci, see accum_q_segmented())
//   grad_k = sum_j a p_kj + b (1 - p_kj)
//   hess_{k,k2} = sum_j wa p_kj p_k2j + wb (1 - p_kj)(1 - p_k2j)   (if HESS)
// grad is K x ni and hess K x K x ni (ni = i1 - i0); both are overwritten.
// gT and hT are scratch accumulators of size K * TILE_I and
// K (K + 1) / 2 * TILE_I. omp = 1 - p.
// (The (1 - p) form comes from the homogeneous likelihood explained at the
// top of the file; Julia: update_q_loop!() and, without the Hessian and
// before multiplying by q, em_q_loop!().)
//
// The block of individuals is fixed and we stream over all J loci. The
// accumulators are kept "transposed" (one vector of length TILE_I per k, or
// per (k, k2) pair) so that each update is a unit-stride loop over the
// individuals of the slice, vectorised across individuals. They are only
// copied into the K-major output layout once, at the end.
template <bool HESS>
static void accum_q_block(const GenoView& geno, const double* qT,
                          const double* p, const double* omp, int K,
                          size_t i0, size_t i1, size_t j0, size_t j1,
                          double* grad, double* hess, Slice& s, double* gT,
                          double* hT) {
  const size_t ni = i1 - i0, KK = (size_t) K * K;
  // number of distinct entries of a symmetric K x K matrix
  const size_t n_pair = (size_t) K * (K + 1) / 2;
  // Accumulator vector k starts at gT + TILE_I * k (and pair idx at
  // hT + TILE_I * idx) whatever ni is: a short last block simply leaves
  // the tail of each vector unused.
  std::fill(gT, gT + K * TILE_I, 0.0);
  if (HESS) std::fill(hT, hT + n_pair * TILE_I, 0.0);
  const double* __restrict a = s.a.data();
  const double* __restrict b = s.b.data();
  const double* __restrict wa = s.wa.data();
  const double* __restrict wb = s.wb.data();
  // loci in order j0..j1-1: the summation order for each individual is fixed
  for (size_t j = j0; j < j1; j++) {
    const double* pj = p + K * j;
    const double* oj = omp + K * j;
    // qp, genotypes and coefficients for this locus and block of individuals
    slice_coef<HESS>(geno, qT, pj, K, j, i0, i1, s);
    // gradient: gT[, k] += a * p_kj + b * (1 - p_kj), elementwise over ii
    for (int k = 0; k < K; k++) {
      const double pk = pj[k], ok = oj[k];
      double* __restrict gk = gT + TILE_I * k;
#pragma omp simd
      for (size_t ii = 0; ii < ni; ii++) gk[ii] += a[ii] * pk + b[ii] * ok;
    }
    if (HESS) {
      // Hessian: only pairs k2 >= k, in the packed order
      // (0,0), (0,1), ..., (0,K-1), (1,1), ..., (K-1,K-1);
      // hk walks through the packed accumulators one TILE_I vector at a time.
      // The products p_kj p_k2j and (1 - p_kj)(1 - p_k2j) do not depend on
      // the individual, so they are hoisted out of the inner loop.
      double* __restrict hk = hT;
      for (int k = 0; k < K; k++) {
        for (int k2 = k; k2 < K; k2++, hk += TILE_I) {
          const double c1 = pj[k] * pj[k2], c2 = oj[k] * oj[k2];
#pragma omp simd
          for (size_t ii = 0; ii < ni; ii++)
            hk[ii] += wa[ii] * c1 + wb[ii] * c2;
        }
      }
    }
  }
  // Scatter the accumulators into the output layout: grad K x ni and hess
  // K x K x ni, column-major, one column (or K x K matrix) per individual.
  // Only the first ni entries of each accumulator vector are meaningful
  // (the last block may be shorter than TILE_I).
  for (size_t ii = 0; ii < ni; ii++) {
    for (int k = 0; k < K; k++) grad[k + K * ii] = gT[ii + TILE_I * k];
    if (HESS) {
      double* H = hess + KK * ii;
      // same packed order as in the accumulation loop above
      size_t idx = 0;
      for (int k = 0; k < K; k++)
        for (int k2 = k; k2 < K; k2++, idx++)
          H[k2 + K * k] = hT[ii + TILE_I * idx];
      mirror_hess(H, K);
    }
  }
}

// ---------------------------------------------------------------------------
// Two-dimensional blocking of the Q side
//
// Q-side sums run over loci for each individual, and are naturally split into
// blocks of TILE_I individuals, each owned by one thread. But with few
// individuals there are few blocks: 379 individuals give 6 blocks, so 4
// threads are poorly balanced, and 50 individuals give a single block, so the
// Q side would run on one thread however many are available. Population
// genetic datasets often have tens to hundreds of individuals and many loci,
// so this matters.
//
// When there are fewer than MIN_Q_TASKS blocks of individuals, we therefore
// also split the loci into segments, and process each (block of individuals,
// segment of loci) pair as a separate task. Each task accumulates partial
// sums into its own buffer (no races), and the partial sums of the segments
// are then added up in segment order. The number of segments depends only
// on I and J, never on the number of threads, so results are identical for
// any n_cores (they differ from a single pass over the loci only by
// rounding, as the sums are grouped differently).
//
// Memory: the partial sums take n_seg * I * (K + K^2) doubles (the Hessians
// are stored in full, K x K, ready to be added up). Segments are only used
// when there are fewer than MIN_Q_TASKS blocks, i.e. I <= 1984, and n_seg is
// ceil(MIN_Q_TASKS / n_blocks), so n_seg * I stays below
// 2 * MIN_Q_TASKS * TILE_I = 4096 (the worst case is I = 1984: 31 blocks,
// 2 segments): under 4 MB for K = 10. Larger datasets keep the single-pass
// path (with block-local buffers only), unchanged.
//
// This is an implementation detail with no counterpart in OpenADMIXTURE.jl,
// which splits the Q side over individuals only.

// target number of Q-side tasks; enough to balance up to a few dozen threads
static const size_t MIN_Q_TASKS = 32;
// minimum number of loci per segment, so that each task has enough work to
// amortise its overheads
static const size_t MIN_SEG_LOCI = 256;

// Number of segments of loci for the Q side (1 = no segmentation). Depends
// only on the dimensions of the data, so that results do not depend on the
// number of threads.
static size_t q_segments(size_t I, size_t J) {
  // number of blocks of individuals; with enough of them, the single-pass
  // path already has enough independent tasks
  const size_t nbi = n_blocks(I, TILE_I);
  // (no individuals: nothing to split, and no division by zero below)
  if (nbi == 0 || nbi >= MIN_Q_TASKS) return 1;
  // otherwise, enough segments to have about MIN_Q_TASKS tasks in total
  // (rounded up: e.g. 6 blocks need 6 segments, giving 36 tasks)
  size_t n_seg = (MIN_Q_TASKS + nbi - 1) / nbi;
  // do not make segments shorter than MIN_SEG_LOCI loci
  n_seg = std::min(n_seg, std::max((size_t) 1, J / MIN_SEG_LOCI));
  return n_seg;
}

// Q-side sums for all individuals, split into (block, segment) tasks:
// grad (K x I) and, if HESS, hess (K x K x I) are overwritten with the sums
// over all loci (see accum_q_block() for the formulas).
template <bool HESS>
static void accum_q_segmented(const GenoView& geno, const double* qT,
                              const double* p, const double* omp, int K,
                              size_t n_seg, double* grad, double* hess,
                              int n_cores) {
  const size_t I = geno.nrow(), J = geno.ncol(), KK = (size_t) K * K;
  // tasks are the (segment, block of individuals) pairs, numbered
  // t = seg * nbi + bl
  // (blocks vary fastest, so threads that start at the same time work on
  // the same segment, i.e. on the same columns of the backing file)
  const size_t nbi = n_blocks(I, TILE_I);
  const ptrdiff_t n_tasks = (ptrdiff_t) (nbi * n_seg);
  // partial sums of each segment: K x I (and K x K x I) per segment, stored
  // one segment after the other (segment seg starts at seg * K * I), so the
  // partial sums of individual i in segment seg are at K * (seg * I + i)
  std::vector<double> grad_part(n_seg * K * I);
  std::vector<double> hess_part(HESS ? n_seg * KK * I : 0);
#pragma omp parallel num_threads(n_cores)
{
  Slice s;
  std::vector<double> gT(K * TILE_I), hT(HESS ? K * (K + 1) / 2 * TILE_I : 0);
  // tasks in any order (dynamic scheduling): each writes only its own
  // (segment, block) part of the partial sums
#pragma omp for schedule(dynamic)
  for (ptrdiff_t t = 0; t < n_tasks; t++) {
    // which segment of loci and which block of individuals this task is
    const size_t seg = t / nbi, bl = t % nbi;
    const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
    // segment boundaries: n_seg contiguous, nearly equal ranges of loci
    const size_t j0 = J * seg / n_seg, j1 = J * (seg + 1) / n_seg;
    // the same accumulation as the single-pass path, restricted to loci
    // [j0, j1), writing the partial sums of individuals [i0, i1) for this
    // segment
    accum_q_block<HESS>(geno, qT, p, omp, K, i0, i1, j0, j1,
                        &grad_part[K * (seg * I + i0)],
                        HESS ? &hess_part[KK * (seg * I + i0)] : NULL, s,
                        gT.data(), HESS ? hT.data() : NULL);
  }
  // add up the segments, always in segment order (deterministic); each
  // thread handles whole individuals
  // (the implicit barrier at the end of the task loop above guarantees that
  // every partial sum has been written before any thread reads it here)
#pragma omp for schedule(static)
  for (ptrdiff_t i = 0; i < (ptrdiff_t) I; i++) {
    for (int k = 0; k < K; k++) {
      double sum = 0.0;
      for (size_t seg = 0; seg < n_seg; seg++)
        sum += grad_part[K * (seg * I + i) + k];
      grad[K * i + k] = sum;
    }
    if (HESS) {
      for (size_t kk = 0; kk < KK; kk++) {
        double sum = 0.0;
        for (size_t seg = 0; seg < n_seg; seg++)
          sum += hess_part[KK * (seg * I + i) + kk];
        hess[KK * i + kk] = sum;
      }
    }
  }
}
}

// P side: for loci [j0, j1), sum over all individuals. Two modes:
// - EM = false: out1 = grad (K x nj), grad_k = sum_i (a - b) q_ki;
//               out2 = hess (K x K x nj),
//               hess_{k,k2} = sum_i (wa + wb) q_ki q_k2i
// - EM = true:  out1 = sum_i a q_ki and out2 = sum_i b q_ki (both K x nj)
// Outputs are overwritten.
// (Julia: update_p_loop!() and em_p_loop!(). The gradient is the plain
// derivative of L with respect to p_kj: d qp / d p_kj = q_ki, so
// d/dp_kj [g log(qp) + h log(1 - qp)] = (a - b) q_ki, and minus the second
// derivative is (wa + wb) q_ki q_k2i.)
//
// Loop order: individuals in tiles of TILE_I (outer), then the nj loci of
// the block (inner). The TILE_I x K corner of qT used by one tile is reused
// for all nj loci, so it stays in cache; each output column is updated once
// per tile of individuals, always in the same order (i0 = 0, TILE_I, ...).
template <bool EM>
static void accum_p_block(const GenoView& geno, const double* qT,
                          const double* p, int K, size_t j0, size_t j1,
                          double* out1, double* out2, Slice& s) {
  const size_t I = geno.nrow(), nj = j1 - j0, KK = (size_t) K * K;
  // out2 holds K values per locus in EM mode, a K x K matrix otherwise
  std::fill(out1, out1 + K * nj, 0.0);
  std::fill(out2, out2 + (EM ? K : KK) * nj, 0.0);
  // a and wa are not const: in gradient mode they are overwritten in place
  // with a - b and wa + wb (see below)
  double* __restrict a = s.a.data();
  double* __restrict b = s.b.data();
  double* __restrict wa = s.wa.data();
  const double* __restrict wb = s.wb.data();
  double* __restrict wq = s.tmp.data();
  for (size_t i0 = 0; i0 < I; i0 += TILE_I) {
    const size_t i1 = std::min(I, i0 + TILE_I), ni = i1 - i0;
    for (size_t jj = 0; jj < nj; jj++) {
      const size_t j = j0 + jj;
      // the EM step needs only a and b; gradient mode also needs wa and wb
      slice_coef<!EM>(geno, qT, p + K * j, K, j, i0, i1, s);
      // output column of this locus (K values, or a K x K matrix)
      double* o1 = out1 + K * jj;
      double* o2 = out2 + (EM ? K : KK) * jj;
      if (EM) {
        // two dot products over the slice per population: sum_i a q_ki and
        // sum_i b q_ki (the factors p_kj and 1 - p_kj are applied in em)
        for (int k = 0; k < K; k++) {
          const double* __restrict qk = qT + I * k + i0;
          double sa = 0.0, sb = 0.0;
#pragma omp simd reduction(+:sa, sb)
          for (size_t ii = 0; ii < ni; ii++) {
            sa += a[ii] * qk[ii];
            sb += b[ii] * qk[ii];
          }
          o1[k] += sa;
          o2[k] += sb;
        }
      } else {
        // a <- a - b and wa <- wa + wb
        // On the P side the coefficients enter only through these
        // combinations, so combining them first halves the work below.
#pragma omp simd
        for (size_t ii = 0; ii < ni; ii++) {
          a[ii] -= b[ii];
          wa[ii] += wb[ii];
        }
        for (int k = 0; k < K; k++) {
          const double* __restrict qk = qT + I * k + i0;
          // gradient: sum_i (a - b) q_ki
          double sg = 0.0;
#pragma omp simd reduction(+:sg)
          for (size_t ii = 0; ii < ni; ii++) sg += a[ii] * qk[ii];
          o1[k] += sg;
          // Hessian column k, rows k2 >= k: sum_i [(wa + wb) q_ki] q_k2i.
          // The bracket depends only on k, so it is formed once in wq and
          // each entry is then a single dot product.
          // Hk[k2] is H[k2 + K * k], i.e. the half with row k2 >= column k,
          // the same half as on the Q side (see mirror_hess()).
#pragma omp simd
          for (size_t ii = 0; ii < ni; ii++) wq[ii] = wa[ii] * qk[ii];
          double* Hk = o2 + K * k;
          for (int k2 = k; k2 < K; k2++) {
            const double* __restrict qk2 = qT + I * k2 + i0;
            double sh = 0.0;
#pragma omp simd reduction(+:sh)
            for (size_t ii = 0; ii < ni; ii++) sh += wq[ii] * qk2[ii];
            Hk[k2] += sh;
          }
        }
      }
    }
  }
  // fill in the other half of each symmetric Hessian once all tiles of
  // individuals have been accumulated
  if (!EM)
    for (size_t jj = 0; jj < nj; jj++) mirror_hess(out2 + KK * jj, K);
}

// x^e for e in {0, 1, 2} (Lagrange interpolation in e), without control flow;
// gives 1 for e = 0 even if x = 0
// The three terms are the Lagrange basis polynomials in e for the nodes
// 0, 1, 2, multiplied by x^0, x^1 and x^2. For other values of e the result
// is meaningless, but loglik_slice() never uses it in that case.
static inline double pow012(double x, double e) {
  return 0.5 * (e - 1.0) * (e - 2.0) + x * e * (2.0 - e) +
    0.5 * x * x * e * (e - 1.0);
}

// Products of up to LOG_BATCH factors of at least LOG_MIN_Z cannot underflow
// (LOG_BATCH^2 must equal TILE_I)
// (8 factors >= 1e-30 give a product >= 1e-240, well above the smallest
// normal double, about 2.2e-308. Each factor is also <= 1, so there is no
// overflow either.)
static const int LOG_BATCH = 8;
static const double LOG_MIN_Z = 1e-30;
// loglik_slice() splits a slice of up to TILE_I individuals into LOG_BATCH
// interleaved products of LOG_BATCH factors each, so it silently assumes
// LOG_BATCH * LOG_BATCH == TILE_I. Check this at compile time, so that
// changing one constant without the other fails to build instead of giving
// wrong log-likelihoods. (Implementation detail with no Julia counterpart.)
static_assert(LOG_BATCH * LOG_BATCH == TILE_I,
              "loglik_slice() requires LOG_BATCH^2 == TILE_I");

// Log-likelihood of one slice (after slice_geno). Each observed genotype
// contributes g log(x) + h log(y), with x = (g > 0) ? qp : 1 and
// y = (h > 0) ? 1 - qp : 1 (h = n - g). For hard calls g and h are 0, 1 or 2,
// so the contribution is log(z) with z = qp^g (1 - qp)^h; the z of LOG_BATCH
// genotypes are multiplied together before taking a single log, which makes
// this several times cheaper than one or two logs per genotype. Any other
// genotype (or one with a z so small that the product could underflow) is
// added with the direct formula.
// Missing genotypes have g = h = 0, hence z = 1 and log(z) = 0.
// The result is the sum over the ni individuals of the slice; s.a and s.b
// are used as scratch for z and c.
static double loglik_slice(size_t ni, Slice& s) {
  const double* __restrict qp = s.qp.data();
  const double* __restrict g = s.g.data();
  const double* __restrict h = s.h.data();
  double* __restrict z = s.a.data();
  double* __restrict c = s.b.data();
  // Pass 1 (vectorised): z = qp^g (1 - qp)^h, assuming integer g, h <= 2,
  // and a flag c that is non-zero when that assumption fails (dosages,
  // e.g. from imputed or scaled genotypes)
#pragma omp simd
  for (size_t ii = 0; ii < ni; ii++) {
    z[ii] = pow012(qp[ii], g[ii]) * pow012(1.0 - qp[ii], h[ii]);
    // c == 0 iff g and h are both 0, 1 or 2
    c[ii] = std::fabs(g[ii] * (g[ii] - 1.0) * (g[ii] - 2.0)) +
      std::fabs(h[ii] * (h[ii] - 1.0) * (h[ii] - 2.0));
  }
  // the exceptions are rare, so this branch is almost always predicted
  // Pass 2 (scalar): genotypes that cannot go into the batched product are
  // added with explicit logs, and their z is replaced by the neutral 1. This
  // includes z = 0 (qp of exactly 0 with g > 0, or 1 with h > 0), which
  // correctly gives -Inf. The | (not ||) avoids a second branch.
  double ll = 0.0;
  for (size_t ii = 0; ii < ni; ii++) {
    if ((z[ii] < LOG_MIN_Z) | (c[ii] != 0)) {
      z[ii] = 1.0;
      if (g[ii] > 0) ll += g[ii] * std::log(qp[ii]);
      if (h[ii] > 0) ll += h[ii] * std::log(1.0 - qp[ii]);
    }
  }
  // LOG_BATCH interleaved products of LOG_BATCH factors each (the slice is
  // padded with ones), so that the products vectorise too
  // Product b collects z[b], z[b + 8], z[b + 16], ..., so the inner loop
  // over b multiplies 8 independent accumulators by 8 contiguous values.
  // Padding a short last slice with ones (log 1 = 0) keeps the loop bounds
  // fixed; z points into s.a, which has TILE_I entries.
  std::fill(z + ni, z + TILE_I, 1.0);
  double prod[LOG_BATCH];
  for (int b = 0; b < LOG_BATCH; b++) prod[b] = 1.0;
  for (int m = 0; m < LOG_BATCH; m++)
#pragma omp simd
    for (int b = 0; b < LOG_BATCH; b++) prod[b] *= z[LOG_BATCH * m + b];
  // 8 logs instead of up to 2 * TILE_I
  for (int b = 0; b < LOG_BATCH; b++) ll += std::log(prod[b]);
  return ll;
}

// Log-likelihood L(q, p). Parallelised over blocks of TILE_J loci; within a
// block, loci are visited in order and, for each locus, individuals in tiles
// of TILE_I. Each block writes its own partial sum, so no atomic or
// reduction clause is needed and the total is deterministic.
// (Julia uses an atomic add of per-task results, whose order can vary.)
double loglik(const GenoView& geno, const double* q, const double* p, int K,
              int n_cores) {
  const size_t I = geno.nrow(), J = geno.ncol();
  std::vector<double> qT = transpose(q, K, I);
  // signed loop counter, as required by older OpenMP versions
  const ptrdiff_t nb = n_blocks(J, TILE_J);
  // one partial sum per block, added up in order: deterministic
  std::vector<double> partial(nb, 0.0);
#pragma omp parallel num_threads(n_cores)
{
  // thread-local slice buffers, allocated once per thread
  Slice s;
  // dynamic scheduling balances the load if some threads are slower
  // (e.g. because the backing file is still being paged in)
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    const size_t j1 = std::min(J, (bl + 1) * TILE_J);
    double ll = 0.0;
    for (size_t j = bl * TILE_J; j < j1; j++) {
      for (size_t i0 = 0; i0 < I; i0 += TILE_I) {
        const size_t i1 = std::min(I, i0 + TILE_I);
        // qp and genotypes only: the log-likelihood needs no coefficients
        slice_geno(geno, qT.data(), p + K * j, K, j, i0, i1, s);
        ll += loglik_slice(i1 - i0, s);
      }
    }
    partial[bl] = ll;
  }
}
  // serial sum in block order, independent of the thread count
  double ll = 0.0;
  for (ptrdiff_t bl = 0; bl < nb; bl++) ll += partial[bl];
  return ll;
}

// Binomial deviance of the genotypes visible in the view, for the
// cross-validation error of ADMIXTURE (Alexander & Lange 2011).
//
// For a held-out genotype g of an individual with n allele copies, the model
// fitted without it predicts mu = n qp (the expected number of alt alleles),
// and the prediction error is measured by the binomial deviance
//   d(g, mu) = g log(g / mu) + (n - g) log((n - g) / (n - mu))
//            = g log(g / (n qp)) + h log(h / (n (1 - qp))),   h = n - g,
// with the convention 0 log 0 = 0. This is the formula given in the ADMIXTURE
// paper and manual for diploids (n = 2), where it is (half) the squared
// deviance residual of the binomial model; using n_i instead of 2 extends it
// to pseudohaploids in the natural way (a deviation from ADMIXTURE, which
// only handles diploids). d is 0 when the prediction is perfect and grows as
// mu moves away from g.
//
// Called with a held-out view (CV_HELD_OUT in radmix_geno_access.h),
// only the genotypes of the scored fold are visible; all the others look
// missing and contribute nothing. Writes the total deviance to *dev and the
// number of visible genotypes to *count, so that the cross-validation error can
// be pooled over folds as sum(dev) / sum(count), i.e. averaged over all
// held-out genotypes, as ADMIXTURE does.
//
// Same structure as loglik(): blocks of loci in parallel, one partial
// sum per block added up serially, so the result does not depend on the
// number of threads.
void deviance(const GenoView& geno, const double* q, const double* p, int K,
              int n_cores, double* dev, double* count) {
  const size_t I = geno.nrow(), J = geno.ncol();
  std::vector<double> qT = transpose(q, K, I);
  const ptrdiff_t nb = n_blocks(J, TILE_J);
  std::vector<double> partial_dev(nb, 0.0), partial_n(nb, 0.0);
#pragma omp parallel num_threads(n_cores)
{
  Slice s;
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    const size_t j1 = std::min(J, (bl + 1) * TILE_J);
    double d = 0.0, n_obs = 0.0;
    for (size_t j = bl * TILE_J; j < j1; j++) {
      for (size_t i0 = 0; i0 < I; i0 += TILE_I) {
        const size_t i1 = std::min(I, i0 + TILE_I);
        // qp and genotypes; hidden and missing genotypes have g = h = 0
        slice_geno(geno, qT.data(), p + K * j, K, j, i0, i1, s);
        for (size_t ii = 0; ii < i1 - i0; ii++) {
          const double g = s.g[ii], h = s.h[ii], qp = s.qp[ii];
          // a visible genotype always has g + h = n_i > 0
          if (g + h <= 0) continue;
          // n_i, recovered as g + h (the slice does not store n_trials)
          const double n = g + h;
          // skip zero coefficients (0 log 0 = 0); the estimates are kept in
          // [1e-5, 1 - 1e-5], so qp is never exactly 0 or 1 here
          if (g > 0) d += g * std::log(g / (n * qp));
          if (h > 0) d += h * std::log(h / (n * (1.0 - qp)));
          n_obs += 1.0;
        }
      }
    }
    partial_dev[bl] = d;
    partial_n[bl] = n_obs;
  }
}
  *dev = 0.0;
  *count = 0.0;
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    *dev += partial_dev[bl];
    *count += partial_n[bl];
  }
}

// Two passes from the same (q, p): a Q-side pass over blocks of individuals
// and a P-side pass over blocks of loci. This reads the genotypes twice but
// keeps both passes free of reductions.
//
// Derivation (FRAPPE EM). Each allele copy of individual i at locus j comes
// from a latent population k. Given an alt allele, the posterior probability
// that it came from k is q_ki p_kj / qp; given a reference allele it is
// q_ki (1 - p_kj) / (1 - qp). The E step therefore assigns
//   g q_ki p_kj / qp + h q_ki (1 - p_kj) / (1 - qp)
//     = q_ki [a p_kj + b (1 - p_kj)]
// expected allele copies to population k, of which g q_ki p_kj / qp =
// q_ki p_kj a are alt alleles. The M step sets q_ki to the share of
// individual i's allele copies assigned to k, and p_kj to the share of
// alt alleles among the copies at locus j assigned to k.
void em(const GenoView& geno, const double* q, const double* p, int K,
        double* q_new, double* p_new, int n_cores) {
  const size_t I = geno.nrow(), J = geno.ncol();
  // shared, read-only inputs prepared before the parallel region
  std::vector<double> omp = one_minus(p, K * J);
  std::vector<double> qT = transpose(q, K, I);
  const ptrdiff_t nbi = n_blocks(I, TILE_I), nbj = n_blocks(J, TILE_J);
  // a single parallel region containing two worksharing loops, so the
  // threads (and their Slice buffers) are created only once
  // Q pass with few individuals: compute the sums over loci with the
  // two-dimensional blocking (accum_q_segmented()) before the parallel
  // region below; the normalisation is then done inside it
  const size_t n_seg = q_segments(I, J);
  if (n_seg > 1)
    accum_q_segmented<false>(geno, qT.data(), p, omp.data(), K, n_seg, q_new,
                             NULL, n_cores);
#pragma omp parallel num_threads(n_cores)
{
  Slice s;
  std::vector<double> gT(K * TILE_I);
  // Q: q_new_ki = q_ki * sum_j (a p_kj + b (1 - p_kj)) / (n_i J_i), where
  // J_i is the number of loci observed in individual i. The denominator equals
  // the column sum of the numerator, so we normalise by that: this keeps q on
  // the simplex with missing data (Julia divides by n_i J, which does not).
  // (Summing the numerator over k gives, for each observed locus,
  // a qp + b (1 - qp) = g + h = n_i, hence n_i J_i in total. Julia's em_q!()
  // divides by 2J for every individual, which equals the column sum only
  // when nothing is missing and everyone is diploid.)
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < nbi; bl++) {
    const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
    // the sums over j are written straight into q_new (columns i0..i1-1),
    // which this thread owns; no Hessian is needed (with segments, they have
    // already been computed above)
    if (n_seg == 1)
      accum_q_block<false>(geno, qT.data(), p, omp.data(), K, i0, i1, 0, J,
                           q_new + K * i0, NULL, s, gT.data(), NULL);
    for (size_t i = i0; i < i1; i++) {
      // multiply by q_ki (taken out of the sum over j) and normalise
      double denom = 0.0;
      for (int k = 0; k < K; k++) {
        q_new[k + K * i] *= q[k + K * i];
        denom += q_new[k + K * i];
      }
      // denom is 0 only if the individual has no observed genotype (or q is
      // degenerate): then there is no information and q is left unchanged
      if (denom > 0)
        for (int k = 0; k < K; k++) q_new[k + K * i] /= denom;
      else
        for (int k = 0; k < K; k++) q_new[k + K * i] = q[k + K * i];
    }
  }
  // (the implicit barrier at the end of the omp for above separates the two
  // passes; they write different outputs, so it is not strictly needed)
  // P: p_new = A / (A + B), A = p sum_i a q, B = (1 - p) sum_i b q
  // (A and B are the expected numbers of alt and reference allele copies at
  // locus j assigned to population k. Note that the P pass uses the old q,
  // not q_new, matching Julia's "base" mode.)
  std::vector<double> sa(K * TILE_J), sb(K * TILE_J);
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < nbj; bl++) {
    const size_t j0 = bl * TILE_J, j1 = std::min(J, j0 + TILE_J);
    // sa = sum_i a q_ki and sb = sum_i b q_ki for the loci of the block
    accum_p_block<true>(geno, qT.data(), p, K, j0, j1, sa.data(), sb.data(),
                        s);
    for (size_t j = j0; j < j1; j++) {
      for (int k = 0; k < K; k++) {
        // x indexes the full K x J matrices, y the block-local K x nj buffers
        const size_t x = k + K * j, y = k + K * (j - j0);
        const double A = p[x] * sa[y], B = omp[x] * sb[y];
        // A + B = 0 means no observed allele copies are assigned to k at
        // this locus (e.g. all genotypes missing): keep the old frequency
        // Unlike the Newton updates, the result is not projected, so a
        // locus where every observed genotype is homozygous reference gets
        // p_new = 0 exactly (A = 0), and one fixed for the alt allele
        // p_new = 1. The kernels that follow cope with this through the
        // guards on qp of exactly 0 or 1 (slice_coef(), loglik_slice()).
        p_new[x] = (A + B > 0) ? A / (A + B) : p[x];
      }
    }
  }
}
}

// Gradient (K x I) and negative Hessian (K x K x I) of L with respect to q,
// for every individual. Parallelised over blocks of TILE_I individuals; each
// block writes its own columns of grad and its own K x K slabs of hess.
void grad_q(const GenoView& geno, const double* q, const double* p, int K,
            double* grad, double* hess, int n_cores) {
  const size_t I = geno.nrow(), J = geno.ncol(), KK = (size_t) K * K;
  std::vector<double> omp = one_minus(p, K * J);
  std::vector<double> qT = transpose(q, K, I);
  const ptrdiff_t nb = n_blocks(I, TILE_I);
  // with few individuals, split the loci into segments as well
  const size_t n_seg = q_segments(I, J);
  if (n_seg > 1) {
    accum_q_segmented<true>(geno, qT.data(), p, omp.data(), K, n_seg, grad,
                            hess, n_cores);
    return;
  }
#pragma omp parallel num_threads(n_cores)
{
  // thread-local slice buffers and transposed accumulators
  Slice s;
  std::vector<double> gT(K * TILE_I), hT(K * (K + 1) / 2 * TILE_I);
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
    accum_q_block<true>(geno, qT.data(), p, omp.data(), K, i0, i1, 0, J,
                        grad + K * i0, hess + KK * i0, s, gT.data(),
                        hT.data());
  }
}
}

// Gradient (K x J) and negative Hessian (K x K x J) of L with respect to p,
// for every locus. Parallelised over blocks of TILE_J loci; each block writes
// its own columns of grad and its own K x K slabs of hess. No 1 - p copy is
// needed: on the P side the coefficients are combined as a - b and wa + wb.
void grad_p(const GenoView& geno, const double* q, const double* p, int K,
            double* grad, double* hess, int n_cores) {
  const size_t I = geno.nrow(), J = geno.ncol(), KK = (size_t) K * K;
  std::vector<double> qT = transpose(q, K, I);
  const ptrdiff_t nb = n_blocks(J, TILE_J);
#pragma omp parallel num_threads(n_cores)
{
  Slice s;
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    const size_t j0 = bl * TILE_J, j1 = std::min(J, j0 + TILE_J);
    accum_p_block<false>(geno, qT.data(), p, K, j0, j1, grad + K * j0,
                         hess + KK * j0, s);
  }
}
}

// For each block, the gradient and Hessian are accumulated into thread-local
// buffers and the per-column QPs are solved straight away, so the full
// K x K x I Hessian is never stored.
//
// For individual i the QP (radmix_qp_sweep.h) finds the step delta that
// maximises the local quadratic model grad' delta - 0.5 delta' H delta, where H
// is the negative Hessian, subject to 0 <= q_.i + delta <= 1 and
// sum_k delta_k = 1 - sum_k q_ki (simplex = true). q_.i + delta is then
// projected onto {sum = 1, q >= eps} to guard against round-off and to keep
// every q_ki strictly positive, so that qp never becomes exactly 0 in later
// iterations. As in Julia, if the QP fails, delta is whatever the solver
// reached and the projection still returns a feasible column.
// Deviation from OpenADMIXTURE.jl: discard a non-finite Newton step.
// Julia adds the QP step to the current column whatever its value. If the
// gradient or Hessian of a column contains Inf/NaN (e.g. a genotype that is
// impossible under the current estimates, with qp of exactly 0 or 1), or
// the QP solver hits a vanishing pivot, delta can contain Inf/NaN, which
// would then spread to the whole column and, through later iterations, to
// the log-likelihood. We instead keep the current value of that column for
// this sweep (delta = 0): the other columns are still updated, and the
// column gets another chance at the next sweep, when its neighbours have
// moved. Returns true if the step was discarded.
static inline bool discard_nonfinite_step(double* delta, int K) {
  for (int k = 0; k < K; k++) {
    if (!std::isfinite(delta[k])) {
      std::fill(delta, delta + K, 0.0);
      return true;
    }
  }
  return false;
}

void update_q(const GenoView& geno, const double* q, const double* p, int K,
              double* q_new, int n_cores) {
  const size_t I = geno.nrow(), J = geno.ncol(), KK = (size_t) K * K;
  std::vector<double> omp = one_minus(p, K * J);
  std::vector<double> qT = transpose(q, K, I);
  const ptrdiff_t nb = n_blocks(I, TILE_I);
  // With few individuals, use the two-dimensional blocking: compute the
  // gradient and Hessian of every individual first (they are small, as I is
  // small), then solve the QPs in parallel over individuals
  const size_t n_seg = q_segments(I, J);
  if (n_seg > 1) {
    std::vector<double> grad(K * I), hess(KK * I);
    accum_q_segmented<true>(geno, qT.data(), p, omp.data(), K, n_seg,
                            grad.data(), hess.data(), n_cores);
#pragma omp parallel num_threads(n_cores)
{
    // per-thread scratch: the step and the QP workspace
    std::vector<double> delta(K);
    QPWork work(K);
    // one small QP per individual; dynamic scheduling in chunks of 16
    // individuals balances QPs that need different numbers of iterations
#pragma omp for schedule(dynamic, 16)
    for (ptrdiff_t i = 0; i < (ptrdiff_t) I; i++) {
      const double* qi = q + K * i;
      // the same steps as in the single-pass path below: constrained Newton
      // step, discard it if not finite, add it and project onto the simplex
      qp_solve_column(&hess[KK * i], &grad[K * i], qi, K, true, delta.data(),
                      work);
      discard_nonfinite_step(delta.data(), K);
      double* qn = q_new + K * i;
      for (int k = 0; k < K; k++) qn[k] = qi[k] + delta[k];
      project_q(qn, K);
    }
}
    return;
  }
#pragma omp parallel num_threads(n_cores)
{
  // per-thread scratch: slice buffers, transposed accumulators, the gradient
  // and Hessian of one block of individuals, the step, and the QP workspace
  Slice s;
  std::vector<double> gT(K * TILE_I), hT(K * (K + 1) / 2 * TILE_I);
  std::vector<double> grad(K * TILE_I), hess(KK * TILE_I), delta(K);
  QPWork work(K);
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    const size_t i0 = bl * TILE_I, i1 = std::min(I, i0 + TILE_I);
    // gradient and negative Hessian for individuals [i0, i1)
    accum_q_block<true>(geno, qT.data(), p, omp.data(), K, i0, i1, 0, J,
                        grad.data(), hess.data(), s, gT.data(), hT.data());
    for (size_t i = i0; i < i1; i++) {
      const size_t ii = i - i0;
      const double* qi = q + K * i;
      // constrained Newton step for individual i (the return value, the
      // number of iterations or 0 on failure, is not used, as in Julia)
      qp_solve_column(&hess[KK * ii], &grad[K * ii], qi, K, true,
                      delta.data(), work);
      // keep the current q_.i if the step is not finite (a deviation from
      // Julia, see discard_nonfinite_step())
      discard_nonfinite_step(delta.data(), K);
      double* qn = q_new + K * i;
      for (int k = 0; k < K; k++) qn[k] = qi[k] + delta[k];
      project_q(qn, K);
    }
  }
}
}

// Same as update_q, for the loci: for each block of TILE_J loci the
// gradient and Hessian are accumulated, then for each locus the QP with box
// constraints 0 <= p_.j + delta <= 1 (simplex = false: the frequencies of
// one locus across populations are not constrained to sum to anything) is
// solved and p_.j + delta is clamped to [eps, 1 - eps], so that qp stays
// strictly inside (0, 1) in later iterations.
void update_p(const GenoView& geno, const double* q, const double* p, int K,
              double* p_new, int n_cores) {
  const size_t I = geno.nrow(), J = geno.ncol(), KK = (size_t) K * K;
  std::vector<double> qT = transpose(q, K, I);
  const ptrdiff_t nb = n_blocks(J, TILE_J);
#pragma omp parallel num_threads(n_cores)
{
  Slice s;
  std::vector<double> grad(K * TILE_J), hess(KK * TILE_J), delta(K);
  QPWork work(K);
#pragma omp for schedule(dynamic)
  for (ptrdiff_t bl = 0; bl < nb; bl++) {
    const size_t j0 = bl * TILE_J, j1 = std::min(J, j0 + TILE_J);
    // gradient and negative Hessian for loci [j0, j1)
    accum_p_block<false>(geno, qT.data(), p, K, j0, j1, grad.data(),
                         hess.data(), s);
    for (size_t j = j0; j < j1; j++) {
      const size_t jj = j - j0;
      const double* pj = p + K * j;
      qp_solve_column(&hess[KK * jj], &grad[K * jj], pj, K, false,
                      delta.data(), work);
      // keep the current p_.j if the step is not finite (a deviation from
      // Julia, see discard_nonfinite_step())
      discard_nonfinite_step(delta.data(), K);
      double* pn = p_new + K * j;
      for (int k = 0; k < K; k++) pn[k] = pj[k] + delta[k];
      project_p(pn, K);
    }
  }
}
}

}  // namespace radmix

using namespace radmix;

// ---------------------------------------------------------------------------
// R interface
//
// The wrappers below are the only code in this file that touches R. Each
// one builds the GenoView from the dataset prepared in R (`data`, see
// radmix_prepare() and radmix_geno_access.h), checks every
// dimension with Rcpp::stop() (safe here, on the main thread, before any
// parallel region starts), allocates the R outputs, and only then calls the
// kernel, which writes into the output memory directly. q and p are passed as
// NumericMatrix: Rcpp wraps a double matrix without copying it (an integer
// matrix is coerced to a new double one), and the kernels only read them.
// Inside, each kernel makes its own C++ copies (the transposed q and 1 - p)
// before the parallel region, so the threads only ever read plain C++ memory
// or the read-only q and p.

// check q (K x I) and p (K x J) against the genotypes; returns K
// The kernels trust their inputs completely (raw pointer arithmetic), so
// every size must be checked here, before entering C++ code that cannot
// report errors.
static int check_dims(const GenoView& geno, const NumericMatrix& q,
                      const NumericMatrix& p, int n_cores) {
  const int K = q.nrow();
  if (K < 1) stop("'q' must have at least one row.");
  if (p.nrow() != K) stop("'q' and 'p' must have the same number of rows (K).");
  if ((size_t) q.ncol() != geno.nrow())
    stop("'q' must have one column per individual (%d), not %d.",
         (int) geno.nrow(), q.ncol());
  if ((size_t) p.ncol() != geno.ncol())
    stop("'p' must have one column per locus (%d), not %d.",
         (int) geno.ncol(), p.ncol());
  if (geno.nrow() == 0 || geno.ncol() == 0)
    stop("There must be at least one individual and one locus.");
  if (n_cores < 1) stop("'n_cores' must be at least 1.");
  return K;
}

// Log-likelihood (a scalar).
// [[Rcpp::export]]
double radmix_loglik_cpp(List data, NumericMatrix q, NumericMatrix p,
                         int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  return loglik(geno, q.begin(), p.begin(), K, n_cores);
}

// Number of segments of loci used by the Q-side kernels for a dataset with
// I individuals and J loci (1 = no segmentation; see q_segments()). Only
// used in tests, to check that both code paths are exercised.
// [[Rcpp::export]]
int radmix_q_segments_cpp(double I, double J) {
  return (int) q_segments((size_t) I, (size_t) J);
}

// Binomial deviance of the visible genotypes, for cross-validation; returns
// c(deviance, count) (see deviance()).
// [[Rcpp::export]]
NumericVector radmix_deviance_cpp(List data, NumericMatrix q,
                                  NumericMatrix p, int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  double dev, count;
  deviance(geno, q.begin(), p.begin(), K, n_cores, &dev, &count);
  return NumericVector::create(dev, count);
}

// One EM step; returns list(q = K x I, p = K x J).
// [[Rcpp::export]]
List radmix_em_cpp(List data, NumericMatrix q, NumericMatrix p,
                   int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  // outputs allocated in R memory up front (never inside the parallel region)
  NumericMatrix q_new(K, q.ncol()), p_new(K, p.ncol());
  em(geno, q.begin(), p.begin(), K, q_new.begin(), p_new.begin(),
          n_cores);
  return List::create(_["q"] = q_new, _["p"] = p_new);
}

// Q-block derivatives; returns list(grad = K x I, hess = K x K x I array).
// [[Rcpp::export]]
List radmix_grad_q_cpp(List data, NumericMatrix q, NumericMatrix p,
                       int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  const int I = q.ncol();
  NumericMatrix grad(K, I);
  // a numeric vector with a dim attribute is an R array; the size is
  // computed in size_t to avoid int overflow for large K * K * I
  NumericVector hess((size_t) K * K * I);
  hess.attr("dim") = IntegerVector::create(K, K, I);
  grad_q(geno, q.begin(), p.begin(), K, grad.begin(), hess.begin(),
              n_cores);
  return List::create(_["grad"] = grad, _["hess"] = hess);
}

// P-block derivatives; returns list(grad = K x J, hess = K x K x J array).
// [[Rcpp::export]]
List radmix_grad_p_cpp(List data, NumericMatrix q, NumericMatrix p,
                       int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  const int J = p.ncol();
  NumericMatrix grad(K, J);
  NumericVector hess((size_t) K * K * J);
  hess.attr("dim") = IntegerVector::create(K, K, J);
  grad_p(geno, q.begin(), p.begin(), K, grad.begin(), hess.begin(),
              n_cores);
  return List::create(_["grad"] = grad, _["hess"] = hess);
}

// Newton (QP) update of q; returns the new K x I matrix.
// [[Rcpp::export]]
NumericMatrix radmix_update_q_cpp(List data, NumericMatrix q,
                                  NumericMatrix p, int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  NumericMatrix q_new(K, q.ncol());
  update_q(geno, q.begin(), p.begin(), K, q_new.begin(), n_cores);
  return q_new;
}

// Newton (QP) update of p; returns the new K x J matrix.
// [[Rcpp::export]]
NumericMatrix radmix_update_p_cpp(List data, NumericMatrix q,
                                  NumericMatrix p, int n_cores) {
  const GenoView geno(data);
  const int K = check_dims(geno, q, p, n_cores);
  NumericMatrix p_new(K, p.ncol());
  update_p(geno, q.begin(), p.begin(), K, p_new.begin(), n_cores);
  return p_new;
}

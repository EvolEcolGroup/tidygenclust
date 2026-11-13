#ifndef RADMIX_GENO_ACCESS_H
#define RADMIX_GENO_ACCESS_H

// Read-only view of a subset of a bigstatsr FBM.code256 (one raw byte per
// genotype, individuals in rows, loci in columns, column-major).
//
// How tidypopgen stores genotypes
// -------------------------------
// A gen_tibble keeps its genotypes in a file-backed matrix (FBM) of type
// code256: a memory-mapped file (.bk) holding one unsigned byte per
// genotype. Rows are individuals and columns are loci, stored column-major,
// so the byte for individual r and locus c (both 0-based, in the full FBM)
// sits at offset r + c * nrow. Each byte is the count of the alternate
// allele (0, 1, 2), and other byte values encode missing data (3 for
// diploids) or imputed genotypes (4, 5, 6 after gt_impute_simple()). The R
// object carries a 256-entry table (code256) that maps bytes to values, but
// that table can be swapped at runtime (e.g. by gt_set_imputed()), so we do
// not rely on it.
//
// What this class does
// --------------------
// A gen_tibble usually shows only a subset of the FBM: some individuals may
// have been filtered out or reordered, and loci may have been selected with
// select_loci(). The view stores the selected rows and columns, so that
// get(i, j) returns the genotype of the i-th individual and the j-th locus
// of the gen_tibble, whatever their position in the FBM.
//
// Bytes are decoded with a 256-entry lookup table built by radmix_prepare()
// (NaN = missing), so the result never depends on the code256 currently
// attached to the R object. Each individual has a number of trials n_i (2
// for diploids, 1 for pseudohaploids) and a scale s_i applied to the decoded
// dosage (1 for diploids, 0.5 for pseudohaploids, which tidypopgen stores as
// 0/2), so that the genotype is g_ij = s_i lut[byte], a count out of n_i.
// In the ADMIXTURE model, g_ij ~ Binomial(n_i, qp_ij), so a pseudohaploid
// individual contributes a single allele (0 or 1) per locus.
// (tidypopgen::gt_pseudohaploid() only checks the first loci for
// heterozygous calls, so a pseudohaploid can still hold a stored 1 further
// on. It decodes to g = 0.5 with n_i = 1, which the kernels treat as a
// fractional count; the log-likelihood takes its slow, exact path for it,
// see loglik_slice() in radmix_kernels.cpp.)
//
// The view is immutable after construction and therefore safe to share
// between OpenMP threads. Do not call any R API from inside parallel regions.
//
// Cross-validation masks
// ----------------------
// ADMIXTURE's cross-validation (Alexander & Lange 2011) splits the observed
// genotypes into v folds, refits the model with one fold masked (treated as
// missing), and measures how well the refitted model predicts the masked
// genotypes. A view can be built in one of three modes:
//   CV_NONE     all observed genotypes are visible (normal fitting);
//   CV_TRAIN    the genotypes of fold `cv_fold` are hidden (fitting a fold);
//   CV_HELD_OUT only the genotypes of fold `cv_fold` are visible (scoring the
//               predictions of a fold).
// Because every kernel reads genotypes through get(), a training view makes
// all of them (log-likelihood, EM, gradients, Newton updates) ignore the
// masked genotypes, exactly as if they were missing, with no change to the
// kernels themselves.
//
// Which fold a genotype belongs to is not stored: it is computed on the fly
// from a hash of (seed, i, j). This needs no memory (a stored fold matrix
// would be as large as the selected genotypes), is thread safe, and gives
// the same folds for the same seed. Each genotype is assigned to fold
// hash mod v, so folds are of roughly (not exactly) equal size, as in
// ADMIXTURE's description ("roughly equally sized"); genotypes that are
// already missing are never assigned to a fold, because get() rejects them
// first. ADMIXTURE draws its random partition with its own random number
// generator (its source code is not available), so these folds, and hence
// the CV errors, cannot be identical to ADMIXTURE's, only statistically
// equivalent.
//
// The hash costs about half a plain pass per pass over the genotypes (a pass
// with a fold hidden took 1.5 times a plain one). A training view is read by
// every pass of the fit of its fold, so the fit precomputes its visibility
// once, as a bit mask (build_cv_mask(): one bit per genotype, 1/8 of the
// size of the genotypes), and the passes then read one 64-bit word per slice
// of 64 individuals instead of hashing every genotype. The bits are those of
// the hash, so the results are identical with or without the mask.

#include <bigstatsr/BMAcc.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

// All the internal C++ code of gt_radmix() lives in the namespace radmix,
// so that its generic names (GenoView, project_q, QPWork, ...) cannot clash
// with other compiled code in tidygenclust. Only the functions exported to R
// (all prefixed radmix_) are global, as Rcpp requires.
namespace radmix {

// cross-validation modes (see above), with the integer codes used by
// radmix_set_cv() in R
enum CvMode { CV_NONE = 0, CV_TRAIN = 1, CV_HELD_OUT = 2 };

// The view of the genotypes of a gen_tibble (see the top of the file): built
// once per call from R, then read by every kernel through get() (one
// genotype) or slice_bytes() (the raw bytes of a block of individuals).
class GenoView {
public:
  // data: the list prepared in R by radmix_prepare(), with elements
  //   fbm: the FBM (an R reference class object, i.e. an environment);
  //   rows, cols: 1-based indices of the selected individuals and loci in
  //     the FBM, in the order in which they appear in the gen_tibble (they
  //     do not need to be sorted);
  //   lut: 256 values, the decoded dosage of each byte (NaN for missing);
  //   n_trials, scale: one value per selected individual (see above);
  //   cv (optional): the cross-validation mask (see above), as created by
  //     radmix_set_cv(): integer mode (0 = none, 1 = train,
  //     2 = held out), 1-based fold, n_folds and seed. Without it, the view
  //     shows all the observed genotypes.
  // The view is built at the start of every C++ function called from R (as
  // tidypopgen does with its accessors), which costs a pass over the
  // indices, negligible compared with any kernel. It reads the backing file
  // through a raw pointer, so `data` (which holds the FBM, and hence the
  // mapping) must stay alive while the view is used, which is always the
  // case within a call from R.
  // (Strictly, the mapping belongs to the external pointer kept in the
  // FBM's `extptr` field, and R code that resets that field, such as
  // FBM$add_columns(), lets the garbage collector unmap it. No R code runs
  // inside the kernels; the fitting loop only calls message() and checks
  // for interrupts between passes.)
  // All checks happen here, in the constructor, which is always called on
  // the main thread before any parallel region, so that Rcpp::stop() can
  // safely be used.
  explicit GenoView(const Rcpp::List& data) {
    const char* needed[] = {"fbm", "rows", "cols", "lut", "n_trials",
                            "scale"};
    for (const char* name : needed) {
      if (!data.containsElementNamed(name))
        Rcpp::stop("'data' must be prepared by radmix_prepare() "
                   "(element '%s' is missing).", name);
    }
    Rcpp::Environment BM = data["fbm"];
    Rcpp::IntegerVector row_ind = data["rows"], col_ind = data["cols"];
    Rcpp::NumericVector lut = data["lut"], n_trials = data["n_trials"],
      scale = data["scale"];
    // the "address" field of the FBM holds an external pointer to the
    // read-only memory mapping of the backing file
    Rcpp::XPtr<FBM> xpBM = BM["address"];
    // type 1 = unsigned char, i.e. code256; other FBM types (e.g. doubles)
    // cannot be decoded with a byte lookup table
    if (xpBM->matrix_type() != 1)
      Rcpp::stop("The backing file must be an FBM.code256 (raw bytes).");
    if (lut.size() != 256) Rcpp::stop("'lut' must have length 256.");
    size_t n = row_ind.size();
    if ((size_t) n_trials.size() != n || (size_t) scale.size() != n)
      Rcpp::stop("'n_trials' and 'scale' must have one value per individual.");
    // raw pointer to the first byte of the mapped file (see above)
    _pMat = static_cast<const unsigned char*>(xpBM->matrix());
    _nrow_fbm = xpBM->nrow();
    // 1-based R indices -> 0-based offsets (bounds checked)
    _rows = vec_int_to_size(row_ind, xpBM->nrow(), 1);
    std::vector<size_t> cols = vec_int_to_size(col_ind, xpBM->ncol(), 1);
    // precompute the offset of the start of each selected column, so that
    // get() needs a single addition instead of a multiplication
    _col_off.resize(cols.size());
    for (size_t j = 0; j < cols.size(); j++) _col_off[j] = cols[j] * _nrow_fbm;
    // split the lookup table into a missing flag and a value (0 for missing
    // bytes), so that get() never has to test for NaN
    for (int b = 0; b < 256; b++) {
      _missing[b] = ISNAN(lut[b]);
      _lut[b] = _missing[b] ? 0.0 : lut[b];
    }
    // copy the per-individual values out of R memory, so that they can be
    // read from any thread
    _n.assign(n_trials.begin(), n_trials.end());
    _scale.assign(scale.begin(), scale.end());
    // cross-validation mask
    int cv_mode = CV_NONE, cv_fold = 0, n_folds = 1;
    double cv_seed = 0;
    if (data.containsElementNamed("cv")) {
      Rcpp::List cv = data["cv"];
      cv_mode = Rcpp::as<int>(cv["mode"]);
      cv_fold = Rcpp::as<int>(cv["fold"]) - 1;  // 1-based in R
      n_folds = Rcpp::as<int>(cv["n_folds"]);
      cv_seed = Rcpp::as<double>(cv["seed"]);
    }
    if (cv_mode != CV_NONE && cv_mode != CV_TRAIN && cv_mode != CV_HELD_OUT)
      Rcpp::stop("'cv_mode' must be 0 (none), 1 (train) or 2 (held out).");
    if (cv_mode != CV_NONE &&
        (n_folds < 2 || cv_fold < 0 || cv_fold >= n_folds))
      Rcpp::stop("'cv_fold' must be between 1 and 'n_folds', with "
                 "'n_folds' >= 2.");
    // the seed is passed from R as a double holding an integer (R has no
    // 64-bit integers), which is enough to give each run its own folds
    // (and below 2^53, the largest integer that a double holds exactly: a
    // larger or infinite value would overflow the cast to uint64_t below;
    // the comparison is also false for NaN)
    if (!(cv_seed >= 0 && cv_seed < 9007199254740992.0) ||
        cv_seed != std::floor(cv_seed))
      Rcpp::stop("'cv_seed' must be a non-negative integer.");
    // whether a training view may precompute its bit mask (see
    // build_cv_mask(); an optional element of cv, TRUE by default, so that
    // the tests and benchmarks can compare the two ways)
    _cv_mask_allowed = true;
    if (data.containsElementNamed("cv")) {
      Rcpp::List cv = data["cv"];
      if (cv.containsElementNamed("mask"))
        _cv_mask_allowed = Rcpp::as<bool>(cv["mask"]);
    }
    _cv_bits = false;
    _mask_words = 0;
    _cv_mode = cv_mode;
    _cv_fold = (uint64_t) cv_fold;
    _cv_hide_fold = cv_mode == CV_TRAIN;
    // n_folds is only checked for masked views; flooring it at 1 keeps
    // fold() (hash % _n_folds) free of a division by zero in any mode
    _n_folds = (uint64_t) (n_folds > 0 ? n_folds : 1);
    _cv_seed = (uint64_t) cv_seed;
    // Fast path for decoding a block of genotypes (see slice_bytes()). It
    // needs the same decoding for every individual, i.e. a uniform number of
    // trials and scale (always true without pseudohaploids). The tables then
    // give g and h = n - g for each byte directly, with g = h = 0 for
    // missing bytes, exactly the values that get() leads to (see
    // slice_geno() in radmix_kernels.cpp); with a cross-validation mask, the
    // caller also multiplies them by cv_visible().
    bool uniform = true;
    for (size_t i = 1; i < n; i++)
      uniform = uniform && _n[i] == _n[0] && _scale[i] == _scale[0];
    _fast = uniform && n > 0;
    if (_fast) {
      for (int b = 0; b < 256; b++) {
        const double obs = _missing[b] ? 0.0 : 1.0;
        const double g = _lut[b] * _scale[0];
        _gtab[b] = obs * g;
        _htab[b] = obs * (_n[0] - g);
      }
    }
    // _contig[i]: the number of individuals from i onwards whose rows follow
    // each other in the FBM (rows i, i + 1, ... are rows r, r + 1, ...),
    // so that the bytes of a block of them are adjacent in each column
    _contig.assign(n, 1);
    for (size_t i = n; i-- > 1;) {
      if (_rows[i] == _rows[i - 1] + 1) _contig[i - 1] = _contig[i] + 1;
    }
  }

  // number of selected individuals (I) and loci (J), i.e. the dimensions of
  // the genotype matrix as the kernels see it
  size_t nrow() const { return _rows.size(); }
  size_t ncol() const { return _col_off.size(); }

  // number of trials (allele copies) n_i for individual i: 2 for a diploid,
  // 1 for a pseudohaploid
  inline double n_trials(size_t i) const { return _n[i]; }

  // The fast path for a block of genotypes: if individuals i0 to
  // i0 + ni - 1 are adjacent in the FBM and the view decodes every
  // individual alike (see the constructor), returns a pointer to the ni
  // bytes of that block at locus j, and sets gtab and htab to the tables
  // that decode them (g and h = n - g for each byte, 0 for missing ones;
  // with a cross-validation mask, the genotypes for which cv_visible() is
  // false must also be set to 0); otherwise returns nullptr, and the
  // genotypes must be read with get().
  inline const unsigned char* slice_bytes(size_t i0, size_t ni, size_t j,
                                          const double*& gtab,
                                          const double*& htab) const {
    if (!_fast || _contig[i0] < ni) return nullptr;
    gtab = _gtab;
    htab = _htab;
    return _pMat + _rows[i0] + _col_off[j];
  }

  // Decode genotype (i, j). Returns whether it is observed (and visible,
  // with a cross-validation mask), and sets g to the (scaled) alt allele
  // count in [0, n_i] (0 for a missing genotype; for a hidden one, its
  // count, which the caller must ignore).
  // This is the innermost operation of every kernel, so it is kept to a
  // memory read, a table lookup and a multiplication, without branches.
  inline bool get(size_t i, size_t j, double& g) const {
    unsigned char b = _pMat[_rows[i] + _col_off[j]];
    bool visible = !_missing[b];
    // cross-validation: hide the genotypes of the masked fold (training
    // view) or all the others (held-out view), i.e. those whose membership
    // of the fold matches _cv_hide_fold. The genotypes of a fold are spread
    // at random, so this is computed without a branch (which would be
    // mispredicted for about one genotype in n_folds); the fold of a missing
    // genotype is also computed, and ignored. For plain views the test on
    // the mode is a single, perfectly predictable branch.
    if (_cv_mode != CV_NONE) {
      visible &= cv_visible(i, j);
    }
    g = _lut[b] * _scale[i];
    return visible;
  }

  // whether the view has a cross-validation mask, and whether it shows
  // genotype (i, j) (if it is observed): the test of get() without the
  // decoding, for the callers of slice_bytes(), which must apply the mask
  // themselves (cv_visible() is only meaningful when cv_masked() is true).
  // With a precomputed bit mask (cv_bits(), see build_cv_mask()), the bit
  // is read instead of computing the hash (the branch is the same for every
  // genotype, so it is perfectly predicted).
  inline bool cv_masked() const { return _cv_mode != CV_NONE; }
  inline bool cv_visible(size_t i, size_t j) const {
    if (_cv_bits)
      return (_cv_mask[j * _mask_words + (i >> 6)] >> (i & 63)) & 1;
    return cv_visible_hash(i, j);
  }
  inline bool cv_visible_hash(size_t i, size_t j) const {
    return (fold(i, j) == _cv_fold) != _cv_hide_fold;
  }

  // Whether the bit mask of a training view has been built, and its word
  // for the 64 individuals [i0, i0 + 64) at locus j (i0 a multiple of 64;
  // bit ii is the visibility of individual i0 + ii, and the bits beyond the
  // last individual are 0), so that the callers of slice_bytes() read one
  // word per slice. Only meaningful when cv_bits() is true.
  inline bool cv_bits() const { return _cv_bits; }
  inline uint64_t cv_word(size_t i0, size_t j) const {
    return _cv_mask[j * _mask_words + (i0 >> 6)];
  }

  // The largest bit mask that build_cv_mask() makes (1 GiB, for 8.6e9
  // genotypes, e.g. 10,000 individuals x 860,000 loci); beyond it, the
  // passes compute the hash, which needs no memory.
  static constexpr double CV_MASK_MAX_BYTES = 1073741824.0;

  // Precompute the visibility of every genotype of a training view as a bit
  // mask (see the top of the file), with n_cores threads; nothing happens
  // for other views, when the mask is not allowed (cv element mask = FALSE)
  // or when it would exceed CV_MASK_MAX_BYTES. Called by the fitting loop
  // on the main thread, before any pass, as the view must not change while
  // it is shared between threads. Each locus is filled by one thread, and the
  // bits are those of cv_visible_hash(), so the results of the passes do not
  // depend on the mask or on the number of threads.
  void build_cv_mask(int n_cores) {
    if (_cv_mode != CV_TRAIN || !_cv_mask_allowed || _cv_bits) return;
    const size_t I = nrow(), J = ncol(), words = (I + 63) / 64;
    if ((double) words * J * 8.0 > CV_MASK_MAX_BYTES) return;
    _mask_words = words;
    _cv_mask.assign(words * J, 0);
    uint64_t* mask = _cv_mask.data();
#pragma omp parallel for num_threads(n_cores) schedule(static)
    for (ptrdiff_t j = 0; j < (ptrdiff_t) J; j++) {
      uint64_t* mj = mask + (size_t) j * words;
      for (size_t i = 0; i < I; i++)
        mj[i >> 6] |= (uint64_t) cv_visible_hash(i, j) << (i & 63);
    }
    _cv_bits = true;
  }

  // Fold (0 to n_folds - 1) of genotype (i, j), where i and j index the
  // selected individuals and loci (so folds follow the gen_tibble, not the
  // layout of the backing file). The splitmix64 finaliser turns the
  // combination of seed, i and j into a well mixed 64-bit value: consecutive
  // genotypes get unrelated folds, so folds are spread evenly across
  // individuals and loci.
  inline uint64_t fold(size_t i, size_t j) const {
    // The input to the finaliser combines the seed with i and j, each
    // multiplied by a different large odd constant (the splitmix64
    // increment and a 64-bit prime of xxHash), so that (i, j) and (j, i)
    // give different values; the + 1 keeps row 0 and column 0 from
    // dropping out of the sum. All arithmetic wraps modulo 2^64.
    uint64_t z = _cv_seed + 0x9E3779B97F4A7C15ULL * ((uint64_t) i + 1) +
      0xC2B2AE3D27D4EB4FULL * ((uint64_t) j + 1);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z = z ^ (z >> 31);
    return mod_n_folds(z);
  }

  // z % _n_folds. A 64-bit division by a value only known at run time takes
  // tens of cycles, and this runs for every genotype in every pass of a
  // cross-validation fit; by a compile-time constant, the compiler computes
  // the same remainder with a multiplication and shifts. The usual numbers
  // of folds get a constant (the switch is the same for every genotype, so
  // it is perfectly predicted); any other number is divided at run time.
  inline uint64_t mod_n_folds(uint64_t z) const {
    switch (_n_folds) {
    case 2: return z % 2;
    case 3: return z % 3;
    case 4: return z % 4;
    case 5: return z % 5;
    case 6: return z % 6;
    case 7: return z % 7;
    case 8: return z % 8;
    case 9: return z % 9;
    case 10: return z % 10;
    default: return z % _n_folds;
    }
  }

private:
  const unsigned char* _pMat;    // start of the memory-mapped FBM
  size_t _nrow_fbm;              // number of rows of the full FBM
  std::vector<size_t> _rows;     // 0-based FBM row of each individual
  std::vector<size_t> _col_off;  // byte offset of the start of each locus
  double _lut[256];              // decoded dosage of each byte (0 if missing)
  bool _missing[256];            // whether each byte means "missing"
  std::vector<double> _n;        // number of trials of each individual
  std::vector<double> _scale;    // dosage scale of each individual
  int _cv_mode;                  // CV_NONE, CV_TRAIN or CV_HELD_OUT
  uint64_t _cv_fold;             // the masked (or scored) fold
  bool _cv_hide_fold;            // whether the fold is hidden (training view)
                                 // or the other genotypes (held-out view)
  uint64_t _n_folds;             // number of folds
  uint64_t _cv_seed;             // seed of the fold assignment
  bool _fast;                    // whether slice_bytes() can be used
  double _gtab[256];             // g of each byte for slice_bytes()
  double _htab[256];             // h = n - g of each byte for slice_bytes()
  std::vector<size_t> _contig;   // adjacent rows from each individual on
  bool _cv_mask_allowed;         // whether build_cv_mask() may build a mask
  bool _cv_bits;                 // whether the bit mask has been built
  size_t _mask_words;            // 64-bit words of the mask per locus
  std::vector<uint64_t> _cv_mask;  // the visibility of each genotype of a
                                   // training view, one bit each
};

}  // namespace radmix

#endif

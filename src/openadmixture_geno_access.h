#ifndef OPENADMIXTURE_GENO_ACCESS_H
#define OPENADMIXTURE_GENO_ACCESS_H

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
// Bytes are decoded with our own 256-entry lookup table (NaN = missing), so the
// result never depends on the code256 currently attached to the R object.
// Each individual has a number of trials n_i (2 for diploids, 1 for
// pseudohaploids) and a scale s_i applied to the decoded dosage (1 for
// diploids, 0.5 for pseudohaploids, which tidypopgen stores as 0/2).
// In the ADMIXTURE model, g_ij ~ Binomial(n_i, qp_ij), so a pseudohaploid
// individual contributes a single allele (0 or 1) per locus.
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
// first. Deviation from ADMIXTURE: ADMIXTURE draws a random partition with
// its own random number generator (its source code is not available), so
// our folds, and hence CV errors, cannot be identical to ADMIXTURE's, only
// statistically equivalent.

#include <bigstatsr/BMAcc.h>
#include <cmath>
#include <cstdint>
#include <vector>

// All the internal C++ code of the OpenADMIXTURE implementation lives in the
// namespace openadmixture, so that its generic names (GenoView, project_q,
// QPWork, ...) cannot clash with other compiled code in tidygenclust. Only
// the functions exported to R (all prefixed openadmixture_) are global, as
// Rcpp requires.
namespace openadmixture {


// cross-validation modes (see above)
enum CvMode { CV_NONE = 0, CV_TRAIN = 1, CV_HELD_OUT = 2 };

class GenoView {
public:
  // data: the list prepared in R by openadmixture_prepare(), with elements
  //   fbm: the FBM (an R reference class object, i.e. an environment);
  //   rows, cols: 1-based indices of the selected individuals and loci in
  //     the FBM, in the order in which they appear in the gen_tibble (they
  //     do not need to be sorted);
  //   lut: 256 values, the decoded dosage of each byte (NaN for missing);
  //   n_trials, scale: one value per selected individual (see above);
  //   cv (optional): the cross-validation mask (see above), as created by
  //     openadmixture_set_cv(): integer mode (0 = none, 1 = train,
  //     2 = held out), 1-based fold, n_folds and seed. Without it, the view
  //     shows all the observed genotypes.
  // The view is built at the start of every C++ function called from R (as
  // tidypopgen does with its accessors), which costs a pass over the
  // indices, negligible compared with any kernel. It reads the backing file
  // through a raw pointer, so `data` (which holds the FBM, and hence the
  // mapping) must stay alive while the view is used, which is always the
  // case within a call from R.
  // All checks happen here, in the constructor, which is always called on
  // the main thread before any parallel region, so that Rcpp::stop() can
  // safely be used.
  explicit GenoView(const Rcpp::List& data) {
    const char* needed[] = {"fbm", "rows", "cols", "lut", "n_trials",
                            "scale"};
    for (const char* name : needed) {
      if (!data.containsElementNamed(name))
        Rcpp::stop("'data' must be prepared by openadmixture_prepare() "
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
    if (cv_seed < 0 || cv_seed != std::floor(cv_seed))
      Rcpp::stop("'cv_seed' must be a non-negative integer.");
    _cv_mode = cv_mode;
    _cv_fold = (uint64_t) cv_fold;
    _n_folds = (uint64_t) (n_folds > 0 ? n_folds : 1);
    _cv_seed = (uint64_t) cv_seed;
  }

  // number of selected individuals (I) and loci (J)
  size_t nrow() const { return _rows.size(); }
  size_t ncol() const { return _col_off.size(); }

  // number of trials (allele copies) for individual i
  inline double n_trials(size_t i) const { return _n[i]; }

  // Decode genotype (i, j). Returns false if missing; otherwise sets g to the
  // (scaled) alt allele count in [0, n_i].
  // This is the innermost operation of every kernel, so it is kept to a
  // memory read, a table lookup and a multiplication.
  inline bool get(size_t i, size_t j, double& g) const {
    unsigned char b = _pMat[_rows[i] + _col_off[j]];
    if (_missing[b]) return false;
    // cross-validation: hide the genotypes of the masked fold (training
    // view) or all the others (held-out view). For plain views this is a
    // single, perfectly predictable branch.
    if (_cv_mode != CV_NONE) {
      const bool in_fold = fold(i, j) == _cv_fold;
      if (in_fold == (_cv_mode == CV_TRAIN)) return false;
    }
    g = _lut[b] * _scale[i];
    return true;
  }

  // Fold (0 to n_folds - 1) of genotype (i, j), where i and j index the
  // selected individuals and loci (so folds follow the gen_tibble, not the
  // layout of the backing file). The splitmix64 finaliser turns the
  // combination of seed, i and j into a well mixed 64-bit value: consecutive
  // genotypes get unrelated folds, so folds are spread evenly across
  // individuals and loci.
  inline uint64_t fold(size_t i, size_t j) const {
    uint64_t z = _cv_seed + 0x9E3779B97F4A7C15ULL * ((uint64_t) i + 1) +
      0xC2B2AE3D27D4EB4FULL * ((uint64_t) j + 1);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z = z ^ (z >> 31);
    return z % _n_folds;
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
  uint64_t _n_folds;             // number of folds
  uint64_t _cv_seed;             // seed of the fold assignment
};

}  // namespace openadmixture

#endif

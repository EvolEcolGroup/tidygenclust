# Interface between gt_adamixture() and ADAMIXTURE (Saurina-i-Ricos, Mas
# Montserrat & Ioannidis 2026), run in the python session of reticulate.
#
# adamixture_fit() reads the genotypes into memory, either from the backing
# file of a gen_tibble (only its individuals and loci, in its order) or from
# PLINK files, and fits ADAMIXTURE to them with adamixture.run() (its
# programmatic entry point, using the CPU solver, BR-QN). It limits the number
# of threads of ADAMIXTURE, restores the environment variables that
# adamixture.run() sets, and returns P as the frequencies of the counted
# allele (ADAMIXTURE flips the coding when that allele is the most common one
# overall, and then returns the frequencies of the other allele).

import logging
import os
import warnings

import numpy as np

# the environment variables set by adamixture.run(): python shares its
# environment with R, so they are restored after each fit, so that they do not
# change the number of threads of the rest of the R session
_THREAD_VARS = [
    "MKL_NUM_THREADS", "MKL_MAX_THREADS", "OMP_NUM_THREADS",
    "OMP_MAX_THREADS", "NUMEXPR_NUM_THREADS", "NUMEXPR_MAX_THREADS",
    "OPENBLAS_NUM_THREADS", "OPENBLAS_MAX_THREADS", "KMP_DUPLICATE_LIB_OK",
]


def _read_fbm(bfile, fbm_dim, rows, cols):
    """Genotypes of a gen_tibble, as a loci x individuals uint8 array.

    The backing file of a gen_tibble has one byte per genotype (0, 1 and 2
    alternate alleles, 3 for missing, larger values for imputed genotypes),
    without a header. It holds the whole FBM, column-major (individuals x
    loci, i.e. a row-major loci x individuals array), of which the gen_tibble
    uses some rows (individuals) and columns (loci), in its own order (R
    passes length-one vectors as scalars, hence atleast_1d).
    """
    rows = np.atleast_1d(np.asarray(rows, dtype=np.intp))
    cols = np.atleast_1d(np.asarray(cols, dtype=np.intp))
    fbm_nrow, fbm_ncol = (int(d) for d in fbm_dim)
    assert os.path.getsize(bfile) == fbm_nrow * fbm_ncol, \
        "Backing file doesn't match the dimensions of the FBM!"
    G_fbm = np.memmap(bfile, dtype=np.uint8, mode="r",
                      shape=(fbm_ncol, fbm_nrow))
    G = np.ascontiguousarray(G_fbm[np.ix_(cols, rows)])
    del G_fbm
    # imputed genotypes are treated as missing
    G[G >= 3] = 3
    return G


def _read_bed(prefix):
    """Genotypes of PLINK files, as a loci x individuals uint8 array.

    The genotypes are the counts of the first allele of the bim file (A1), as
    in the reader of ADAMIXTURE, with 3 for missing.
    """
    with open(prefix + ".fam") as f:
        n = sum(1 for _ in f)
    with open(prefix + ".bim") as f:
        m = sum(1 for _ in f)
    bytes_per_locus = (n + 3) // 4
    raw = np.fromfile(prefix + ".bed", dtype=np.uint8)
    assert raw[:3].tolist() == [0x6C, 0x1B, 0x01], \
        "Not a PLINK bed file in SNP-major mode!"
    assert raw.size == 3 + m * bytes_per_locus, \
        "The bed file doesn't match the bim and fam files!"
    raw = raw[3:].reshape(m, bytes_per_locus)
    # the four 2-bit codes of each byte, lowest bits first: 00 homozygous for
    # A1, 01 missing, 10 heterozygous, 11 homozygous for A2
    shifts = np.array([0, 2, 4, 6], dtype=np.uint8)
    codes = ((raw[:, :, None] >> shifts) & 3).reshape(m, 4 * bytes_per_locus)
    return np.array([2, 3, 1, 0], dtype=np.uint8)[codes[:, :n]]


def _counted_mean(G):
    """Mean frequency of the counted allele, as used by ADAMIXTURE to decide
    whether to flip the coding."""
    try:
        from adamixture.src.utils_c import get_mean_unpacked
        return get_mean_unpacked(G)
    except ImportError:
        return G[G != 3].mean() / 2.0


def adamixture_fit(bfile, plink, fbm_dim, rows, cols, K, seed, threads, tol,
                   max_iter, cv, verbose):
    """Fit ADAMIXTURE once, returning a dict with Q, P and the CV error.

    'bfile' is the prefix of PLINK files (with 'plink' True), or the backing
    file of a gen_tibble, with 'fbm_dim' (the numbers of rows and columns of
    its FBM), 'rows' and 'cols' (the 0-based rows and columns of the FBM used
    by the gen_tibble). 'cv' is the number of folds of the cross-validation
    (None for none).
    """
    import adamixture

    if plink:
        G = _read_bed(bfile)
    else:
        G = _read_fbm(bfile, fbm_dim, rows, cols)
    # adamixture.run() flips the coding (of G, in place) if the counted allele
    # is the most common one overall
    flipped = G.shape[0] > 0 and _counted_mean(G) >= 0.5

    def fit():
        return adamixture.run(
            data=G, K=int(K), threads=int(threads), seed=int(seed),
            tol=float(tol), max_iter=int(max_iter),
            cv=None if cv is None else int(cv), no_freqs=False,
            chrom_mode="all", verbose=bool(verbose))

    saved = {var: os.environ.get(var) for var in _THREAD_VARS}
    # the modules of ADAMIXTURE log their progress with loggers named after
    # them (adamixture.*), whatever 'verbose' (which only applies to reading
    # the data), so their parent logger is silenced without verbose
    logger = logging.getLogger("adamixture")
    saved_level = logger.level
    if not verbose:
        logger.setLevel(logging.WARNING)
    try:
        # adamixture.run() only sets environment variables such as
        # OMP_NUM_THREADS, which are read when a library is initialised (and
        # its OpenMP library is loaded when it is imported); threadpoolctl
        # limits the libraries that are loaded at run time, and restores
        # their previous limits on exit
        try:
            from threadpoolctl import threadpool_limits
        except ImportError:
            warnings.warn("threadpoolctl is not installed, so the number of "
                          "threads cannot be limited; reinstall ADAMIXTURE "
                          "with tgc_adamixture_install()")
            res = fit()
        else:
            with threadpool_limits(limits=int(threads)):
                res = fit()
    finally:
        logger.setLevel(saved_level)
        for var, value in saved.items():
            if value is None:
                os.environ.pop(var, None)
            else:
                os.environ[var] = value

    if cv is None:
        Q, P = res
        cv_error = None
    else:
        Q, P, cv_res = res
        cv_error = cv_res["avg"]
    if flipped:
        P = 1.0 - P
    return {"Q": Q, "P": P, "cv": cv_error}

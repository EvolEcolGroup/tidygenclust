# this is derived from https://github.com/Rosemeis/fastmixture/blob/main/fastmixture/main.py
# @author of the original file Jonas Meisner
# modified by Andrea Manica for use with tidygenclust
#
# fastmixture_run() follows the main function of fastmixture (version 1.3.0,
# fastmixture/fastmixture.py), but it returns the estimates to R instead of
# writing them to files, and it can read the genotypes directly from the
# backing file of a gen_tibble. It is called by gt_fastmixture(), which
# passes the arguments as an argparse Namespace with the same fields as the
# command line options of fastmixture (see the parser below), plus 'plink'
# (whether 'bfile' is the prefix of PLINK files or the backing file of a
# gen_tibble), 'fbm_dim' (the numbers of rows and columns of the FBM of the
# gen_tibble), 'rows' and 'cols' (the 0-based rows and columns of the FBM
# used by the gen_tibble, i.e. its individuals and loci, in its order).

# Libraries
import argparse
import os
import sys
import warnings
from datetime import datetime
from time import time

from fastmixture import __version__ as VERSION

### Argparse
# the parser of the command line version, kept for reference: the arguments
# are now created in R
# parser = argparse.ArgumentParser(prog="fastmixture")
# parser.add_argument("--version", action="version",
# 	version=f"v{VERSION}")
# parser.add_argument("-b", "--bfile", metavar="PLINK",
# 	help="Prefix for PLINK files (.bed, .bim, .fam)")
# parser.add_argument("-k", "--K", metavar="INT", type=int,
# 	help="Number of ancestral components")
# parser.add_argument("-t", "--threads", metavar="INT", type=int, default=1,
# 	help="Number of threads (1)")
# parser.add_argument("-s", "--seed", metavar="INT", type=int, default=42,
# 	help="Set random seed (42)")
# parser.add_argument("-o", "--out", metavar="OUTPUT", default="fastmixture",
# 	help="Prefix output name (fastmixture)")
# parser.add_argument("--iter", metavar="INT", type=int, default=1000,
# 	help="Maximum number of iterations in EM (1000)")
# parser.add_argument("--tole", metavar="FLOAT", type=float, default=1e-9,
# 	help="Tolerance in scaled log-likelihood units (1e-9)")
# parser.add_argument("--batches", metavar="INT", type=int, default=32,
# 	help="Number of initial mini-batches (32)")
# parser.add_argument("--supervised", metavar="FILE",
# 	help="Path to population assignment file")
# parser.add_argument("--projection", metavar="FILE",
# 	help="Path to ancestral allele frequencies file")
# parser.add_argument("--cv", metavar="INT", type=int,
# 	help="Number of folds in cross-validation to evaluate model fit")
# parser.add_argument("--cv-tole", metavar="FLOAT", type=float, default=1e-7,
# 	help="Tolerance for CV in scaled log-likelihood units (1e-7)")
# parser.add_argument("--check", metavar="INT", type=int, default=5,
# 	help="Number of iterations between convergence checks (5)")
# parser.add_argument("--subsample", metavar="FLOAT", type=float, default=0.7,
# 	help="Fraction of SNPs to subsample in SVD/ALS (0.7)")
# parser.add_argument("--min-subsample", metavar="INT", type=int, default=50000,
# 	help="Minimum number of SNPs to subsample in SVD/ALS (50000)")
# parser.add_argument("--max-subsample", metavar="INT", type=int, default=500000,
# 	help="Maximum number of SNPs to subsample in SVD/ALS (500000)")
# parser.add_argument("--power", metavar="INT", type=int, default=11,
# 	help="Number of power iterations in randomized SVD (11)")
# parser.add_argument("--chunk", metavar="INT", type=int, default=4096,
# 	help="Number of SNPs in chunk operations (4096)")
# parser.add_argument("--als-iter", metavar="INT", type=int, default=1000,
# 	help="Maximum number of iterations in ALS (1000)")
# parser.add_argument("--als-tole", metavar="FLOAT", type=float, default=1e-4,
# 	help="Tolerance for RMSE of P between iterations (1e-4)")
# parser.add_argument("--no-freqs", action="store_true",
# 	help="Do not save P-matrix")
# parser.add_argument("--random-init", action="store_true",
# 	help="Random initialization of parameters")



##### fastmixture #####
def fastmixture_run(args):
	"""Run fastmixture with at most args.threads threads.

	fastmixture parallelises with OpenMP (and numpy with its BLAS library).
	Their numbers of threads cannot be set with environment variables such as
	OMP_NUM_THREADS, which are only read when a library is initialised (R
	has usually started OpenMP already), and which would stay set in the R
	session, as python shares its environment. Instead, threadpoolctl limits
	the libraries that are loaded at run time, and restores their previous
	limits on exit.
	"""
	try:
		from threadpoolctl import threadpool_limits
	except ImportError:
		warnings.warn("threadpoolctl is not installed, so the number of "
			"threads cannot be limited; reinstall fastmixture with "
			"tgc_fastmixture_install()")
		return _fastmixture_run(args)
	with threadpool_limits(limits=args.threads):
		return _fastmixture_run(args)


#def main():
def _fastmixture_run(args):
	# args = parser.parse_args()
	# if len(sys.argv) < 2:
	# 	parser.print_help()
	# 	sys.exit()
	print("-------------------------------------------------")
	print(f"fastmixture v{VERSION}")
	print("C.G. Santander, A. Refoyo-Martinez and J. Meisner")
	print(f"K={args.K}, seed={args.seed}, batches={args.batches}, threads={args.threads}")
	print("-------------------------------------------------\n")

	# Check input
#	assert args.bfile is not None, "No input data (--bfile)!"
	assert args.K > 1, "Please select K > 1!"
	assert args.threads > 0, "Please select a valid number of threads!"
	assert args.seed >= 0, "Please select a valid seed!"
	assert args.batches > 1, "Please select a valid number of batches!"
	assert args.iter > 0, "Please select a valid number of iterations!"
	assert args.tole >= 0.0, "Please select a valid tolerance!"
	assert args.check > 0, "Please select a valid value for convergence check!"
	assert args.min_subsample > 0, "Please select a valid number of SNPs!"
	assert args.max_subsample > 0, "Please select a valid number of SNPs!"
	assert args.power > 0, "Please select a valid number of power iterations!"
	assert args.chunk > 0, "Please select a valid value for chunk size in SVD!"
	assert args.als_iter > 0, "Please select a valid number of iterations in ALS!"
	assert args.als_tole >= 0.0, "Please select a valid tolerance in ALS!"
	assert (args.subsample > 0.0) and (args.subsample <= 1.0), "Please select a valid fraction!"
	if args.cv is not None:
		assert args.cv > 1, "Please select a valid number of folds for cross-validation!"
		assert (args.supervised is None) and (args.projection is None), \
			"Only unsupervised mode works with cross-validation!"
	start = time()

	# Create log-file of used arguments
	# full = vars(parser.parse_args())
	# deaf = vars(parser.parse_args([]))
	# mand = ["seed", "batches"]
	# with open(f"{args.out}.K{args.K}.s{args.seed}.log", "w") as log:
	# 	log.write(f"fastmixture v{VERSION}\n")
	# 	log.write(f"Time: {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}\n")
	# 	log.write(f"Directory: {os.getcwd()}\n")
	# 	log.write("Options:\n")
	# 	for key in full:
	# 		if full[key] != deaf[key]:
	# 			if type(full[key]) is bool:
	# 				log.write(f"\t--{key}\n")
	# 			else:
	# 				log.write(f"\t--{key} {full[key]}\n")
	# 		elif key in mand:
	# 			log.write(f"\t--{key} {full[key]}\n")
	# del full, deaf, mand

	# (the command line version sets the number of threads of the numerical
	# libraries here, through environment variables; this is done by
	# fastmixture_run() instead, see above)

	# Load numerical libraries
	import numpy as np
	from fastmixture import utils
	from fastmixture import shared

	rng = np.random.default_rng(args.seed) # Set up random number generator

	### Read data
	# G is the M x N matrix of genotypes (allele counts, 9 for missing), with
	# the SNPs in random order (s_ord), and q_nrm is the number of
	# non-missing alleles of each individual
	if args.plink:
		assert os.path.isfile(f"{args.bfile}.bed"), "bed file doesn't exist!"
		assert os.path.isfile(f"{args.bfile}.bim"), "bim file doesn't exist!"
		assert os.path.isfile(f"{args.bfile}.fam"), "fam file doesn't exist!"
		print("Reading data...", end="", flush=True)
		G, q_nrm, s_ord, M, N = utils.readPlink(args.bfile, rng)
		print(f"\rLoaded {N} samples and {M} SNPs.")
		assert not np.any(q_nrm == 0), "Sample(s) with zero information!"
	else:
		# the backing file of a gen_tibble has one byte per genotype (0, 1
		# and 2 alternate alleles, 3 for missing, larger values for imputed
		# genotypes), without a header. It holds the whole FBM, of which the
		# gen_tibble uses some rows (individuals) and columns (loci), in its
		# own order (R passes length-one vectors as scalars, hence
		# atleast_1d)
		rows = np.atleast_1d(np.asarray(args.rows, dtype=np.intp))
		cols = np.atleast_1d(np.asarray(args.cols, dtype=np.intp))
		N = rows.shape[0]
		M = cols.shape[0]
		fbm_nrow, fbm_ncol = (int(d) for d in args.fbm_dim)
		assert os.path.getsize(args.bfile) == fbm_nrow*fbm_ncol, \
			"Backing file doesn't match the dimensions of the FBM!"
		# the FBM is column-major (individuals x loci), i.e. a row-major
		# loci x individuals array: select the loci (rows of this array) and
		# the individuals (columns) of the gen_tibble, reading only those
		# from the memory-mapped file
		G_fbm = np.memmap(args.bfile, dtype=np.uint8, mode="r",
			shape=(fbm_ncol, fbm_nrow))
		G_raw = G_fbm[np.ix_(cols, rows)]
		del G_fbm
		# shuffle SNPs as utils.readPlink does, so that mini-batches are random
		s_ord = np.arange(M, dtype=np.uint32)
		rng.shuffle(s_ord)
		G = np.empty((M, N), dtype=np.uint8)
		# row s_ord[j] of G holds SNP j, the same permutation applied by
		# shared.expandShuf() to PLINK files, so that shared.reorderP()
		# restores the original order of P below
		G[s_ord] = G_raw
		del G_raw
		# set missing/imputed values to 9
		G[G>=3] = 9
		# q_nrm is the number of non-missing alleles per individual
		q_nrm = 2.0*np.sum(G != 9, axis=0, dtype=np.double)
		assert not np.any(q_nrm == 0), "Sample(s) with zero information!"

	if int(M/args.batches) < 10000:
		print("\nWARNING: Few SNPs per mini-batch!\n")
		
	# Set up parameters
	if args.supervised is not None: # Supervised mode
		# Check input of ancestral sources
		print("Ancestry estimation in supervised mode!")
		y = np.loadtxt(args.supervised, dtype=np.uint8).reshape(-1)
		assert y.shape[0] == N, "Number of samples differ between files!"
		assert np.max(y) <= args.K, "Wrong number of ancestral sources!"
		assert np.min(y) >= 0, "Wrong format in population assignments!"
		print(f"{np.sum(y > 0)}/{N} individuals with fixed ancestry.")

		# Count ancestral sources
		z = np.unique(y[y > 0])
		z -= 1
		z = np.sort(z)

		# Initialize parameters
		P = rng.random(size=(M, args.K))
		Q = rng.random(size=(N, args.K))
		P[:,z] = 0.0
		shared.initP(G, P, y)
		shared.initQ(Q, y)
		del z
	elif args.projection is not None: # Projection mode
		# Check input of ancestral allele frequencies
		print("Ancestry estimation in projection mode!")
		P_raw = np.loadtxt(args.projection, dtype=float).clip(min=1e-5, max=1.0-(1e-5))
		assert P_raw.shape[0] == M, "Number of SNPs differ between files!"
		assert P_raw.shape[1] == args.K, "Wrong number of ancestral sources!"
		
		# Shuffle SNPs according to genotype matrix
		P = np.zeros_like(P_raw)
		shared.shuffleP(P_raw, P, s_ord)
		del P_raw

		# Initialize Q matrix
		Q = rng.random(size=(N, args.K)).clip(min=1e-5, max=1.0-(1e-5))
		Q /= np.sum(Q, axis=1, keepdims=True)
	else: # Unsupervised mode
		if args.random_init: # Random initialization
			print("Random initialization.")
			P = rng.random(size=(M, args.K)).clip(min=1e-5, max=1.0-(1e-5))
			Q = rng.random(size=(N, args.K)).clip(min=1e-5, max=1.0-(1e-5))
			Q /= np.sum(Q, axis=1, keepdims=True)
		else: # SVD-based initialization
			f = np.zeros(M, dtype=np.float32)
			shared.estimateFreq(G, f)
			# (the message of fastmixture, "Please perform MAF filtering!", is
			# expanded to explain the options in R)
			assert (np.min(f) > 0.0) & (np.max(f) < 1.0), \
				"Please perform MAF filtering! The SVD-based initialisation " \
				"(random_init = FALSE) needs all loci to be polymorphic: " \
				"remove the monomorphic loci, or use random_init = TRUE"

			# Initialize P and Q matrices from SVD and ALS
			print("SVD/ALS initialization.", end="", flush=True)
			ts = time()
			if (args.subsample < 1.0) and (args.min_subsample < M): # Subsampling mode
				B = int(max(args.min_subsample, min(M*args.subsample, args.max_subsample)))
				U_sub, S, V = utils.randomSVD(G, f, args.K-1, B, args.chunk, args.power, rng)
				U_rem = utils.projectSVD(G, S, V, f, B, args.chunk)
				P, Q = utils.factorSub(U_sub, U_rem, S, V, f, args.als_iter, args.als_tole, rng)
				del U_sub, U_rem
			else: # Standard mode
				U, S, V = utils.randomSVD(G, f, args.K-1, M, args.chunk, args.power, rng)
				P, Q = utils.factorALS(U, S, V, f, args.als_iter, args.als_tole, rng)
				del U
			print(f"\rSVD/ALS initialization.\t\t({time()-ts:.1f}s)")
			del S, V, f
		y = None

	# Run options dictionary
	run = {
		"iter":args.iter,
		"tole":args.tole,
		"check":args.check,
		"batches":args.batches,
	}

	# Fastmixture
	if args.projection is not None: # Projection mode
		from fastmixture import projection
		res = projection.fastRun(G, P, Q, q_nrm, rng, run)
	else: # Unsupervised or supervised mode
		from fastmixture import functions
		res = functions.fastRun(G, P, Q, q_nrm, y, rng, run)

		# Cross-validation mode in unsupervised mode
		if args.cv is not None:
			from fastmixture import cross
			run["tole"] = args.cv_tole
			run["cross"] = args.cv
			res_crv = cross.crossRun(G, P, Q, q_nrm, rng, run)

	# Print elapsed time for estimation
	t_tot = time()-start
	t_min = int(t_tot//60)
	t_sec = int(t_tot - t_min*60)
	print(f"Total elapsed time: {t_min}m{t_sec}s")
	
	# return Q, followed by P (unless no_freqs) and the cross-validation
	# results (a dict with the average and standard deviation of the error
	# across folds) if cv was requested; reticulate converts the tuple to an
	# unnamed list
	if not args.no_freqs:
		# Return P in the original SNP order
		if args.projection is None:
			P_ord = np.zeros_like(P)
			shared.reorderP(P, P_ord, s_ord)
			P = P_ord
		if args.cv is not None:
			return Q, P, res_crv
		else:
			return Q, P
	else:
		if args.cv is not None:
			return Q, res_crv
		else:
			return Q


	# the output files of the command line version, kept for reference
	# ### Save estimates and write output to log-file
	# np.savetxt(f"{args.out}.K{args.K}.s{args.seed}.Q", Q, fmt="%.6f")
	# print(f"Saved Q matrix as {args.out}.K{args.K}.s{args.seed}.Q")
	# if not args.no_freqs and (args.projection is None): # Save frequencies
	# 	P_ord = np.zeros_like(P)
	# 	shared.reorderP(P, P_ord, s_ord)
	# 	np.savetxt(f"{args.out}.K{args.K}.s{args.seed}.P", P_ord, fmt="%.6f")
	# 	print(f"Saved P matrix as {args.out}.K{args.K}.s{args.seed}.P")
	# 
	# # Write to log-file
	# with open(f"{args.out}.K{args.K}.s{args.seed}.log", "a") as log:
	# 	log.write(f"\nFinal log-likelihood: {res["like"]:.1f}\n")
	# 	if res["conv"]:
	# 		log.write(f"Converged in {res["iter"]+1} iterations.\n")
	# 	else:
	# 		log.write(f"EM algorithm did not converge in {args.iter} iterations!\n")
	# 	log.write(f"Total elapsed time: {t_min}m{t_sec}s\n")
	# 	log.write(f"Saved Q matrix as {args.out}.K{args.K}.s{args.seed}.Q\n")
	# 	if not args.no_freqs and (args.projection is None):
	# 		log.write(f"Saved P matrix as {args.out}.K{args.K}.s{args.seed}.P\n")



##### Main exception #####
assert __name__ != "__main__", "Please use the 'fastmixture' command!"

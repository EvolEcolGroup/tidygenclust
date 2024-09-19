#' fastmixture algorithm for population genetics clustering
#'
#' This function implements the fastmixture algorithm for population genetics clustering
#' by calling the relevant python module.
#'
#' @param bfile is the name of the binary plink file (without the .bed extension)
#' @param k the number of ancestral components (clusters)
#' @param threads the number of threads to use (1)
#' @param seed the random seed (42)
#' @param outprefix the prefix of the output files (fastmixture)
#' @param iter the maximum number of iterations (1000)
#' @param tole the tolerance in log-likelihood units between iterations (0.5)
#' @param batches the number of maximum mini-batches (32)
#' @param supervised the name fo the file with the supervised labels (NULL)
#' @param check the number of iterations to check for convergence (5)
#' @param power number of power iterations in randomised SVD (11)
#' @param chunk the number of SPs in chunk operations (8192)
#' @param als_iter the maximum number of iterations in the ALS algorithm (1000)
#' @param als_tole the tolerance for the RMSE of P between iterations (1e-4)
#' @param no_freqs do not save P-matrix (TRUE)
#' @param random_init random initialisation of parameters (TRUE)
#' @return NULL (for the moment, in the future a list)
#' @export

fastmixture <- function(bfile, k, threads=1, seed=42,
                        outprefix="fastmixture", iter=1000, tole=0.5,
                        batches=32, supervised=NULL, check=5, power=11,
                        chunk=8192, als_iter=1000, als_tole=1e-4,
                        no_freqs=TRUE, random_init=TRUE) {

  # create a namespace object with all the inputs
  argparse <- reticulate::import("argparse")
  rfastmixture_args <- argparse$Namespace(bfile = bfile, K = as.integer(k), threads = as.integer(threads),
                                          seed = as.integer(seed),
                                    out = outprefix, iter = as.integer(iter), tole = tole,
                                    batches = as.integer(batches), supervised = supervised, check = as.integer(check),
                                    power = as.integer(power), chunk = as.integer(chunk), als_iter = as.integer(als_iter),
                                    als_tole = als_tole, no_freqs = no_freqs,
                                    random_init = random_init)
  fastmixture_res<-.py_rfastmixture$fastmixture_run(args = rfastmixture_args)
 # system(paste("fastmixture --bfile", bfile, "--k", k, "--threads", threads, "--seed", seed, "--outprefix", outprefix, "--iter", iter, "--tole", tole, "--batches", batches, "--supervised", supervised, "--check", check, "--power", power, "--chunk", chunk, "--als_iter", als_iter, "--als_tole", als_tole, "--no_frequs", no_frequs, "--random_init", random_init))
}

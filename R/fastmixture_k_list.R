#' Run fastmixture for multiple K values
#'
#' This function runs the `fastmixture` algorithm for a range of K values and returns
#' a list of `q_matrix` objects for each K.
#'
#' @param gt_object A gentibble object that will be analyzed by the fastmixture algorithm.
#' @param k_values A vector of integers specifying the different K values for which to run the `fastmixture`.
#' @param threads Number of threads to use for parallelization. Default is 1.
#' @param seed Random seed for reproducibility. Default is 42.
#' @param iter Number of iterations for the fastmixture algorithm. Default is 1000.
#' @param tole Convergence tolerance value for the fastmixture algorithm. Default is 0.5.
#' @param batches Number of batches for the fastmixture algorithm. Default is 32.
#' @param supervised Supervised parameter for fastmixture. Default is NULL.
#' @param check Number of checks to be performed during the process. Default is 5.
#' @param power Power parameter used for the fastmixture algorithm. Default is 11.
#' @param chunk Size of data chunk used in processing. Default is 8192.
#' @param als_iter Number of ALS (alternating least squares) iterations. Default is 1000.
#' @param als_tole Tolerance value for ALS convergence. Default is 1e-4.
#' @param no_freqs Logical parameter to indicate whether to use allele frequencies. Default is TRUE.
#' @param random_init Logical parameter to indicate whether to use random initialization. Default is TRUE.
#' @return A list of `q_matrix` objects, one for each K value provided.
#' @examples
#' \dontrun{
#' k_values <- c(2, 3, 4)
#' result <- fastmixture_k_list(agt_object, k_values)
#' }
#' @export

fastmixture_k_list <- function(gt_object, k_values,
                               threads = 1, seed = 42, 
                               iter = 1000, tole = 0.5, 
                               batches = 32, supervised = NULL, 
                               check = 5, power = 11, 
                               chunk = 8192, als_iter = 1000, 
                               als_tole = 1e-4, no_freqs = TRUE, 
                               random_init = TRUE) {
  
  # Run fastmixture for each K in k_values and store the results in a list
  qmat_list <- lapply(k_values, function(k) {
    fastmixture(gt_object, k = k, threads = threads, seed = seed, 
                outprefix = paste0("fastmixture_k", k), iter = iter, tole = tole, 
                batches = batches, supervised = supervised, check = check, 
                power = power, chunk = chunk, als_iter = als_iter, als_tole = als_tole, 
                no_freqs = no_freqs, random_init = random_init)
  })
  
  # Return the list of q_matrix objects
  return(qmat_list)
}

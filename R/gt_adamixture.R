#' Run Adamixture clustering
#'
#' Run the Adamixture ancestry estimation algorithm and return results as a
#' `gt_admix` object.
#'
#' This function provides an R interface to the Python implementation of
#' Adamixture via `reticulate`.
#'
#' The function follows the same conventions and output structure as
#' [gt_fastmixture()], allowing interchangeable downstream analysis and
#' visualisation.
#'
#' @param x A `gen_tibble`.
#' @param k Integer vector of K values.
#' @param n_runs Number of replicate runs per K.
#' @param threads Number of CPU threads.
#' @param seed Random seed.
#' @param lr Adam learning rate.
#' @param beta1 Adam β₁ parameter.
#' @param beta2 Adam β₂ parameter.
#' @param reg_adam Adam epsilon parameter.
#' @param lr_decay Learning-rate decay factor.
#' @param min_lr Minimum learning rate.
#' @param patience_adam Number of evaluations without improvement before
#'   reducing the learning rate.
#' @param tol_adam Convergence tolerance.
#' @param max_iter Maximum number of Adam-EM iterations.
#' @param check Frequency of likelihood evaluation.
#' @param chunk_size Number of SNPs processed per chunk.
#'
#' @return A `gt_admix` object.
#'
#' @details
#' Adamixture is a maximum-likelihood ancestry estimation method that uses
#' an Adam-EM optimisation procedure.
#'
#' Results are returned in the standard `gt_admix` format used throughout
#' tidygenclust.
#'
#' @export
gt_adamixture <- function(
    x,
    k,
    n_runs = 1,
    threads = 1,
    seed = 42,
    lr = 0.005,
    beta1 = 0.80,
    beta2 = 0.88,
    reg_adam = 1e-8,
    lr_decay = 0.5,
    min_lr = 1e-4,
    patience_adam = 3,
    tol_adam = 0.1,
    max_iter = 10000,
    check = 5,
    chunk_size = 4096,
    python_env = NULL
) {
  
  # ---------------------------------------------------------------------------
  # Checks
  # ---------------------------------------------------------------------------
  
  if (!inherits(x, "gen_tibble")) {
    stop("x must be a gen_tibble.")
  }
  
  if (!is.numeric(k) || length(k) < 1) {
    stop("k must contain at least one value.")
  }
  
  k <- as.integer(k)
  
  bfile <- tidypopgen::gt_get_file_names(x)[2]
  n_indiv <- nrow(x)
  n_loci <- nrow(tidypopgen::show_loci(x))
  
  

  # ---------------------------------------------------------------------------
  # Load module
  # ---------------------------------------------------------------------------
  # TODO this could be moved to on load likein gt_fastmixture.R
  adamixture <- reticulate::import(
    "adamixture",
    delay_load = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # get metadata from the gentibble
  # ---------------------------------------------------------------------------
  
  ids <- NULL
  groups <- NULL
  
  if ("id" %in% names(x)) {
    ids <- x$id
  }
  
  if ("group" %in% names(x)) {
    groups <- x$group
  }
  
  # ---------------------------------------------------------------------------
  # Storage
  # ---------------------------------------------------------------------------
  
  results_k <- list()
  results_Q <- list()
  results_P <- list()
  results_loglik <- list()
  results_log <- list()
  
  idx <- 1
  
  # ---------------------------------------------------------------------------
  # Main loop
  # ---------------------------------------------------------------------------
  
  for (k_i in k) {
    
    for (run_i in seq_len(n_runs)) {
      
      run_seed <- seed + run_i - 1
      
      # -----------------------------------------------------------------------
      # Build argument object expected by .py_adamixture()
      # -----------------------------------------------------------------------
      
      args <- list(
        X = bfile,
        K = as.integer(k_i),
        
        lr = lr,
        beta1 = beta1,
        beta2 = beta2,
        reg_adam = reg_adam,
        
        lr_decay = lr_decay,
        min_lr = min_lr,
        patience_adam = as.integer(patience_adam),
        
        tol_adam = tol_adam,
        max_iter = as.integer(max_iter),
        check = as.integer(check),
        
        chunk_size = as.integer(chunk_size),
        
        threads = as.integer(threads),
        seed = as.integer(run_seed),
        n_indiv = n_indiv,
        n_loci = n_loci
      )
      
      # -----------------------------------------------------------------------
      # Run Adamixture
      # -----------------------------------------------------------------------
      
      fit <- adamixture$.py_adamixture(args)
      
      fit <- reticulate::py_to_r(fit)
      
      # -----------------------------------------------------------------------
      # Store results
      # -----------------------------------------------------------------------
      
      results_k[[idx]] <- k_i
      
      results_Q[[idx]] <- fit$Q
      
      results_P[[idx]] <- fit$P
      
      results_loglik[[idx]] <- fit$loglik
      
      results_log[[idx]] <- fit$log
      
      idx <- idx + 1
    }
  }
  
  # ---------------------------------------------------------------------------
  # Build output
  # ---------------------------------------------------------------------------
  
  structure(
    list(
      k = results_k,
      Q = results_Q,
      P = results_P,
      log = results_log,
      loglik = results_loglik,
      id = ids,
      group = groups
    ),
    class = "gt_admix"
  )
}
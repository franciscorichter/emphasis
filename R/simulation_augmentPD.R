#' Augment Phylogenetic Diversity using MCEM
#'
#' This function performs augmentation of phylogenetic diversity (PD) on a given phylogenetic tree using the Monte Carlo Expectation-Maximization (MCEM) algorithm. It aims to estimate the maximum likelihood parameters for the given phylo object based on the provided parameters and constraints.
#'
#' @param phylo A phylogenetic tree of class 'phylo'.
#' @param pars A numeric vector of initial parameters for the MCEM algorithm.
#' @param maxN The maximum number of iterations for the MCEM algorithm.
#' @param max_missing The maximum amount of missing data allowed.
#' @param lower_bound The lower bound constraint for parameter optimization.
#' @param upper_bound The upper bound constraint for parameter optimization.
#' @param num_threads The number of threads to be used for parallel computation. Defaults to 1.
#' @param sample_size The size of the sample to be used in the MCEM algorithm. Defaults to 1.
#' @param soc The second order correction for the MCEM algorithm. Defaults to 2.
#' @param xtol_rel The relative tolerance for convergence in the optimization algorithm. Defaults to 0.00001.
#' @param verbose A logical value indicating if progress messages should be printed. Defaults to FALSE.
#'
#' @details The function uses the \code{\link[ape:branching.times]{branching.times}} function from the 'ape' package to calculate the branching times of the phylogenetic tree. It then calls a C++ function through the .Call interface for the MCEM algorithm. The parameters `max_lambda`, `brts`, and `result` are used internally within the function.
#'
#' @return The result of the MCEM algorithm, which could be parameters estimates, log-likelihood values, or other relevant metrics depending on the implementation of the C++ function `_emphasis_rcpp_mce`.
#'
#' @examples
#' \dontrun{
#'   data(bird.orders) # assuming bird.orders is a phylo object
#'   initial_pars <- c(0.1, 0.1)
#'   result <- augmentPD(phylo = bird.orders, pars = initial_pars, 
#'                       maxN = 1000, max_missing = 0.1, 
#'                       lower_bound = 0.001, upper_bound = 10)
#' }
#' 
#' @export
augmentPD <- function(phylo, pars, maxN, max_missing, lower_bound, upper_bound, 
                      num_threads = 0, sample_size = 1, soc = 2, xtol_rel = 0.00001, 
                      verbose = FALSE) {
  # Validate inputs
  if(!inherits(phylo, "phylo")) stop("The 'phylo' argument must be of class 'phylo'.")

  # Ensure numeric inputs are all numeric
  if(!all(sapply(list(maxN, max_missing, lower_bound, upper_bound, num_threads, sample_size, soc, xtol_rel), is.numeric))) {
    stop("All parameters defining numerical values must be of type numeric.")
  }
  
  max_lambda = max_missing * 100
  brts = ape::branching.times(phylo)
  
  result <- .Call('_emphasis_rcpp_mce', PACKAGE = 'emphasis', brts, pars, 
                  sample_size, maxN, soc, max_missing, max_lambda, 
                  lower_bound, upper_bound, xtol_rel, num_threads)
  
  if(verbose) {
    cat("Augmentation completed.\n")
  }
  
  return(result)
}

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
                      num_threads = 1, sample_size = 1, soc = 2, xtol_rel = 0.00001,
                      verbose = FALSE) {
  if(!inherits(phylo, "phylo")) stop("The 'phylo' argument must be of class 'phylo'.")

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

#' Augment Multiple Phylogenetic Trees with Parameter Diversity
#'
#' This function simulates and augments phylogenetic data across multiple trees
#' using parameter sampling and the augmentPD function.
#'
#' @param phylo A phylogenetic tree object.
#' @param n_trees Number of trees to simulate.
#' @param mu_interval Interval for sampling mu parameter.
#' @param lambda_interval Interval for sampling lambda parameter.
#' @param betaN_interval Interval for sampling betaN parameter.
#' @param betaP_interval Interval for sampling betaP parameter.
#' @param max_lin Maximum number of lineages (default 1e+6).
#' @param max_tries Maximum number of tries for augmentation (default 1).
#' @return A list containing the generated trees, parameters, rejection reasons,
#'         timings, and log-likelihood estimations.
#' @examples
#' # Example usage
#' result <- AugmentMultiplePhyloPD(my_phylo, 10, c(0.1, 0.5), c(0.2, 0.6),
#'                                  c(0.3, 0.7), c(0.4, 0.8))
#' @export
AugmentMultiplePhyloPD <- function(phylo,
                           n_trees,
                           mu_interval,
                           lambda_interval,
                           betaN_interval,
                           betaP_interval,
                           max_lin = 1e+6,
                           max_tries = 1){

  trees <- vector("list", n_trees)
  rejected_overruns <- vector("list", n_trees)
  rejected_lambda <- vector("list", n_trees)
  rejected_zero_weights <- vector("list", n_trees)
  times <- vector("list", n_trees)
  loglik_estimation <- vector("list", n_trees)
  logf <- vector("list", n_trees)
  logg <- vector("list", n_trees)

  name.param <- c("mu","lambda", "betaN","betaP")
  true.param <- vector(mode='list', length=4)
  names(true.param) <- name.param

  for(j in 1:n_trees){

    lambda_sample <- runif(1, lambda_interval[1], lambda_interval[2])
    mu_sample <- runif(1, mu_interval[1], mu_interval[2])
    betaN_sample <- runif(1, betaN_interval[1], betaN_interval[2])
    betaP_sample <- runif(1, betaP_interval[1], betaP_interval[2])
    sim.param <- c(mu_sample, lambda_sample, betaN_sample, betaP_sample)

    outputs <- try({
      augmentPD(phylo = phylo,
            pars = sim.param,
            maxN = max_tries,
            max_missing = max_lin,
            lower_bound = -1+c(lambda_interval[1],
                            mu_interval[1],betaN_interval[1],betaP_interval[1]),
            upper_bound = 1+c(lambda_interval[2],mu_interval[2],betaN_interval[2],
                            betaP_interval[2]))
    },silent = TRUE)

    if (!inherits(outputs, "try-error")) {
      trees[[j]] <- outputs$trees
      rejected_overruns[[j]] <- outputs$rejected_overruns
      rejected_lambda[[j]] <- outputs$rejected_lambda
      rejected_zero_weights[[j]] <- outputs$rejected_zero_weights
      times[[j]] <- outputs$time
      loglik_estimation[[j]] <- outputs$fhat
      logf[[j]] <- outputs$logf
      logg[[j]] <- outputs$logg

      for (i in 1:4) {
        true.param[[i]] <- c(true.param[[i]], sim.param[i])
      }
    }else{
      trees[[j]] <- outputs
      rejected_overruns[[j]] <- outputs
      rejected_lambda[[j]] <- outputs
      rejected_zero_weights[[j]] <- outputs
      times[[j]] <- outputs
      loglik_estimation[[j]] <- outputs
      logf[[j]] <- outputs
      logg[[j]] <- outputs
    }

    message(sprintf("Progress: %d/%d", j, n_trees))
  }

  results = list(
    trees = trees,
    param = true.param,
    rejected_overruns = rejected_overruns,
    rejected_lambda = rejected_lambda,
    rejected_zero_weights = rejected_zero_weights,
    times = times,
    loglik_estimation = loglik_estimation,
    logf = logf,
    logg = logg
  )
  return(results)
}

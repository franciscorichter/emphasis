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

  # Preallocate lists
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
    
    # Randomly sample parameter values
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
    
    if (class(outputs) != "try-error") {
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
    
    # Print progress
    svMisc::progress(j, n_trees, progress.bar = TRUE, init = (j == 1))
  }
  
  # Package and return results
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

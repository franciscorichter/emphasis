AugmentMultiplePhyloPD <- function(phylo,
                           n_trees,
                           mu_interval,
                           lambda_interval,
                           betaN_interval,
                           betaP_interval,
                           max_lin = 1e+6,
                           max_tries = 1){
  
  brts = sort(ape::branching.times(phylo))
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
      augmentPD(brts = brts,
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

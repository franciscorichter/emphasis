#' Train a Generalized Additive Model on Simulation Data
#'
#' This function fits a Generalized Additive Model (GAM) to the provided simulation
#' results. It models the log-likelihood of simulation success (`loglik`) as a function
#' of the simulation parameters, using smooth terms for each parameter. Only simulations
#' that completed successfully are used for the model fitting.
#' 
#' @param results A list containing simulation results with the following components:
#'   - `loglik_estimation`: A list of log-likelihood values or `try-error` objects.
#'   - `trees`: A list of data frames, each containing tree data.
#'   - `param`: A list of named vectors containing the simulation parameters `mu`,
#'     `lambda`, `betaN`, and `betaP`.
#' 
#' @return A `gam` object representing the fitted model.
#' @export
#'
#' @examples
#' # Assuming 'results' is your list of simulation results:
#' gam_model <- train_GAM(results)
#' summary(gam_model)
#'
train_GAM <- function(results){
  
  count_extinct_species <- function(tree) {
    extinction_times = as.data.frame(tree)["t_ext"]
    sum(extinction_times < 1e+10)
  }
  
  completed_indices <- which(!sapply(results$loglik_estimation, inherits, "try-error"))
  trees = results$trees[completed_indices]
  num_extinct_per_tree <- sapply(trees, count_extinct_species)
  
  sim_data <- data.frame(mu = results$param$mu,
                         lambda = results$param$lambda,
                         betaN = results$param$betaN,
                         betaP = results$param$betaP,
                         loglik = unlist(results$loglik_estimation[completed_indices]),
                         nspecies =  num_extinct_per_tree
  )
  
  
  cat("Training GAM...")
  gam_loglik = mgcv::gam(loglik~s(mu)+s(lambda)+s(betaN)+s(betaP),data=sim_data,family= mgcv::scat(link="identity"))
  
  # Check if all parameters are significant
  pvals <- summary(gam_loglik)$p.values
  significant <- pvals < 0.05
  if(all(significant)) {
    cat("All parameters are significant.\n")
  } else {
    cat("Not all parameters are significant.\n")
  }
  
  
  return(gam_loglik)
}

#' Optimize GAM Model Parameters with Fixed betaP
#'
#' This function optimizes the parameters of a given Generalized Additive Model (GAM) 
#' with the constraint that betaP is fixed at zero. It applies penalties to discourage 
#' selecting parameter values near the boundaries of their respective intervals, aiming 
#' to find optimal values of mu, lambda, and betaN that minimize the estimated log likelihood.
#'
#' @param pgam A GAM model object for which the log likelihood is to be optimized. 
#'             The model should include mu, lambda, betaN, and betaP as parameters.
#'
#' @return A list containing the optimized parameter values (pars) and the minimized 
#'         log likelihood value (ml).
#'
#' @examples
#' # Assuming you have a GAM model 'pgam' fitted with parameters including mu, lambda, betaN, betaP
#' result <- dd_ML_est(pgam)
#' print(result$pars)  # Optimized parameters
#' print(result$ml)    # Minimized log likelihood
#'
#' @export
#'
#' @importFrom stats optim
#' @importFrom mgcv predict
dd_ML_est <- function(pgam){
  
  # Adjusted objective function with betaP fixed at zero
  objective_function_fixed_betaP <- function(params) {
    # Extract the parameters, but fix betaP to 0
    mu <- params[2]
    lambda <- params[1]
    betaN <- (mu-lambda)/params[3]
    betaP <- 0  # Fixed
    
    
    betaN_interval <- range(pgam$model$betaN)
    betaP_interval <-  range(pgam$model$betaP)
    lambda_interval <-  range(pgam$model$lambda)
    mu_interval <-  range(pgam$model$mu)
    
    # Penalty factors for each parameter near the boundaries
    penalty_factor <- -1000  # Adjust the penalty factor based on the scale of your log likelihood values
    
    # Calculate penalty for being near the boundary of each parameter
    # Assuming mu_interval, lambda_interval, and betaN_interval are defined globally
    penalty <- 0
    if (mu <= mu_interval[1] + (mu_interval[2] - mu_interval[1]) * 0.05 || mu >= mu_interval[2] - (mu_interval[2] - mu_interval[1]) * 0.05) {
      penalty <- penalty + penalty_factor
    }
    if (lambda <= lambda_interval[1] + (lambda_interval[2] - lambda_interval[1]) * 0.05 || lambda >= lambda_interval[2] - (lambda_interval[2] - lambda_interval[1]) * 0.05) {
      penalty <- penalty + penalty_factor
    }
    # Adding penalty for betaN based on its interval
    if (betaN <= betaN_interval[1] + (betaN_interval[2] - betaN_interval[1]) * 0.05 || betaN >= betaN_interval[2] - (betaN_interval[2] - betaN_interval[1]) * 0.05) {
      penalty <- penalty + penalty_factor
    }
    
    # Create a data frame with these parameters
    newdata <- data.frame(mu = mu, lambda = lambda, betaN = betaN, betaP = betaP)
    
    # Use predict() to get the estimated log likelihood
    loglik_pred <- predict(pgam, newdata = newdata, type = "response")
    
    # Return the negative log likelihood
    return(loglik_pred - penalty)
  }
  
  # Initial values for the parameters
  init_params_fixed_betaP <- c(0.8, 0.1, 120)
  
  # Adjust optim to use the new objective function and initial parameters
  # Also adjust the bounds to exclude betaP
  optim_results_fixed_betaP <- optim(par = init_params_fixed_betaP, fn = objective_function_fixed_betaP, 
                                     method = "L-BFGS-B", 
                                     lower = c(0, 0, 60), 
                                     upper = c(1, 0.5, 1000))
  
  
  return(list(pars=optim_results_fixed_betaP$par,ml = optim_results_fixed_betaP$value))
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

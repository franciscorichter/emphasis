#' Generate Phylogenetic Diversity Data
#'
#' This function simulates phylogenetic trees based on given parameters.
#' @param n_trees Number of trees to simulate.
#' @param mu_interval Interval for mu parameter.
#' @param lambda_interval Interval for lambda parameter.
#' @param betaN_interval Interval for betaN parameter.
#' @param betaP_interval Interval for betaP parameter.
#' @param max_lin Maximum number of lineages.
#' @param max_tries Maximum number of tries for simulation.
#' @return A list containing the simulated trees and parameters.
#' @examples
#' generatePhyloPD(n_trees = 10, 
#'                mu_interval = c(0.1, 0.5), 
#'                lambda_interval = c(0.1, 0.5),
#'                betaN_interval = c(0.1, 0.5), 
#'                betaP_interval = c(0.1, 0.5),
#'                max_lin = 1e6,
#'                max_tries = 1)
#' @export
generatePhyloPD <- function(n_trees,
                            mu_interval,
                            lambda_interval,
                            betaN_interval,
                            betaP_interval,
                            max_lin = 1e+6,
                            max_tries = 1){
  
  # Preallocate lists
  trees <- vector("list", n_trees)
  extrees <- vector("list", n_trees)
  Lmats <- vector("list", n_trees)
  brds_s <- vector("list", n_trees)
  
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
    
    outputs <- tryCatch({
      sim_tree_pd_cpp(pars = sim.param, 
                      max_t = 1, 
                      max_lin = max_lin, 
                      max_tries = max_tries,useDDD = TRUE)
    }, error = function(e) NULL)
    
    if (is.list(outputs) && max(outputs$brts) == 1) {
      trees[[j]] <- outputs[[1]]
      extrees[[j]] <- outputs[[2]]
      Lmats[[j]] <- outputs[[3]]
      brds_s[[j]] <- outputs[[4]]
      for (i in 1:4) {
        true.param[[i]] <- c(true.param[[i]], sim.param[i])
      }
    }
    
    # Print progress
    svMisc::progress(j, n_trees, progress.bar = TRUE, init = (j == 1))
  }
  
  # Package and return results
  results = list(trees = trees, param = true.param, tas = extrees, L = Lmats, brts = brds_s)
  return(results)
}


generatePhyloPD_old <- function(n_trees,
                            mu_interval,
                            lambda_interval,
                            betaN_interval,
                            betaP_interval,
                            max_lin = 1e+6,
                            max_tries = 1){
  
  trees <- list()
  extrees <- list()
  Lmats <- list()
  brds_s <- list()
  
  name.param <- c("mu","lambda", "betaN","betaP") 
  true.param <- vector(mode='list', length=4)
  names(true.param) <- name.param
  j=0
  while(length(trees) < n_trees){
    j=j+1
    
    lambda_sample <- runif(1, min = lambda_interval[1], max = lambda_interval[2])
    mu_sample <- runif(1, min = mu_interval[1], max = mu_interval[2])
    betaN_sample <- runif(1, min = betaN_interval[1], max = betaN_interval[2])
    betaP_sample <- runif(1, min = betaP_interval[1], max = betaP_interval[2])

    sim.param <- c(mu_sample,lambda_sample,betaN_sample,betaP_sample)
    
    key = 1 
    while(key){
      
      outputs <- try(sim_tree_pd_cpp(pars = sim.param,
                                                max_t = 1,
                                                max_lin = max_lin,
                                                max_tries = max_tries),silent = TRUE)

      if(is.list(outputs)){
        if(max(outputs$brts)==1) key = 0
      } 
    }
    
    tree  <- outputs[[1]]
    extree <- outputs[[2]]
    Lmat  <- outputs[[3]]
    brds  <- outputs[[4]]
    
    trees <- append(trees, list(tree))                    # save tree
    extrees <- append(extrees, list(extree))                    # Additional Battery
    Lmats <- append(Lmats, list(Lmat))                    #
    brds_s <- append(brds_s, list(brds))                    #
    
    for (i in 1:4){
      true.param[[i]] <- c(true.param[[i]], sim.param[i]) # save param.
    }
    
    svMisc::progress(length(trees), n_trees, progress.bar = TRUE, # print
             init = (length(trees)==1))                   # progression
    
  }
  
  out <- list("trees"    = trees, "param"    = true.param, "tas" = extrees, "L" = Lmats,"brts"= brds_s )
  
  return(out)
}

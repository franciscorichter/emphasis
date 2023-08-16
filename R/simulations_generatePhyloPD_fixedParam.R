generatePhyloPD_fixedParam <- function(n_trees, sim.param, max_t=1, max_lin = 1e+6, max_tries = 100){
  
  trees <- list()
  extrees <- list()
  Lmats <- list()
  brds_s <- list()
  
  name.param <- c("mu","lambda", "betaN","betaP") 
  
  # Ensuring the sim.param vector has the right number of parameters
  if (length(sim.param) != 4) {
    stop("sim.param must have 4 values corresponding to mu, lambda, betaN, and betaP.")
  }
  
  j = 0
  while(length(trees) < n_trees){
    j = j + 1
    
    key = 1 
    tries = 0
    while(key && tries < max_tries){
      tries = tries + 1
      
      outputs <- try(sim_tree_pd_cpp(pars = sim.param,
                                                max_t = max_t,
                                                max_lin = max_lin,
                                                max_tries = max_tries, 
                                                useDDD=TRUE), 
                     silent = TRUE)
      
      if(is.list(outputs)){
        if(max(outputs$brts) == 1) key = 0
      } 
    }
    
    # Check if simulation was successful
    if(is.list(outputs)) {
      tree  <- outputs[[1]]
      extree <- outputs[[2]]
      Lmat  <- outputs[[3]]
      brds  <- outputs[[4]]
      
      trees <- append(trees, list(tree))
      extrees <- append(extrees, list(extree))
      Lmats <- append(Lmats, list(Lmat))
      brds_s <- append(brds_s, list(brds))
    }
    
    # Assuming 'progress' is a function available in your environment
    svMisc::progress(length(trees), n_trees, progress.bar = TRUE, init = (length(trees) == 1))
    
  }
  
  out <- list("trees" = trees, "param" = sim.param, "tas" = extrees, "L" = Lmats, "brts" = brds_s)
  
  return(out)
}

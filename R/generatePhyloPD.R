generatePhyloPD <- function(n_trees,
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

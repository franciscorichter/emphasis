generatePhyloDDD <- function(n_trees,
                             lambda_interval,
                             K_interval){
  
  trees <- list()
  extrees <- list()
  Lmats <- list()
  brds_s <- list()
  
  name.param <- c("lambda", "mu","K") 
  true.param <- vector(mode='list', length=3)
  names(true.param) <- name.param
  j=0
  while(length(trees) < n_trees){
    j=j+1
    
    lambda0_sample <- runif(1, min = lambda_interval[1], max = lambda_interval[2])
    mu_sample <- runif(1, min = 0, max = lambda0_sample)
    k_sample <- runif(1, min = K_interval[1], max = K_interval[2])
    
    beta_sample = (mu_sample-lambda0_sample)/k_sample

    vec.param <- c(lambda0_sample,mu_sample,k_sample)
    sim.param <- c(mu_sample,lambda0_sample,beta_sample,0)
    
    key = 1 
    while(key){
      
      outputs <- try(emphasis:::sim_tree_pd_cpp(pars = pars*max_t,
                          max_t = 1,
                          max_lin = max_lin,
                          max_tries = max_tries),silent = TRUE)
      #outputs <- emphasis:::sim.tree(pars = sim.param,
       #                        max_t = 1,
        #                       max_lin = 1e+6,
         #                      max_tries = 1)
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
      
    for (i in 1:3){
        true.param[[i]] <- c(true.param[[i]], vec.param[i]) # save param.
    }
      
    progress(length(trees), n_trees, progress.bar = TRUE, # print
               init = (length(trees)==1))                   # progression
      
    }
  
  out <- list("trees"    = trees, "param"    = true.param, "tas" = extrees, "L" = Lmats,"brts"= brds_s )
  
  return(out)
}


generatePhyloDDD_minimal <- function(n_trees,
                                     lambda_interval,
                                     K_interval){
  
  trees <- list()
  
  name.param <- c("lambda", "mu","K") 
  true.param <- vector(mode='list', length=3)
  names(true.param) <- name.param
  j=0
  while(length(trees) < n_trees){
    j=j+1
    
    lambda0_sample <- runif(1, min = lambda_interval[1], max = lambda_interval[2])
    mu_sample <- runif(1, min = 0, max = lambda0_sample)
    k_sample <- runif(1, min = K_interval[1], max = K_interval[2])
    
    beta_sample = (mu_sample-lambda0_sample)/k_sample
    
    vec.param <- c(lambda0_sample,mu_sample,k_sample)
    sim.param <- c(mu_sample,lambda0_sample,beta_sample,0)
    
    
    outputs <- emphasis:::sim_tree_pd_cpp(pars = sim.param,
                                          max_t = 1,
                                          max_lin = 1e+6,
                                          max_tries = 100)
    
    trees <- append(trees, outputs)                    # save tree                   #
    
    for (i in 1:3){
      true.param[[i]] <- c(true.param[[i]], vec.param[i]) # save param.
    }
    
    svMisc:::progress(length(trees), n_trees, progress.bar = TRUE, # print
                      init = (length(trees)==1))                   # progression
    
  }
  
  out <- list("trees"    = trees, "param"    = true.param)
  
  return(out)
}





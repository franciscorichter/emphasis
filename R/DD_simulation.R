generatePhyloDDD <- function(n_trees, 
                             mu_interval,
                             lambda_interval,
                             K_interval, 
                             ss_check = TRUE){
  
  trees <- list()
  extrees <- list()
  Lmats <- list()
  brds_s <- list()
  
  name.param <- c("mu", "lambda","K") 
  true.param <- vector(mode='list', length=3)
  names(true.param) <- name.param
  
  while(length(trees) < n_trees){
    
    lambda0_sample <- runif(1, min = lambda_interval[1], max = lambda_interval[2])
    mu_sample <- runif(1, min = mu_interval[1], max = mu_interval[2])
    k_sample <- runif(1, min = K_interval[1], max = K_interval[2])
    
    beta_sample = (mu_sample-lambda0_sample)/k_sample
    #crown_time_sample <- runif(1, min = age1[1], max = age1[2])
    
    vec.param <- c(mu_sample,lambda0_sample,beta_sample)
    sim.param <- c(mu_sample,lambda0_sample,beta_sample,0)
    
    #outputs <-  DDD::dd_sim(pars = sim.param, age = crown_time_sample, ddmodel = ddmodel1)
    outputs <- simulate_single_pd_tree_cpp(pars = sim.param,
                                           max_t = 1,
                                           max_lin = 1e+7,
                                           max_tries = 1)
    
    
    tree  <- outputs[[1]]
    extree  <- outputs[[2]]
    Lmat  <- outputs[[3]]
    brds  <- outputs[[4]]
    
    # Checking that summary statistics have no NA
    if (ss_check){
      ss <- get_ss(tree) # compute summary statistics
      no_NA_ss <- !any(is.na(ss)) # does SS have any NA values?
    }
    
    if (no_NA_ss || !ss_check){
      
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
  }
  
  out <- list("trees"    = trees, "param"    = true.param, "tas" = extrees, "L" = Lmats,"brts"= brds_s )
  
  return(out)
}

fname_ddd<- function(n_trees, lambda0,mu,beta,age,dd_mod){
  
  lambd_text <- ifelse(length(lambda0) == 2, paste(lambda0[1], lambda0[2], sep="-"), 
                       as.character(lambda0))
  mu_text <- ifelse(length(mu) == 2, paste(mu[1], mu[2], sep="-"), 
                    as.character(mu))
  
  k_text <- ifelse(length(beta) == 2, paste(beta[1], beta[2], sep="-"), 
                   as.character(beta))
  
  age_text <- as.character(age)
  
  dd_mod_text <- as.character(dd_mod)
  
  
  fname <- paste("nt", n_trees,"la0",lambd_text,"mu",mu_text,"beta",k_text,"age",age_text,"ddmod",dd_mod_text ,sep="-")
  
  return(fname)
  
}


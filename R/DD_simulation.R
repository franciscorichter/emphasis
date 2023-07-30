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
      outputs <- emphasis:::sim.tree(pars = sim.param,
                               max_t = 1,
                               max_lin = 1e+6,
                               max_tries = 1)
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


generate_ltt_dataframe <- function(trees, n_taxa, true.param){
  
  n_trees  <- length(trees) # number of trees 
  n_row <- ifelse(length(n_taxa) == 1, n_taxa, n_taxa[2])
  df.ltt <- data.frame("tree1" = rep(NA, n_row))
  
  df.rates <- as.data.frame(do.call(cbind, true.param))
  
  cat("Creating LTT dataframe...\n")
  
  for (i in 1:n_trees){
    tree <- trees[[i]] # get tree 
    ltt.coord <- ape::ltt.plot.coords(tree) # get ltt coordinates 
    ltt.coord <- as.data.frame(ltt.coord)
    ltt.coord.time <- ltt.coord$time
    n <- length(ltt.coord.time)
    df.ltt[1:n,paste("tree", i, sep = "")] <- ltt.coord$time
    progress(i, n_trees, progress.bar = TRUE, init = (i==1))
  }
  
  cat("\nCreating LTT dataframe... Done.")
  
  out <- list("ltt" = df.ltt, "rates" = df.rates) # function output
  
  return(out)
  
}

convert_ltt_dataframe_to_dataset <- function(df.ltt, true.param, nn_type){
  
  if (nn_type == "cnn-ltt"){
    ds.ltt <- torch::dataset(
      name <- "ltt_dataset", 
      initialize = function(df.ltt, true.param){
        
        # input
        df.ltt[is.na(df.ltt)] <- 0
        
        array.ltt <- df.ltt %>% 
          as.matrix() %>% 
          torch_tensor()
        self$x <- array.ltt
        
        # target 
        self$y <- torch_tensor(do.call(cbind, true.param)) # target
      }, 
      
      .getitem = function(i) {list(x = self$x[,i]$unsqueeze(1), y = self$y[i, ])},
      
      .length = function() {self$y$size()[[1]]}
    )
  }
  
  else{
    ds.ltt <- torch::dataset(
      name <- "ltt_dataset", 
      initialize = function(df.ltt, true.param){
        
        # input
        df.ltt[is.na(df.ltt)] <- 0
        
        array.ltt <- df.ltt %>% 
          as.matrix() %>% 
          torch_tensor()
        self$x <- array.ltt
        
        # target 
        self$y <- torch_tensor(do.call(cbind, true.param)) # target
      }, 
      
      .getitem = function(i) {list(x = self$x[,i], y = self$y[i, ])},
      
      .length = function() {self$y$size()[[1]]}
    )
  }
  
  return(ds.ltt)
}

extract_elements <- function(list_of_vectors, indices_to_extract){
  l <- as.list(do.call(cbind, list_of_vectors)[indices_to_extract,] 
               %>% as.data.frame())
  return(l)
}


compute_dim_ouput_flatten_cnn <- function(n_input, n_layer, kernel_size = 2){
  
  k <- kernel_size - 1 
  
  for (i in 1:n_layer){
    n_input <- as.integer((n_input - k)/2)
  }
  
  return(n_input)
  
}

DDD_par_est <- function(tree,cnn_ltt=cnn_ltt){
  max_nodes_rounded <- 550
  
  ltt.coord <- ape::ltt.plot.coords(tree) # get ltt coordinates 
  ltt.coord <- as.data.frame(ltt.coord)
  ltt.coord.time <- ltt.coord$time
  n <- length(ltt.coord.time)
  
  df.ltt <- data.frame("tree1" = rep(NA, max_nodes_rounded))
  df.ltt[1:n,1] <- ltt.coord$time
  
  crown_time = max(abs(df.ltt[!is.na(df.ltt)]))
  df.ltt = df.ltt / crown_time
  
  ds.ltt_new <- new_convert_ltt_dataframe_to_dataset(df.ltt,  "cnn-ltt")
  
  ds_eval  <- ds.ltt_new(df.ltt)
  
  data_loader_ltt <- ds_eval  %>% dataloader(batch_size=1, shuffle=FALSE)
  
  device<-"cpu"
  
  cnn_ltt$eval()
  n_out<-3
  nn.pred <- vector(mode = "list", length = n_out)
  names(nn.pred) <- names(c("lambda0", "mu"    ,  "K" ))
  
  p_dropout <- 0.01
  
  b = data_loader_ltt
  
  
  coro::loop(for (b in data_loader_ltt) {
    out <- cnn_ltt(b$x$to(device = device))
    pred <- as.numeric(out$to(device = "cpu")) # move the tensor to CPU 
    for (i in 1:n_out){nn.pred[[i]] <- c(nn.pred[[i]], pred[i])}
  })
  
  par_estim = pred/crown_time
  
  return(par_estim)
}

new_generate_ltt_dataframe<-
  function(trees, n_taxa){
    
    n_trees  <- length(trees) # number of trees 
    n_row <- ifelse(length(n_taxa) == 1, n_taxa, n_taxa[2])
    df.ltt <- data.frame("tree1" = rep(NA, n_row))
    
    #df.rates <- as.data.frame(do.call(cbind, true.param))
    
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

    out <- df.ltt # function output
    
    return(out)
    
  }

new_convert_ltt_dataframe_to_dataset<-function(df.ltt, nn_type){
  
  if (nn_type == "cnn-ltt"){
    ds.ltt <- torch::dataset(
      name <- "ltt_dataset", 
      initialize = function(df.ltt){
        
        # input
        df.ltt[is.na(df.ltt)] <- 0
        
        array.ltt <- df.ltt %>% 
          as.matrix() %>% 
          torch_tensor()
        self$x <- array.ltt

      }, 
      
      .getitem = function(i) {list(x = self$x[,i]$unsqueeze(1))},
      
      .length = function() {self$x$size()[[2]]}
    )
  }
  
  else{
    ds.ltt <- torch::dataset(
      name <- "ltt_dataset", 
      initialize = function(df.ltt){
        
        # input
        df.ltt[is.na(df.ltt)] <- 0
        
        array.ltt <- df.ltt %>% 
          as.matrix() %>% 
          torch_tensor()
        self$x <- array.ltt

      }, 
      
      .getitem = function(i) {list(x = self$x[,i])},
      
      .length = function() {self$x$size()[[2]]}
    )
  }
  
  return(ds.ltt)
}

dm_fun_all <- function(pars,
                       input,
                       num_threads) {
  
  res <- rcpp_mce_grid(as.matrix(pars),
                       brts = input$brts,
                       sample_size = input$sample_size,
                       maxN = input$maxN,
                       soc = 2,
                       max_missing = input$max_missing,
                       max_lambda = input$max_lambda,
                       lower_bound = input$lower_bound,
                       upper_bound = input$upper_bound,
                       xtol_rel = 0.1,
                       num_threads = num_threads)

  return(res)
}



#' @keywords internal
dm_fun_ext <- function(pars, input){
  pars = as.numeric(pars)
  it <- try(e_cpp(brts = input$brts,
                  init_pars = pars,
                  sample_size = input$sample_size,
                  maxN = input$sample_size,
                  soc = 2,
                  max_missing = input$max_missing,
                  max_lambda = input$max_lambda,
                  lower_bound = input$lower_bound,
                  upper_bound = input$upper_bound,
                  xtol_rel = 0.1,                   
                  num_threads = 8),
            silent = TRUE)
  
  out <- list()
  if (!inherits(it, "try-error")) {
    
    out <- list(fhat = it$fhat,
                rej_lam = it$rejected_lambda,
                rej_ov = it$rejected_overruns,
                rej_zw = it$rejected_zero_weights,
                it = it)
  }else{
    str = it
    v1 <- unlist(gregexpr("*reasons*", str)) + 9
    v2 <- unlist(gregexpr("*lamba*", str)) - 2
    
    nlamb <- as.numeric(substr(str, v1, v2))
    
    v1 <- unlist(gregexpr("*lamba*", str)) + 7
    v2 <- unlist(gregexpr("*overruns*", str)) - 2
    
    nover <- as.numeric(substr(str, v1, v2))
    
    v1 <- unlist(gregexpr("*overruns*", str)) + 10
    v2 <- unlist(gregexpr("*zero weights*", str)) - 2
    
    zerow <- as.numeric(substr(str, v1, v2))
    out <- list(fhat = NA,
                rej_lam = nlamb,
                rej_ov = nover,
                rej_zw = zerow,
                it = it)
  }
  return(out)
}




#' @export
get_random_grid <- function(num_points,
                            lower_bound = c(0, 0, -0.4, 0),
                            upper_bound = c(1, 7, 0.01, 0)){
  num_dim = length(lower_bound)
  
  total_num <- num_dim * num_points
  
  pars <- matrix(stats::runif(total_num,
                       min = lower_bound,
                       max = upper_bound),
                 byrow = TRUE,
                 ncol = num_dim)
  
  pars <- as.data.frame(pars)
  names(pars) <- paste0("par", 1:num_dim)
  return(pars)
}

which_val <- function(val, dmval1){
  check <- FALSE
  s <- 0
  while (check == FALSE) {
    s <- s + 1 
    if (!is.null(dmval1[[s]]$fhat))  check = (dmval1[[s]]$fhat == val)
  }
  return(s)
}




#' perform emphasis analysis using DE method
#' @param brts branching times of tree to fit on
#' @param num_iterations number of iterations of the DE algorithm
#' @param num_points number of particles per iteration
#' @param max_missing maximum number of missing trees
#' @param sd_vec vector of initial values of standard deviation for perturbation
#' @param lower_bound vector of lower bound values for parameters,
#' used to populate the particles
#' @param upper_bound vector of upper bound values for parameters, used to
#' populate the particles
#' @param maxN maximum number of tries per parameter combination before giving
#' up
#' @param max_lambda maximum value of lambda
#' @param disc_prop proportion of particles retained per iteration
#' @param verbose verbose output if TRUE
#' @param num_threads number of threads
#' @export
#' @rawNamespace useDynLib(emphasis)
#' @rawNamespace import(nloptr)
#' @rawNamespace import(Rcpp)
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)
emphasis_de <- function(brts,
                        num_iterations,
                        num_points,
                        max_missing,
                        sd_vec,
                        lower_bound,
                        upper_bound,
                        maxN = 10,
                        max_lambda,
                        disc_prop = 0.5,
                        verbose = FALSE,
                        num_threads = 1) {
  
  if (length(upper_bound) != length(lower_bound)) {
    stop("lower bound and upper bound vectors need to be same length")
  }
  
  if (length(upper_bound) != length(sd_vec)) {
    stop("sd vector is not adequate length")
  }
  
  init_time <- proc.time()
  
  alpha = sd_vec / num_iterations
  
  input <- list(brts = brts,
                max_missing = max_missing,
                lower_bound = lower_bound,
                upper_bound =  upper_bound,
                sample_size = 1,
                maxN = maxN,
                max_lambda = max_lambda,
                disc_prop = disc_prop)
  # get first grid
  pars <- get_random_grid(num_points,
                          lower_bound = lower_bound,
                          upper_bound = upper_bound)
  
  # initialize
  pv <- list()
  rejl_count <- rejo_count <- NULL
  logf_DE1 <- NULL
  pars_DE1 <- NULL
  
  logf_DE1b <- NULL
  pars_DE1b <- NULL
  
  if (verbose) {
    pb <- progress::progress_bar$new(format = "Progress: [:bar] :percent", total = num_iterations)
  }
  fhatdiff <- c()
  for (k in 1:num_iterations) {
   
    pv[[k]] <- pars  # Save parameter grid
    dmval = rcpp_mce_grid(as.matrix(pars),
                          brts = input$brts,
                          sample_size = input$sample_size,
                          maxN = input$maxN,
                          soc = 2,
                          max_missing = input$max_missing,
                          max_lambda = input$max_lambda,
                          lower_bound = input$lower_bound,
                          upper_bound = input$upper_bound,
                          xtol_rel = 0.1,
                          num_threads = num_threads)
    
    # rcpp_mce_grid returns a matrix with columns:
    # 1 = fhat
    # 2 = rejected_lambda
    # 3 = rejected_overruns
    # 4 = rejected_zero_weights
    
    vals = dmval[, 1]  #as.numeric(dmval[seq(from = 1, to = length(dmval), by = 4)])
    fails <- which(is.na(vals)) # failed runs return fhat = -1
    rejl = dmval[, 2] # as.numeric(dmval[seq(from = 2, to = length(dmval), by = 4)])
    rejo = dmval[, 3] # as.numeric(dmval[seq(from = 3, to = length(dmval), by = 4)])
    
    # if there were any stopped simulations due to max_lambda constrain
    if (sum(rejl[fails]) > 0) {
      input$max_lambda = input$max_lambda + 100
      rejl_count = c(rejl_count, k)
    }
    
    # 
    if (sum(rejo[fails]) > 0) {
      input$max_missing = input$max_missing + 1000
      rejo_count = c(rejo_count, k)
    }
    
    
    # we no longer need to do the lookup stuff,
    # because dmval is now a matrix of results, not a list
    # wi <- NULL
    # for (i in 1:length(vals)) {
    # wi = c(wi, which_val(vals[i], dmval1))
    #   #wi = c(wi, which_val2(vals[i], dmval))
    # }

    # still, we remove the NA values 
    wi <- which(!is.na(vals))
    
    pars <- pars[wi, ]
    vals <- -vals[wi]
    
    # Methods 
    
    ## DE1
    minval = which.min(vals)
    logf_DE1 = c(logf_DE1, vals[minval])
    pars_DE1 = rbind(pars_DE1, as.numeric(pars[minval,]))
    
    ## DE1b
    logf_DE1b = c(logf_DE1b, mean(vals))
    pars_DE1b = rbind(pars_DE1b, colMeans(pars[1:4]))
    
    ## Saving variation in the estimations 
    sd_pars <- apply(pars, 2, stats::sd)

    fhatdiff <- c(fhatdiff, stats::sd(vals))

    # if we have more than 20% of the original particles, drop half of the particles
    if (nrow(pars) > num_points * 0.2)  pars = pars[vals < stats::quantile(vals, probs = disc_prop),]

    # Increase the number of particles until having at least the initial number of particles
    
    num_to_add <- num_points - nrow(pars) 
    indices <- sample(x = seq_len(nrow(pars)), size = num_to_add, replace = TRUE)
    pars_to_add <- pars[indices, ]
    for (index in 1:nrow(pars_to_add)) {
      for (param in 1:length(upper_bound)) {
        if (upper_bound[param] == lower_bound[param]) {
          pars_to_add[index, param] <- upper_bound[param]
        } else {
          new_val <- stats::rnorm(n = 1, mean = pars_to_add[index, param], sd = sd_vec[param])
          
          # here, we can either clamp the value, or decide to redraw it
          clamp_value <- FALSE
          if (clamp_value) {
            new_val <- min(new_val, upper_bound[param])
            new_val <- max(new_val, lower_bound[param])
          } else {
            within_bounds <- new_val < upper_bound[param] && new_val < lower_bound[param]
            while (within_bounds) {
              new_val <- stats::rnorm(n = 1, mean = pars_to_add[index, param], sd = sd_vec[param])
              within_bounds <- new_val < upper_bound[param] && new_val < lower_bound[param]
            }
          }
        
          pars_to_add[index, param] <- new_val
        }
      }
    }
    
    pars <- rbind(pars, pars_to_add)

    # Decrease variation
    sd_vec = sd_vec - alpha
    
    if (verbose) {
      cat("iteration: ", k, "sd: ", sd_pars, "\n")
      pb$tick()
    }
    
  }
  total_time = proc.time() - init_time

  out <- list("parameters" = pv,
              "time" = total_time,
              "fhatdiff" = fhatdiff)
  return(out)
}
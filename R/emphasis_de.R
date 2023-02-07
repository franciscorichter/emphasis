#' @keywords internal
dm_fun_ext <- function(pars, input){
  pars = as.numeric(pars)
  it <- try(e_cpp(brts = input$brts,
                  init_pars = pars,
                  sample_size = input$sample_size,
                  maxN = input$sample_size*100000,
                  plugin = locate_plugin("remphasisrpd5c"),
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




#' @keywords internal
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
#' @param max_lambda maximum value of lambda
#' @param disc_prop proportion of particles retained per iteration
#' @param verbose verbose output if TRUE
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
                        max_lambda,
                        disc_prop = 0.5,
                        verbose = FALSE){
  init_time <- proc.time()
  
  alpha = sd_vec / num_iterations
  
  input <- list(brts = brts,
                max_missing = max_missing,
                lower_bound = lower_bound,
                upper_bound =  upper_bound,
                sample_size = 1,
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
    pb <- progress::progress_bar$new(format = "Progress: [:bar] :percent", total = 100)
  }
  fhatdiff <- c()
  for (k in 1:num_iterations) {
    if (verbose) {
      
      sd_vals <- as.numeric(apply(pars, 2, stats::sd))
      
      message(paste0("iteration: ", k, "sd: ", sd_vals, " ", sd(vals), "\n"))
     # pb$tick()
    }
    
    pv[[k]] <- pars  # Save parameter grid
    dmval = NULL
    
    # Evaluate / Simulate  loglikelihood at every point
    while (is.null(dmval)) {
      
      if (verbose) {
        dmval1 = pbmcapply::mclapply(X = split(pars, seq(nrow(pars))), 
                                    FUN = dm_fun_ext,input = input,
                                    mc.cores = parallel::detectCores())
      } else {
        dmval1 = parallel::mclapply(X = split(pars, seq(nrow(pars))), 
                                  FUN = dm_fun_ext,input = input,
                                  mc.cores = parallel::detectCores())
      }
      
      dmval = unlist(dmval1, use.names = FALSE)
      if (is.null(dmval)) {
        message("No valid values, trying a new grid with looser
                 max_lambda and max_missing")
        input$max_missing = input$max_missing * 10
        input$max_lambda = input$max_lambda * 10
        pars <- get_random_grid(num_points, lower_bound, upper_bound)
        rejo_count = c(rejo_count,Inf)
        rejl_count = c(rejl_count,Inf)
      }
    }
    
    vals = as.numeric(dmval[seq(from = 1, to = length(dmval), by = 5)])
    rejl = as.numeric(dmval[seq(from = 2, to = length(dmval), by = 5)])
    rejo = as.numeric(dmval[seq(from = 3, to = length(dmval), by = 5)])
    
    # if there were any stopped simulations due to max_lambda constrain
    if (max(rejl) > 0) {
      input$max_lambda = input$max_lambda + 100
      rejl_count = c(rejl_count, k)
    }
    
    # 
    if (max(rejo) > 0) {
      input$max_missing = input$max_missing + 1000
      rejo_count = c(rejo_count, k)
    }
    
    wi <- NULL
    for (i in 1:length(vals)) {
       wi = c(wi, which_val(vals[i], dmval1))
    }
    
    pars = pars[wi,]
    vals = -vals 
    
    # Methods 
    
    ## DE1
    minval = which.min(vals)
    logf_DE1 = c(logf_DE1, vals[minval])
    pars_DE1 = rbind(pars_DE1, as.numeric(pars[minval,]))
    
    ## DE1b
    logf_DE1b = c(logf_DE1b, mean(vals))
    pars_DE1b = rbind(pars_DE1b, colMeans(pars[1:4]))
    
    ## Saving variation in the estimations 
    par_vars <- c("par1", "par2", "par3", "par4")
    for (par_var in par_vars) {
      var_diff_name <- paste0(par_var, "diff")
      assign(var_diff_name, c(get(var_diff_name), stats::sd(pars[[par_var]])))
    }

    fhatdiff <- c(fhatdiff, stats::sd(vals))

    # if we have more than 20% of the original particles, drop half of the particles
    if (nrow(pars) > num_points * 0.2)  pars = pars[vals < stats::quantile(vals, probs = disc_prop),]

    # Increase the number of particles until having at least the initial number of particles
    
    num_to_add <- num_points - nrow(pars) 
    indices <- sample(x = seq_len(nrow(pars)), size = num_to_add, replace = TRUE)
    pars_to_add <- pars[indices, ] + data.frame(par1 = stats::rnorm(num_to_add, mean = 0, sd = sd_vec[1]),
                                                par2 = stats::rnorm(num_to_add, mean = 0, sd = sd_vec[2]),
                                                par3 = stats::rnorm(num_to_add, mean = 0, sd = sd_vec[3]),
                                                par4 = stats::rnorm(num_to_add, mean = 0, sd = sd_vec[4]))
    pars <- rbind(pars, pars_to_add)

    # Decrease variation
    sd_vec = sd_vec - alpha
  }
  total_time = proc.time() - init_time

  out <- list("parameters" = pv,
              "time" = total_time,
              "fhatdiff" = fhatdiff)
  return(out)
}
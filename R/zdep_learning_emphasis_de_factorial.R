
#' @keywords internal
get_results_factorial <- function(pars, input, num_threads, num_points) {
  dmval <- rcpp_mce_grid_factorial(as.matrix(pars),
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
  
  local_max_miss <- input$max_missing
  local_max_lambda <- input$max_lambda
  local_max_N <- input$maxN
  while (sum(is.na(dmval$results[, 1])) == num_points) {
    local_max_miss <- local_max_miss * 10
    local_max_lambda <- local_max_lambda * 10
    local_max_N <- local_max_N * 10
    if (local_max_N > 10000) {
      stop("exceeded maxN")
    }
    dmval <- rcpp_mce_grid_factorial(as.matrix(pars),
                                     brts = input$brts,
                                     sample_size = input$sample_size,
                                     maxN = local_max_N,
                                     soc = 2,
                                     max_missing = local_max_miss,
                                     max_lambda = local_max_lambda,
                                     lower_bound = input$lower_bound,
                                     upper_bound = input$upper_bound,
                                     xtol_rel = 0.1,
                                     num_threads = num_threads)
  }
  return(dmval)
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
emphasis_de_factorial <- function(brts,
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
  
  alpha <- sd_vec / num_iterations
  
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
  logf_grid <- list()
  logg_grid <- list()
  rejl_count <- rejo_count <- NULL
  
  min_loglik <- c()
  min_pars <- c()
  mean_loglik <- c()
  mean_pars <- c()
  
  
  if (verbose) {
    pb <- progress::progress_bar$new(format = "Progress: [:bar] :percent",
                                     total = num_iterations)
  }
  fhatdiff <- c()
  for (k in 1:num_iterations) {
    
    dmval <- get_results_factorial(pars, input, num_threads, num_points)
    logf_grid[[k]] <- dmval$logf
    logg_grid[[k]] <- dmval$logg
    dmval <- dmval$results # this is the old results matrix.
    
    
    to_store <- cbind(pars, dmval[, 5], dmval[, 6], dmval[, 1])
    colnames(to_store) <- c("par1", "par2", "par3", "par4",
                            "logf", "logg", "fhat")
    pv[[k]] <- to_store
    
    # dmval is a matrix with columns:
    # 1 = fhat
    # 2 = rejected_lambda
    # 3 = rejected_overruns
    # 4 = rejected_zero_weights
    
    vals <- dmval[, 1] #as.numeric(dmval[seq(from = 1, to = length(dmval), by = 4)])
    fails <- which(is.na(vals)) # failed runs return fhat = -1
    rejl <- dmval[, 2] # as.numeric(dmval[seq(from = 2, to = length(dmval), by = 4)])
    rejo <- dmval[, 3] # as.numeric(dmval[seq(from = 3, to = length(dmval), by = 4)])
    
    # if there were any stopped simulations due to max_lambda constrain
    if (sum(rejl[fails]) > 0) {
      input$max_lambda <- input$max_lambda + 100
      rejl_count <- c(rejl_count, k)
    }
    
    if (sum(rejo[fails]) > 0) {
      input$max_missing <- input$max_missing + 1000
      rejo_count <- c(rejo_count, k)
    }
    
    # still, we remove the NA values
    wi <- which(!is.na(vals))
    pars <- pars[wi, ]
    vals <- -vals[wi]
    # Methods
    
    ## DE1
    minval <- which.min(vals)
    min_loglik <- c(min_loglik, vals[minval])
    min_pars <- rbind(min_pars, as.numeric(pars[minval, ]))
    
    ## DE1b
    mean_loglik <- c(mean_loglik, mean(vals))
    mean_pars <- rbind(mean_pars, colMeans(pars[1:4]))
    
    ## Saving variation in the estimations
    # sd_pars <- apply(pars, 2, stats::sd)
    fhatdiff <- c(fhatdiff, stats::sd(vals))
    
    pars <- update_pars(pars, num_points, disc_prop, vals,
                        lower_bound, upper_bound, sd_vec)
    
    # Decrease variation
    sd_vec <- sd_vec - alpha
    
    if (verbose) {
      #   cat("iteration: ", k, " par estim: ")
      #    last_min_pars <- min_pars[nrow(min_pars), ]
      #    cat(last_min_pars)
      #    cat(" sd: ", sd_pars, "\n")
      pb$tick()
    }
  }
  total_time <- proc.time() - init_time
  
  obtained_estim <- colMeans(min_pars)
  
  out <- list("parameters" = pv,
              "time" = total_time,
              "fhatdiff" = fhatdiff,
              "minloglik" = min_loglik,
              "meanloglik" = mean_loglik,
              "min_pars" = min_pars,
              "mean_pars" = mean_pars,
              "obtained_estim" = obtained_estim,
              "logf_grid" = logf_grid,
              "logg_grid" = logg_grid)
  return(out)
}

#' Parameter Discovery using Monte Carlo EM Algorithm
#'
#' This function implements a Monte Carlo Expectation Maximization (MCEM) algorithm
#' to find suitable parameters for a model given branching times (brts) and parameter bounds.
#'
#' @param brts Numeric vector of branching times.
#' @param lower_bound Numeric vector specifying the lower bounds for each parameter.
#' @param upper_bound Numeric vector specifying the upper bounds for each parameter.
#' @param verbose Logical; if TRUE, detailed progress and diagnostics are printed.
#' @return A list containing the final parameter estimates (`pars`) and a dataframe (`MCEM`)
#'         with iteration history including parameter estimates and their corresponding fhat values.
#' @examples
#' # Example usage:
#' brts_example <- runif(10, 0, 1)  # Example branching times
#' lower_bound_example <- c(0, 0, 0, 0)
#' upper_bound_example <- c(1, 1, 1, 1)
#' result <- pd_ML(brts_example, lower_bound_example, upper_bound_example, verbose = TRUE)
#'
#' @export
pd_ML <- function(brts, lower_bound, upper_bound, verbose = TRUE) {
  
  #brts <- ape::branching.times(phylo)
  num_points <- 1000

  pars <- matrix(nrow = num_points, ncol = length(lower_bound))
  for (i in seq_along(lower_bound)) {
    pars[, i] <- runif(num_points, min = lower_bound[i], max = upper_bound[i])
  }
  
  cat("Searching suitable initial parameters...\n")
  
  dmval <- rcpp_mce_grid(pars, brts = brts, 
                                   sample_size = 1, 
                                   maxN = 1, 
                                   soc = 2, 
                                   max_missing = 1e+4, 
                                   max_lambda = 1e+4, 
                                   lower_bound = lower_bound, 
                                   upper_bound = upper_bound, 
                                   xtol_rel = 0.01, 
                                   num_threads = 4)
  
  vals <- dmval[, 1]
  wi <- which(!is.na(vals))
  
  if (length(wi) == 0) {
    cat("Cannot find initial parameters for the model\n")
    return(NULL)
  }
  
  pars <- pars[wi, ]
  vals <- -vals[wi]
  minval <- which.min(vals)
  init_pars <- pars[minval, ]
  
  cat("Initial parameters:", init_pars, "\n")
  
  if (verbose) cat("Initial parameter search completed. Initializing MCEM algorithm...\n")
  
  mcem <- data.frame(iteration = integer(), par1 = numeric(), par2 = numeric(), 
                     par3 = numeric(), par4 = numeric(), fhat = numeric())  
  improvement <- TRUE
  j <- 1 
  max_iterations <- 1000 
  
  # Main loop
  while (improvement && j <= max_iterations) {
    results <- em_cpp(brts = brts, pars = init_pars, sample_size = 2, soc = 2, maxN = 1e+7, 
                      xtol_rel = 1e-4, max_missing = 1e+6, max_lambda = 1e+6, 
                      lower_bound = lower_bound, upper_bound = upper_bound, num_threads = 4,copy_trees = F)
    
    init_pars <- results$estimates
    mcem <- rbind(mcem, data.frame(iteration = j, par1 = init_pars[1], par2 = init_pars[2], 
                                   par3 = init_pars[3], par4 = init_pars[4], fhat = results$fhat))
    
    if (verbose) {
      cat(sprintf("Epoch: %d, Parameters: %s, fhat: %f\n", j, toString(init_pars), results$fhat))
    }
    
    if (!is.null(previous_mean_params)) {
      if (check_convergence_in_parameters(init_pars, previous_mean_params)) {
        cat("Convergence achieved based on parameter stability.\n")
        improvement <- FALSE
      }
    }
    
    if (j > min_iterations) {  # Start checking for convergence after the minimum number of iterations
      if (j > burn_in) {  # Start calculating means after the burn-in period
        recent_pars <- mcem[(j-burn_in):j, 2:5]  # Adjust indices as per your actual dataframe structure
        recent_fhat <- mcem$fhat[(j-burn_in):j]
        mean_pars <- colMeans(recent_pars)
        mean_fhat <- mean(recent_fhat)
        
        # You can now use `mean_pars` and `mean_fhat` for convergence checking
        # Example: Check if the mean of parameters has changed less than a tolerance
        if (check_convergence_in_parameters(mean_pars, previous_mean_pars, tolerance)) {
          if (verbose) cat("Convergence achieved based on mean parameter stability.\n")
          improvement <- FALSE
        }
        previous_mean_pars <- mean_pars
      }
    }
  
  j <- j + 1  # Increment iteration counter
  }
  
  return(list(pars = init_pars, MCEM = mcem))
}


pd_ML <- function(brts, lower_bound, upper_bound, verbose = TRUE, burn_in = 50, min_iterations = 100, tolerance = 0.001) {
  
  num_points <- 1000
  pars <- matrix(nrow = num_points, ncol = length(lower_bound))
  for (i in seq_along(lower_bound)) {
    pars[, i] <- runif(num_points, min = lower_bound[i], max = upper_bound[i])
  }
  
  cat("Searching suitable initial parameters...\n")
  
  dmval <- rcpp_mce_grid(pars, brts = brts, sample_size = 1, maxN = 1, soc = 2, max_missing = 1e+5, max_lambda = 1e+4, lower_bound = lower_bound, upper_bound = upper_bound, xtol_rel = 0.01, num_threads = 4)
  
  vals <- dmval[, 1]
  wi <- which(!is.na(vals))
  
  if (length(wi) == 0) {
    cat("Cannot find initial parameters for the model\n")
    return(NULL)
  }
  
  pars <- pars[wi, ]
  vals <- -vals[wi]
  minval <- which.min(vals)
  init_pars <- pars[minval, ]
  
  cat("Initial parameters:", init_pars, "\n")
  
  if (verbose) cat("Initial parameter search completed. Initializing MCEM algorithm...\n")
  
  mcem <- data.frame(iteration = integer(), par1 = numeric(), par2 = numeric(), par3 = numeric(), par4 = numeric(), fhat = numeric())  
  improvement <- TRUE
  j <- 1 
  max_iterations <- 1000 
  previous_mean_pars <- rep(NA, length(lower_bound))
  
  # Main loop
  while (improvement && j <= max_iterations) {
    results <- em_cpp(brts = brts, pars = init_pars, sample_size = 2, soc = 2, maxN = 1e+7, xtol_rel = 1e-4, max_missing = 1e+6, max_lambda = 1e+6, lower_bound = lower_bound, upper_bound = upper_bound, num_threads = 4, copy_trees = FALSE)
    
    init_pars <- results$estimates
    mcem <- rbind(mcem, data.frame(iteration = j, par1 = init_pars[1], par2 = init_pars[2], par3 = init_pars[3], par4 = init_pars[4], fhat = results$fhat))
    
    if (verbose) {
      cat(sprintf("Iteration: %d, Parameters: %s, fhat: %f\n", j, toString(init_pars), results$fhat))
    }
    
    if (j >= min_iterations && j > burn_in) {
      recent_pars <- mcem[(j-burn_in+1):j, 2:5]  # Adjust indices as per your actual dataframe structure
      recent_fhat <- mcem$fhat[(j-burn_in+1):j]
      mean_pars <- colMeans(recent_pars, na.rm = TRUE)
      mean_fhat <- mean(recent_fhat, na.rm = TRUE)
      
      if (verbose) cat(sprintf("Mean Parameters: %s, Mean fhat: %f\n", toString(mean_pars), mean_fhat))
      
      # Check convergence based on mean parameter stability
      if (all(abs(mean_pars - previous_mean_pars) < tolerance)) {
        if (verbose) cat("Convergence achieved based on mean parameter stability.\n")
        improvement <- FALSE
      }
      
      previous_mean_pars <- mean_pars
    }
    
    j <- j + 1  # Increment iteration counter
  }
  
  return(list(pars = init_pars, MCEM = mcem))
}



#' Check Convergence in Parameters
#'
#' Determines if the parameter estimates have stabilized (converged) based on a specified tolerance.
#'
#' @param current_params Numeric vector of current parameter estimates.
#' @param previous_params Numeric vector of previous parameter estimates.
#' @param tolerance Numeric value specifying the convergence tolerance.
#' @return Logical; TRUE if the change in all parameters is less than or equal to the tolerance, FALSE otherwise.
#' @examples
#' current_params <- c(0.5, 0.5, 0.5, 0.5)
#' previous_params <- c(0.49, 0.51, 0.49, 0.51)
#' tolerance <- 0.05
#' has_converged <- check_convergence_in_parameters(current_params, previous_params, tolerance)
#' print(has_converged)  # Should return TRUE if the change in parameters is within the tolerance
#' 
#' @export
check_convergence <- function(current_means, previous_means, tolerance = 1e-4) {
  if (length(current_means) != length(previous_means)) {
    stop("Current means and previous means must have the same length.")
  }
  
  # Calculate the absolute difference between current and previous means
  differences <- abs(current_means - previous_means)
  
  # Check if all differences are within the specified tolerance
  all_within_tolerance <- all(differences <= tolerance)
  
  return(all_within_tolerance)
}





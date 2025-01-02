#' function to perform one step of the E-M algorithm
#' @param brts vector of branching times
#' @param init_pars vector of initial parameter files
#' @param sample_size number of samples
#' @param maxN maximum number of failed trees
#' @param soc number of lineages at the root/crown (1/2)
#' @param max_missing maximum number of species missing
#' @param max_lambda maximum speciation rate
#' @param lower_bound vector of lower bound values for optimization, should
#' be equal in length to the vector of init_pars
#' @param upper_bound vector of upper bound values for optimization, should
#' be equal in length to the vector of init_pars
#' @param xtol_rel relative tolerance for optimization
#' @param num_threads number of threads used.
#' @return a list with the following components: 
#' \itemize{
#'  \item{trees}{list of trees}
#'  \item{rejected}{number of rejected trees}
#'  \item{rejected_overruns}{number of trees rejected due to too large size}
#'  \item{rejected_lambda}{number of trees rejected due to lambda errors}
#'  \item{rejected_zero_weights}{number of trees rejected to zero weight}
#'  \item{time_elapsed}{time used}
#'  \item{weights}{vector of weights}
#'  \item{fhat}{vector of fhat values}
#'  \item{logf}{vector of logf values}
#'  \item{logg}{vector of logg values}
#' }
#' @export
e_cpp <- function(brts, init_pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads=1) {
  .Call('_emphasis_rcpp_mce', PACKAGE = 'emphasis', brts, init_pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads)
}

em_cpp <- function(brts, pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads, copy_trees, rconditional = NULL) {
  .Call('_emphasis_rcpp_mcem', PACKAGE = 'emphasis', brts, pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads, copy_trees, rconditional)
}

augmentPD_old <- function(phylo, pars, maxN, max_missing, lower_bound, upper_bound) {
    num_threads=1
    sample_size = 1
    init_pars = pars
    soc = 2 
    max_lambda = max_missing*100
    xtol_rel = 0.00001
    brts = ape::branching.times(phylo)
    # brts = sort(max(brts) - brts)
    # brts = brts[-1]
  .Call('_emphasis_rcpp_mce', PACKAGE = 'emphasis', brts, init_pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads)
}


#' Monte Carlo Grid Calculation
#'
#' This function performs a Monte Carlo Expectation-Maximization (MCEM) grid calculation
#' for given parameters, using Rcpp for improved performance. It is designed to work
#' with the 'emphasis' package.
#'
#' @param pars_R Numeric matrix of parameters for the MCEM algorithm.
#' @param brts Numeric vector of branching times.
#' @param sample_size Integer, the number of samples to draw in the Monte Carlo simulation.
#' @param maxN Integer, the maximum number of iterations for the EM algorithm.
#' @param soc Numeric, scale parameter for optimization control. Defaults to 2.
#' @param max_missing Numeric, the maximum allowed missing information.
#' @param max_lambda Numeric, the maximum allowed lambda value.
#' @param lower_bound Numeric vector, the lower bounds for each parameter.
#' @param upper_bound Numeric vector, the upper bounds for each parameter.
#' @param xtol_rel Numeric, the relative tolerance for optimization convergence. Defaults to 0.00001.
#' @param num_threads Integer, the number of threads to use for parallel computation.
#'
#' @return Returns the result of the MCEM grid calculation as a numeric matrix.
#' @export
#' @examples
#' # Define parameters for the MCEM grid calculation
#' pars_R <- matrix(runif(20), nrow = 5) # Example parameter matrix
#' brts <- runif(10, 0, 1) # Example branching times
#' sample_size <- 100
#' maxN <- 1000
#' lower_bound <- rep(0, ncol(pars_R))
#' upper_bound <- rep(1, ncol(pars_R))
#' result <- mcemGridCalculation(pars_R, brts, sample_size, maxN, 2, 1e4, 1e4, lower_bound, upper_bound, 0.00001, 1)
#' print(result)
mcGrid <- function(pars_R, brts, sample_size, maxN, soc=2, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel=0.00001, num_threads) {
  .Call('_emphasis_rcpp_mce_grid', PACKAGE = 'emphasis', pars_R, brts, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads)
}

rcpp_mce_grid <- function(pars_R, brts, sample_size, maxN, soc=2, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel=0.00001, num_threads) {
  .Call('_emphasis_rcpp_mce_grid', PACKAGE = 'emphasis', pars_R, brts, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads)
}

rcpp_mce_grid_factorial <- function(pars_R, brts, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads) {
  .Call('_emphasis_rcpp_mce_grid_factorial', PACKAGE = 'emphasis', pars_R, brts, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads)
}

#' function to perform one step of the E-M algorithm
#' @param brts vector of branching times
#' @param init_pars vector of initial parameter files
#' @param sample_size number of samples
#' @param maxN maximum number of failed trees
#' @param soc number of lineages at the root/crown (1/2)
#' @param max_missing maximum number of species missing
#' @param max_lambda maximum speciation rate
#' @param lower_bound vector of lower bound values for optimization, should
#' be equal in length to the vector of init_pars
#' @param upper_bound vector of upper bound values for optimization, should
#' be equal in length to the vector of init_pars
#' @param xtol_rel relative tolerance for optimization
#' @param num_threads number of threads used.
#' @param copy_trees if set to true, the trees generated are returned as well
#' @param rconditional R function that evaluates the GAM function.
#' @return a list with the following components: 
#' \itemize{
#'  \item{trees}{list of trees}
#'  \item{rejected}{number of rejected trees}
#'  \item{rejected_overruns}{number of trees rejected due to too large size}
#'  \item{rejected_lambda}{number of trees rejected due to lambda errors}
#'  \item{rejected_zero_weights}{number of trees rejected to zero weight}
#'  \item{estimates}{vector of estimates}
#'  \item{nlopt}{nlopt status}
#'  \item{fhat}{vector of fhat values}
#'  \item{time}{time elapsed}
#'  \item{weights}{vector of weights}
#' }
#' @export
em_cpp <- function(brts, pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads, copy_trees, rconditional = NULL) {
  .Call('_emphasis_rcpp_mcem', PACKAGE = 'emphasis', brts, pars, sample_size, maxN, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads, copy_trees, rconditional)
}

#' function to perform one step of the E-M algorithm
#' @param e_step result of e_step function, as a list
#' @param init_pars vector of initial parameter values
#' @param plugin string indicating plugin used, currently available: 'rpd1' and
#' 'rpd5c'
#' @param lower_bound vector of lower bound values for optimization, should
#' be equal in length to the vector of init_pars
#' @param upper_bound vector of upper bound values for optimization, should
#' be equal in length to the vector of init_pars 
#' @param xtol_rel relative tolerance for optimization
#' @param num_threads number of threads used.
#' @param rconditional R function that evaluates the GAM function.
#' @return list with the following entries:
#' \itemize{
#'  \item{estimates}{vector of estimates}
#'  \item{nlopt}{nlopt status}
#'  \item{time}{used computation time} 
#' }
#' @export
m_cpp <- function(e_step, init_pars, plugin, lower_bound, upper_bound, xtol_rel, num_threads, rconditional = NULL) {
  .Call('_emphasis_rcpp_mcm', PACKAGE = 'emphasis', e_step, init_pars, plugin, lower_bound, upper_bound, xtol_rel, num_threads, rconditional)
}





#' @keywords internal
get_results <- function(pars, input, num_threads=1) {
  dmval <- emphasis:::rcpp_mce_grid(as.matrix(pars),
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
  
  # dmval is a matrix with columns:
  # 1 = fhat
  # 2 = rejected_lambda
  # 3 = rejected_overruns
  # 4 = rejected_zero_weights
  num_points = nrow(dmval)
  rejl <- dmval[, 2] # as.numeric(dmval[seq(from = 2, to = length(dmval), by = 4)])
  rejo <- dmval[, 3] # as.numeric(dmval[seq(from = 3, to = length(dmval), by = 4)])
  rejw <- dmval[, 4] #as.numeric(dmval[seq(from = 4, to = length(dmval), by = 4)])
  rejected = cbind(rejl,rejo,rejw)
  colnames(rejected) = c("rate","overruns","weights")
  
  statistics <- cbind(dmval[, 5], dmval[, 6], dmval[, 1])
  colnames(statistics) <- c("logf", "logg", "fhat")
  
  empty = (sum(is.na(statistics[,"fhat"])) == num_points)
  
  return(list(statistics=statistics,rejected=rejected, empty = empty))
}


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
#' 
learn_em <- function(phylo,
                     max_iterations=100,
                     max_missing = 1e+5,
                     lower_bound = c(0,0,-0.1,-0.1),
                     upper_bound = c(1,5,0.1,0.1),
                     num_threads = 1){
  verbose = TRUE
  cat("Searching suitable initial parameters...\n")
  # get first grid
  
  num_points = 100
  pars <- get_random_grid(num_points=num_points,
                          lower_bound = lower_bound,
                          upper_bound = upper_bound)
  
  
  times <- Sys.time()
  max_lambda = 1e+5
  brts = ape::branching.times(phylo)
  input <- list(brts = brts,
                max_missing = max_missing,
                lower_bound = lower_bound,
                upper_bound =  upper_bound,
                sample_size = 1,
                maxN = 1,
                max_lambda = max_lambda)
  
  simulation <- get_results(pars, input, num_threads = 1)
  
  wi <- which(!is.na(simulation$statistics[,"fhat"]))
  
  if (length(wi) == 0) {
    cat("Cannot find initial parameters for the model\n")
    return(NULL)
  }
  
  pars_valid <- pars[wi, ]
  vals <- -simulation$statistics[wi,"fhat"]
  minval <- which.min(vals)
  init_pars <- as.numeric(pars_valid[minval, ])
  
  cat("Initial parameters:", init_pars, "\n")
  
  if (verbose) cat("Initial parameter search completed. Initializing MCEM algorithm...\n")
  
  mcem <- data.frame(iteration = integer(), par1 = numeric(), par2 = numeric(), 
                     par3 = numeric(), par4 = numeric(), loglik = numeric())  
  improvement <- TRUE
  j <- 1 
  max_iterations <- 100
  
  # Main loop
  while (j <= max_iterations) {
    results <- em_cpp(brts = brts, pars = init_pars, sample_size = 10, soc = 2, maxN = 1e+5, 
                      xtol_rel = 1e-4, max_missing = max_missing, max_lambda = max_lambda, 
                      lower_bound = lower_bound, upper_bound = upper_bound, num_threads = 2,copy_trees = T)
    
    e_cpp
    
    init_pars <- results$estimates
    mcem <- rbind(mcem, data.frame(iteration = j, par1 = init_pars[1], par2 = init_pars[2], 
                                   par3 = init_pars[3], par4 = init_pars[4], fhat = results$fhat))
    
    if (verbose) {
      cat(sprintf("Epoch: %d, Parameters: %s, fhat: %f\n", j, toString(init_pars), results$fhat))
    }
    
    burn_in = 30
    if (j > burn_in) {  # Start checking for convergence after the minimum number of iterations
      recent_pars <- mcem[(j-burn_in):j, 2:5]  # Adjust indices as per your actual dataframe structure
      recent_fhat <- mcem$fhat[(j-burn_in):j]
      mean_pars <- colMeans(recent_pars)
      mean_fhat <- mean(recent_fhat)
      
      
      previous_mean_pars <- mean_pars
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
    results <- em_cpp(brts = brts, pars = init_pars, sample_size = 2, soc = 2, maxN = 1e+6, xtol_rel = 1e-4, max_missing = 1e+6, max_lambda = 1e+6, lower_bound = lower_bound, upper_bound = upper_bound, num_threads = 4, copy_trees = FALSE)
    
    init_pars <- results$estimates
    mcem <- rbind(mcem, data.frame(iteration = j, par1 = init_pars[1], par2 = init_pars[2], par3 = init_pars[3], par4 = init_pars[4], fhat = results$fhat))
    
    if (verbose) {
      cat(sprintf("Iteration: %d, Parameters: %s, fhat: %f\n", j, toString(init_pars), results$fhat))
    }
    
    if (j >= 30 && j > burn_in) {
      recent_pars <- mcem[(j-burn_in+1):j, 2:5]  # Adjust indices as per your actual dataframe structure
      recent_fhat <- mcem$fhat[(j-burn_in+1):j]
      mean_pars <- colMeans(recent_pars, na.rm = TRUE)
      mean_fhat <- mean(recent_fhat, na.rm = TRUE)
      
      if (verbose) cat(sprintf("Mean Parameters: %s, Mean fhat: %f\n", toString(mean_pars), mean_fhat))
      
      # Check convergence based on mean parameter stability
      
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


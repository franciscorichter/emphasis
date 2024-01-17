#' Generate Variates for a Non-Homogeneous Exponential Process
#'
#' This function generates random variates for a non-homogeneous exponential
#' distribution using a thinning algorithm.
#'
#' @param num_variates The number of variates to generate.
#' @param covariates A numeric matrix of covariates where each row represents a 
#'                   set of covariates and the last column is expected to be 
#'                   time-dependent.
#' @param parameters Parameters for the exponential rate function, including the bias.
#' @param start_time The start time for the generation of variates.
#' @param max_time The maximum time limit for the generation of variates.
#'
#' @return A numeric vector of generated variates.
#'
#' print(generated_variates)
#'
#' @export
generateNonHomogeneousExp <- function(num_variates, rate_func, start_time, max_time) {
  if (max_time <= start_time) {
    stop("'max_time' must be greater than 'start_time'.")
  }
  
  variates <- numeric(0)  # Initialize an empty vector for the variates
  
  while (length(variates) < num_variates) {
    # Find the maximum rate over the interval
    time_points <- seq(start_time, max_time, length.out = 100)
    rate_values <- sapply(time_points, function(t) rate_func(t))
    max_rate <- max(rate_values)
    
    # Generate a candidate event time
    candidate_time <- start_time + rexp(1, max_rate)
    print(length(variates))
    if (candidate_time < max_time) {
      # Evaluate the true rate at the candidate time
      true_rate <- rate_func(candidate_time)
      
      # Accept or reject the candidate time
      if (runif(1) < true_rate / max_rate) {
        variates <- c(variates, candidate_time)
      }
    }
  }
  
  return(variates)
}


nhExpRand <- function(n, rate_func, now = 0, tMax = Inf) {
  if (!is.function(rate_func)) {
    stop("rate_func must be a function.")
  }
  if (tMax == Inf) {
    stop("Need a valid tMax for computation")
  }
  if (now < 0) {
    stop("now must be greater than or equal to 0")
  }
  
  vars <- numeric(n)  # Initialize the output vector
  
  for (i in 1:n) {
    p <- runif(1)  # Generate a uniform random number
    f <- Vectorize(function(t) {
      1 - p - exp(-integrate(Vectorize(function(x) rate_func(x)), lower = now, upper = t, 
                             subdivisions = 2000, stop.on.error = FALSE)$value)
    })
    vars[i] <- suppressWarnings(uniroot(f, c(now, tMax), extendInt = "yes", tol = 0.0001)$root)
  }
  
  return(vars)
}


rate_t <- function(t, params, cov_funcs, use_exponential = FALSE) {
  cov_values <- sapply(cov_funcs, function(f) f(t))
  linear_combination <- sum(params * c(1, cov_values))
  
  if (use_exponential) {
    return(exp(linear_combination))
  } else {
    return(linear_combination)
  }
}



#' Compute the Rate of a Non-Homogeneous Exponential Process
#'
#' This function calculates the rate of a non-homogeneous exponential process
#' given a set of covariates and corresponding parameters.
#'
#' @param covariates A numeric vector or matrix of covariates.
#' @param parameters A numeric vector of parameters corresponding to the covariates,
#'                   where the first element is the bias term.
#'
#' @return The computed rate value.
#'
#'
#' @export
ExponentialRate <- function(covariates, parameters) {
  if (ncol(covariates) + 1 != length(parameters)) {
    stop("The length of 'parameters' must be one more than the number of 'covariates'.")
  }
  
  # First element of parameters is the bias
  bias = parameters[1]
  
  # The rest of the parameters correspond to covariates
  rate_parameters = parameters[-1]
  
  # Calculate the rate
  rate = sum(sapply(X = covariates, FUN = function(covariates) exp(bias + sum(rate_parameters * covariates))))
  return(rate)
}

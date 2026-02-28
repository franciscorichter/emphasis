#' Generate Phylogenetic Diversity Data
#'
#' This function simulates phylogenetic trees based on given parameters.
#' @param n_trees Number of trees to simulate.
#' @param mu_interval Interval for mu parameter.
#' @param lambda_interval Interval for lambda parameter.
#' @param betaN_interval Interval for betaN parameter.
#' @param betaP_interval Interval for betaP parameter.
#' @param max_lin Maximum number of lineages.
#' @param max_tries Maximum number of tries for simulation.
#' @return A list containing the simulated trees and parameters.
#' @examples
#' generatePhyloPD(n_trees = 10,
#'                mu_interval = c(0.1, 0.5),
#'                lambda_interval = c(0.1, 0.5),
#'                betaN_interval = c(0.1, 0.5),
#'                betaP_interval = c(0.1, 0.5),
#'                max_lin = 1e6,
#'                max_tries = 1)
#' @export
generatePhyloPD <- function(n_trees,
                            mu_interval,
                            lambda_interval,
                            betaN_interval,
                            betaP_interval,
                            max_lin = 1e+6,
                            max_tries = 1){

  trees <- vector("list", n_trees)
  extrees <- vector("list", n_trees)
  Lmats <- vector("list", n_trees)
  brds_s <- vector("list", n_trees)

  name.param <- c("mu","lambda", "betaN","betaP")
  true.param <- vector(mode='list', length=4)
  names(true.param) <- name.param

  for(j in 1:n_trees){

    lambda_sample <- runif(1, lambda_interval[1], lambda_interval[2])
    mu_sample <- runif(1, mu_interval[1], mu_interval[2])
    betaN_sample <- runif(1, betaN_interval[1], betaN_interval[2])
    betaP_sample <- runif(1, betaP_interval[1], betaP_interval[2])
    sim.param <- c(mu_sample, lambda_sample, betaN_sample, betaP_sample)

    outputs <- tryCatch({
      sim_tree_pd_cpp(pars = sim.param,
                      max_t = 1,
                      max_lin = max_lin,
                      max_tries = max_tries,useDDD = TRUE)
    }, error = function(e) NULL)

    if (is.list(outputs) && max(outputs$brts) == 1) {
      trees[[j]] <- outputs[[1]]
      extrees[[j]] <- outputs[[2]]
      Lmats[[j]] <- outputs[[3]]
      brds_s[[j]] <- outputs[[4]]
      for (i in 1:4) {
        true.param[[i]] <- c(true.param[[i]], sim.param[i])
      }
    }

    message(sprintf("Progress: %d/%d", j, n_trees))
  }

  results = list(trees = trees, param = true.param, tas = extrees, L = Lmats, brts = brds_s)
  return(results)
}

#' Generate Variates for a Non-Homogeneous Exponential Process
#'
#' This function generates random variates for a non-homogeneous exponential
#' distribution using a thinning algorithm.
#'
#' @param num_variates The number of variates to generate.
#' @param rate_func A function representing the rate of the non-homogeneous process.
#' @param start_time The start time for the generation of variates.
#' @param max_time The maximum time limit for the generation of variates.
#' @return A numeric vector of generated variates.
#' @examples
#' rate_function_example <- function(t) { 2 * t }
#' variates <- generateNonHomogeneousExp(5, rate_function_example, 0, 10)
#' print(variates)
#' @export
generateNonHomogeneousExp <- function(num_variates, rate_func, start_time, max_time) {
  if (max_time <= start_time) {
    stop("'max_time' must be greater than 'start_time'.")
  }

  variates <- numeric(0)

  while (length(variates) < num_variates) {
    time_points <- seq(start_time, max_time, length.out = 100)
    rate_values <- sapply(time_points, function(t) rate_func(t))
    max_rate <- max(rate_values)

    candidate_time <- start_time + rexp(1, max_rate)

    if (candidate_time < max_time) {
      true_rate <- rate_func(candidate_time)

      if (runif(1) < true_rate / max_rate) {
        variates <- c(variates, candidate_time)
      }
    }
  }

  return(variates)
}

#' Generate Non-Homogeneous Exponential Random Variables
#'
#' This function generates random variables from a non-homogeneous exponential
#' distribution based on a provided rate function.
#'
#' @param n The number of variables to generate.
#' @param rate_func The rate function of the non-homogeneous process.
#' @param now The current time, default is 0.
#' @param tMax The maximum time, default is Inf.
#' @return A numeric vector of generated variables.
#' @examples
#' rate_function_example <- function(t) { 1 + sin(t) }
#' nh_exp_vars <- nhExpRand(10, rate_function_example, now = 0, tMax = 20)
#' print(nh_exp_vars)
#' @export
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

  vars <- numeric(n)

  for (i in 1:n) {
    p <- runif(1)
    f <- Vectorize(function(t) {
      1 - p - exp(-integrate(Vectorize(function(x) rate_func(x)), lower = now, upper = t,
                             subdivisions = 2000, stop.on.error = FALSE)$value)
    })
    vars[i] <- suppressWarnings(uniroot(f, c(now, tMax), extendInt = "yes", tol = 0.0001)$root)
  }

  return(vars)
}

#' Compute the Rate at a Given Time for a Non-Homogeneous Process
#'
#' This function calculates the rate at a specific time point for a non-homogeneous
#' process, based on a set of covariates and their corresponding parameters. The rate
#' can be computed either as a linear combination of the parameters and covariates or
#' as an exponential of this linear combination.
#'
#' @param t Numeric, the time at which the rate is to be calculated.
#' @param params Numeric vector, parameters for the rate function including the baseline rate.
#' @param cov_funcs List of functions, each representing a covariate as a function of time.
#' @param use_exponential Logical, if TRUE, the exponential of the linear combination
#'        is returned, otherwise, the linear combination itself is returned.
#'
#' @return Numeric, the calculated rate at time `t`.
#'
#' @examples
#' # Define example covariate functions
#' cov_func1 <- function(t) { sin(t) }
#' cov_func2 <- function(t) { cos(t) }
#'
#' # Example parameters (baseline and coefficients for covariates)
#' params <- c(0.5, 1.2, -0.8)
#'
#' # Compute the rate at a specific time
#' rate_at_time_5 <- rate_t(t = 5, params = params,
#'                          cov_funcs = list(cov_func1, cov_func2),
#'                          use_exponential = TRUE)
#' print(rate_at_time_5)
#'
#' @export
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
#' @return The computed rate value.
#' @examples
#' covariates_example <- matrix(c(1, 2, 3, 4), ncol = 2)
#' parameters_example <- c(0.5, 1, -0.5)
#' rate <- ExponentialRate(covariates_example, parameters_example)
#' print(rate)
#' @export
ExponentialRate <- function(covariates, parameters) {
  if (ncol(covariates) + 1 != length(parameters)) {
    stop("The length of 'parameters' must be one more than the number of 'covariates'.")
  }

  bias = parameters[1]
  rate_parameters = parameters[-1]

  rate = sum(sapply(X = covariates, FUN = function(covariates) exp(bias + sum(rate_parameters * covariates))))
  return(rate)
}

#' Generate Phylogenetic Diversity Data
#'
#' Simulates multiple phylogenetic trees under a PD-dependent diversification
#' model.  Parameters are drawn uniformly from user-supplied intervals.
#'
#' @param n_trees Number of trees to simulate.
#' @param mu_interval Length-2 numeric vector: lower and upper bounds for the
#'   extinction rate (\code{mu}).  Lower must be \eqn{\le} upper.
#' @param lambda_interval Interval for the baseline speciation rate.
#' @param betaN_interval Interval for the diversity-dependence coefficient.
#' @param betaP_interval Interval for the PD-dependence coefficient.
#' @param max_t Crown age for simulated trees.  Default \code{1}.
#' @param max_lin Maximum number of lineages.  Default \code{1e6}.
#' @param max_tries Maximum number of retries per tree.  Default \code{1}.
#' @param seed Optional integer seed for reproducibility.
#' @return A list with elements \code{trees}, \code{param}, \code{tas},
#'   \code{L}, \code{brts}.  Failed trees are filtered out; all vectors have
#'   the same (possibly shorter) length.
#' @examples
#' \dontrun{
#' generatePhyloPD(
#'   n_trees = 10,
#'   mu_interval = c(0.1, 0.5),
#'   lambda_interval = c(0.1, 0.5),
#'   betaN_interval = c(0.1, 0.5),
#'   betaP_interval = c(0.1, 0.5),
#'   max_t = 5
#' )
#' }
#' @export
generatePhyloPD <- function(n_trees,
                            mu_interval,
                            lambda_interval,
                            betaN_interval,
                            betaP_interval,
                            max_t = 1,
                            max_lin = 1e6,
                            max_tries = 1,
                            seed = NULL) {
  # --- Input validation --------------------------------------------------
  .validate_interval <- function(x, nm) {
    if (!is.numeric(x) || length(x) != 2L) {
      stop(sprintf("'%s' must be a length-2 numeric vector.", nm))
    }
    if (x[1] > x[2]) {
      stop(sprintf("'%s' lower bound (%.4g) > upper bound (%.4g).", nm, x[1], x[2]))
    }
  }
  .validate_interval(mu_interval, "mu_interval")
  .validate_interval(lambda_interval, "lambda_interval")
  .validate_interval(betaN_interval, "betaN_interval")
  .validate_interval(betaP_interval, "betaP_interval")
  if (!is.numeric(n_trees) || n_trees < 1L) {
    stop("'n_trees' must be a positive integer.")
  }

  if (!is.null(seed)) set.seed(seed)

  trees <- vector("list", n_trees)
  extrees <- vector("list", n_trees)
  Lmats <- vector("list", n_trees)
  brds_s <- vector("list", n_trees)

  name.param <- c("mu", "lambda", "betaN", "betaP")
  true.param <- vector(mode = "list", length = 4L)
  names(true.param) <- name.param

  pb <- progress::progress_bar$new(
    format = "  simulating [:bar] :current/:total (:percent) eta: :eta",
    total = n_trees, clear = FALSE, width = 60
  )

  for (j in seq_len(n_trees)) {
    lambda_sample <- stats::runif(1, lambda_interval[1], lambda_interval[2])
    mu_sample <- stats::runif(1, mu_interval[1], mu_interval[2])
    betaN_sample <- stats::runif(1, betaN_interval[1], betaN_interval[2])
    betaP_sample <- stats::runif(1, betaP_interval[1], betaP_interval[2])
    sim.param <- c(mu_sample, lambda_sample, betaN_sample, betaP_sample)

    outputs <- tryCatch(
      sim_tree_pd_cpp(
        pars = sim.param,
        max_t = max_t,
        max_lin = max_lin,
        max_tries = max_tries,
        useDDD = TRUE
      ),
      error = function(e) NULL
    )

    if (is.list(outputs) && !is.null(outputs$brts) && max(outputs$brts) == max_t) {
      trees[[j]] <- outputs[[1]]
      extrees[[j]] <- outputs[[2]]
      Lmats[[j]] <- outputs[[3]]
      brds_s[[j]] <- outputs[[4]]
      for (i in seq_len(4L)) {
        true.param[[i]] <- c(true.param[[i]], sim.param[i])
      }
    }

    pb$tick()
  }

  # Filter out NULL entries from failed simulations
  keep <- !vapply(trees, is.null, logical(1))
  trees <- trees[keep]
  extrees <- extrees[keep]
  Lmats <- Lmats[keep]
  brds_s <- brds_s[keep]

  if (sum(keep) < n_trees) {
    message(sprintf("%d/%d trees succeeded.", sum(keep), n_trees))
  }

  list(trees = trees, param = true.param, tas = extrees, L = Lmats, brts = brds_s)
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
#' @param grid_size Number of grid points for rate integration. Default \code{200L}.
#' @param max_iter Maximum thinning iterations before stopping. Default \code{1e6}.
#' @return A numeric vector of generated variates.
#' @examples
#' rate_function_example <- function(t) {
#'   2 * t
#' }
#' variates <- generateNonHomogeneousExp(5, rate_function_example, 0, 10)
#' print(variates)
#' @export
generateNonHomogeneousExp <- function(num_variates,
                                      rate_func,
                                      start_time,
                                      max_time,
                                      grid_size = 200L,
                                      max_iter = 1e6) {
  if (!is.function(rate_func)) {
    stop("'rate_func' must be a function.")
  }
  if (max_time <= start_time) {
    stop("'max_time' must be greater than 'start_time'.")
  }
  if (num_variates < 1L) {
    stop("'num_variates' must be >= 1.")
  }

  variates <- numeric(num_variates)
  count <- 0L
  iter <- 0L

  # Pre-compute max rate over the interval
  time_points <- seq(start_time, max_time, length.out = grid_size)
  rate_values <- vapply(time_points, rate_func, numeric(1))
  max_rate <- max(rate_values)

  if (max_rate <= 0) {
    stop("The rate function is non-positive over [start_time, max_time].")
  }

  while (count < num_variates) {
    iter <- iter + 1L
    if (iter > max_iter) {
      stop(sprintf(
        "generateNonHomogeneousExp: exceeded max_iter (%d) with only %d/%d variates generated.",
        as.integer(max_iter), count, num_variates
      ))
    }

    candidate_time <- start_time + stats::rexp(1, max_rate)

    if (candidate_time < max_time) {
      true_rate <- rate_func(candidate_time)
      if (stats::runif(1) < true_rate / max_rate) {
        count <- count + 1L
        variates[count] <- candidate_time
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
#' rate_function_example <- function(t) {
#'   1 + sin(t)
#' }
#' nh_exp_vars <- nhExpRand(10, rate_function_example, now = 0, tMax = 20)
#' print(nh_exp_vars)
#' @export
nhExpRand <- function(n, rate_func, now = 0, tMax = Inf) {
  if (!is.function(rate_func)) {
    stop("rate_func must be a function.")
  }
  if (!is.finite(tMax) || tMax <= now) {
    stop("Need a valid finite tMax > now for computation.")
  }
  if (now < 0) {
    stop("now must be greater than or equal to 0.")
  }
  if (n < 1L) {
    stop("'n' must be >= 1.")
  }

  vars <- numeric(n)

  for (i in seq_len(n)) {
    p <- stats::runif(1)
    f <- Vectorize(function(t) {
      1 - p - exp(-stats::integrate(
        Vectorize(function(x) rate_func(x)),
        lower = now, upper = t,
        subdivisions = 2000, stop.on.error = FALSE
      )$value)
    })
    result <- tryCatch(
      stats::uniroot(f, c(now, tMax), extendInt = "yes", tol = 1e-4),
      warning = function(w) {
        # Retry without extendInt on warnings (e.g. "values at end points not
        # of opposite sign")
        tryCatch(
          stats::uniroot(f, c(now, tMax * 1.5), tol = 1e-4),
          error = function(e2) NULL
        )
      },
      error = function(e) NULL
    )
    root <- if (!is.null(result)) result$root else NA_real_
    vars[i] <- if (!is.na(root) && root > tMax) tMax else root
  }

  n_na <- sum(is.na(vars))
  if (n_na > 0L) {
    warning(sprintf("nhExpRand: %d/%d variates could not be computed (returned NA).", n_na, n))
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
#' cov_func1 <- function(t) {
#'   sin(t)
#' }
#' cov_func2 <- function(t) {
#'   cos(t)
#' }
#'
#' # Example parameters (baseline and coefficients for covariates)
#' params <- c(0.5, 1.2, -0.8)
#'
#' # Compute the rate at a specific time
#' rate_at_time_5 <- rate_t(
#'   t = 5, params = params,
#'   cov_funcs = list(cov_func1, cov_func2),
#'   use_exponential = TRUE
#' )
#' print(rate_at_time_5)
#'
#' @export
rate_t <- function(t, params, cov_funcs, use_exponential = FALSE) {
  if (!is.list(cov_funcs)) {
    stop("'cov_funcs' must be a list of functions.")
  }
  if (length(params) != 1L + length(cov_funcs)) {
    stop(sprintf(
      "'params' must have length %d (1 baseline + %d covariates), got %d.",
      1L + length(cov_funcs), length(cov_funcs), length(params)
    ))
  }

  cov_values <- vapply(cov_funcs, function(f) f(t), numeric(1))
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
  if (!is.matrix(covariates)) {
    covariates <- as.matrix(covariates)
  }
  if (ncol(covariates) + 1L != length(parameters)) {
    stop("The length of 'parameters' must be one more than the number of covariate columns.")
  }

  bias <- parameters[1]
  rate_parameters <- parameters[-1]

  # For each row of covariates, compute exp(bias + beta . x)
  rate <- sum(apply(covariates, 1, function(row) {
    exp(bias + sum(rate_parameters * row))
  }))
  return(rate)
}

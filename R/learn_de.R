#' @keywords internal
get_results <- function(pars, input, num_threads) {
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



#' @keywords internal
get_random_grid <- function(num_points,
                            lower_bound,
                            upper_bound) {
  
  dim <- length(lower_bound)
  pars = NULL
  for(i in 1:dim){
    pa <- runif(num_points,
                  min = lower_bound[i],
                  max = upper_bound[i])
    pars = cbind(pars,pa)

  }
  pars = as.data.frame(pars)
  names(pars) <- paste0("par", 1:dim)
  
  return(pars)
}

#' @keywords internal
update_parameters <- function(pars, vals, Cr = 0.9, num_points) {
  n <- nrow(pars)
  d <- ncol(pars)
  
  # Default mutation factor
  mutation_factor <- 0.8
  
  # Convert objective values to weights inversely (smaller values should have larger weights)
  weights <- 1 / vals
  weights <- weights / sum(weights)
  
  new_pars <- pars
  
  for (i in 1:n) {
    # Select three candidates using weights
    indices <- sample(1:n, 3, prob = weights, replace = TRUE)
    r1 <- indices[1]
    r2 <- indices[2]
    r3 <- indices[3]
    
    # Mutation
    mutant_vector <- pars[r1, ] + mutation_factor * (pars[r2, ] - pars[r3, ])
    
    # Crossover
    trial_vector <- pars[i, ]
    j_rand <- sample(1:d, 1)
    for (j in 1:d) {
      if (runif(1) <= Cr || j == j_rand) {
        trial_vector[j] <- mutant_vector[j]
      }
    }
    
    # Selection using weights
    candidate_vectors <- rbind(trial_vector, pars[i, ])
    chosen_index <- sample(c(1, 2), 1, prob = c(weights[i], 1 - weights[i]))
    new_pars[i, ] <- candidate_vectors[chosen_index, ]
  }
  
  # Fill in the population to maintain the original size with slight perturbations
  current_size <- nrow(new_pars)
  if (current_size < num_points) {
    # Using the best parameters of the last generation
    additional_indices <- sample(1:n, num_points - current_size, prob = weights, replace = TRUE)
    additional_pars <- pars[additional_indices, , drop = FALSE]
    
    # Apply slight perturbations to the additional parameters
    perturbation_factor <- 0.001
    perturbations <- matrix(rnorm(nrow(additional_pars) * d, mean = 0, sd = perturbation_factor), nrow = nrow(additional_pars), ncol = d)
    additional_pars <- additional_pars + perturbations
    
    new_pars <- rbind(new_pars, additional_pars)
  }
  
  return(new_pars)
}



#' perform emphasis analysis using DE method
#' @param brts branching times of tree to fit on
#' @param max_iterations number of iterations of the DE algorithm
#' @param num_points number of particles per iteration
#' @param max_missing maximum number of missing trees
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
learn_de <- function(phylo,
                     max_iterations=100,
                     num_points = 1000,
                     max_missing = 100000,
                     lower_bound = c(0,0,-0.1,-0.1),
                     upper_bound = c(1,5,0.1,0.1),
                     verbose = TRUE,
                     num_threads = 1) {
  
  if (!is.numeric(lower_bound) || !is.numeric(upper_bound)) {
    stop("lower_bound and upper_bound must be numeric vectors.")
  }
  
  if (length(upper_bound) != length(lower_bound)) {
    stop("lower bound and upper bound vectors need to be same length")
  }
  
  
  
  init_time <- proc.time()
  

  

  # get first grid
  pars <- get_random_grid(num_points,
                          lower_bound = lower_bound,
                          upper_bound = upper_bound)
  
  # initialize
  pv <- list()

  min_loglik <- c()
  mean_loglik <- c()
  
  min_pars <- matrix(nrow=max_iterations, ncol=length(lower_bound))
  mean_pars <- matrix(nrow=max_iterations, ncol=length(lower_bound))
  
  fhatdiff <- c()
  times <- Sys.time()
  
  input <- list(brts = ape::branching.times(phylo),
                max_missing = max_missing,
                lower_bound = lower_bound,
                upper_bound =  upper_bound,
                sample_size = 1,
                maxN = 1,
                max_lambda = 10000)
  
  
  k = 0
  while (k < max_iterations) {
    k = k + 1 
    print("Simulating...")
    simulation <- get_results(pars, input, num_threads)
    
    if(simulation$empty){
      print("Not possible to simulate trees with this (meta)paramters")
    }

    pv[[k]] <- list(pars = pars, simulation = simulation)
    

    vals <- simulation$statistics[,"fhat"]

    wi <- which(!is.na(vals) & vals <0)
    pars <- pars[wi, ]
    vals <- -vals[wi]
    print(paste("Number of legit trees:",nrow(pars)))
   # plot(pars[,1],pars[,2])
   # pairs(pars)
    # Methods
    
    ## DE1
    minval <- which.min(vals)
    min_loglik <- c(min_loglik, vals[minval])
    min_pars[k,] <- as.numeric(pars[minval, ])
    
    ## DE1b
    mean_loglik <- c(mean_loglik, mean(vals))
    mean_pars[k,] <- colMeans(pars)
    
    pars <- update_parameters(pars, vals, num_points=num_points)
    
    if (verbose) {
      # Display iteration number
      cat(sprintf("Iteration: %d/%d\n", k, max_iterations))
      
      # Display summary statistics
      cat(sprintf("Minimum Log-Likelihood: %f\n", min_loglik[length(min_loglik)]))
      cat(sprintf("Mean Log-Likelihood: %f\n", mean_loglik[length(mean_loglik)]))
      
      # Display time taken for the current iteration
      cat(sprintf("Time for iteration: %f seconds\n", difftime(Sys.time(), times[length(times)], units="secs")))
      
      # Update progress bar
      
      flush.console()
    }
    
    
    times  <- c(times,Sys.time())
  }
  total_time <- proc.time() - init_time
  
  obtained_estim <- colMeans(min_pars)
  
  out <- list("parameters" = pv,
              "time" = total_time,
              "minloglik" = min_loglik,
              "meanloglik" = mean_loglik,
              "min_pars" = min_pars,
              "mean_pars" = mean_pars,
              "obtained_estim" = obtained_estim,
              "times" = times)
  return(out)
}



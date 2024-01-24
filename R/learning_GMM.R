nhExpRand



# GRADIENT DESCENT FOR EVOLUTIONARY MODEL WITH CONSTANT D MATRIX
gradDescentEvolutionModel <- function(N_max, stats_obs, n_stats, n_pars, learn_rate, patience,
                                      iters, initial_params, n_trees_D, n_trees_sgd, max_t, max_attempts,print_info=TRUE) {
  iteration_counter <- 0
  patience_counter <- 0
  pars_matrix <- matrix(nrow=0, ncol=length(initial_params))
  pars_matrix <- rbind(pars_matrix, t(initial_params))
  updated_params <- initial_params
  stats_diff_matrix <- matrix(nrow=0, ncol=n_stats)
  stats_diff_abs_vector <- c()
  execution_times <- c()
  
  previous_params <- updated_params
  
  # Initial scaling
  scaling <- abs(updated_params)
  D_matrix_info <- list()
  stats_diff_history <- NULL
  
  convergence_threshold <- 1e-3  
  
  while (iteration_counter < iters) {
    iteration_start_time <- Sys.time()
    
    # Update D matrix periodically
    
    if (iteration_counter == 0 | isStatsDiffIncreasing(stats_diff_history)) {
      if (print_info) cat("Calculating D matrix...")
      D_matrix <- calculateDMatrix(n_stats = n_stats, 
                                   n_pars = n_pars, 
                                   initial_params = initial_params,
                                   n_trees_D =  n_trees_D,
                                   max_attempts =  max_attempts,
                                   max_t =  max_t,
                                   use_parallel = FALSE)
      if (length(D_matrix) == 0) { D_matrix <- D_matrix_info[[length(D_matrix_info)]] }
      
      MD_matrix <- D_matrix / rowSums(abs(D_matrix))
      D_matrix_info <- c(D_matrix_info, list(D_matrix=D_matrix,iteration=iteration_counter))
    }
    
    generated_stats <- generateStatistics(N_max, updated_params, n_trees_sgd, max_attempts, max_t)
    
    # Check for errors in tree generation
    if (is.null(generated_stats)) {
      return(list(pars=pars_matrix, stats_diff=stats_diff_matrix, stats_diff_abs=stats_diff_abs_vector, D=D_matrix_info, times=execution_times))
    }
    
    # Calculate difference in statistics
    stats_diff <- colMedians(generated_stats) - stats_obs
    
    # Update parameters
    gradient_adjustment <- learn_rate * (MD_matrix %*% (as.matrix(stats_diff))) * scaling
    updated_params <- initial_params - gradient_adjustment
    
    # Update condition
    parameter_change <- updated_params - initial_params
    stats_diff_matrix <- rbind(stats_diff_matrix, stats_diff)
    stats_diff_abs_vector <- c(stats_diff_abs_vector, sum(abs(stats_diff)))
    stats_diff_history <- c(stats_diff_history, sum(abs(stats_diff)))
    
    # Update initial parameters for next iteration
    pars_matrix <- rbind(pars_matrix, t(initial_params))
    initial_params <- updated_params
    iteration_counter <- iteration_counter + 1
    execution_times <- c(execution_times, Sys.time() - iteration_start_time)
    
    previous_params <- updated_params
    
    # Logging
    if (print_info) {
      print(iteration_counter)
      cat("Updated Params =", initial_params, "\n")
      cat("Moving average of stats_diff_abs =", mean(tail(stats_diff_abs_vector, 20)), "\n")
    }
    
    if (sum(abs(updated_params - previous_params)) < convergence_threshold) {
      #  break  # Exit the loop
    }
    
  }
  pars_matrix <- rbind(pars_matrix, t(initial_params))
  total_runtime <- Sys.time() - iteration_start_time
  return(list(pars=pars_matrix, stats_diff=stats_diff_matrix, stats_diff_abs=stats_diff_abs_vector, D=D_matrix_info, times=execution_times))
}

isStatsDiffIncreasing <- function(history) {
  if (length(history) < 10) return(FALSE)
  return(mean(tail(history, 10)) > mean(tail(history, 9)) & 
           mean(tail(history, 9)) > mean(tail(history, 8))& 
           mean(tail(history, 8)) > mean(tail(history, 7))) 
}


generateStatistics <- function(N_max, params, n_trees_sgd, max_attempts, max_t) {
  stats_matrix <- matrix(nrow=0, ncol=n_stats)
  
  for (i in 1:n_trees_sgd) {
    result <- create_tree_mat_fixed_N(N_max, params[1:3], params[4:length(params)], 0, max_attempts, max_t)
    if (length(result) == 0 || result$attempt >= max_attempts) {
      return(NULL)
    }
    
    tree_mat <- result$tree_mat
    tree <- L2phylo(unname(tree_mat), dropextinct=TRUE)
    stats_matrix <- rbind(stats_matrix, stats_PDM(tree))
  }
  
  return(stats_matrix)
}



calc_stats_for_perturbed_params <- function(pars, N_max, n_stats, n_pars, n_trees_D, max_attempts, max_t) {
  s = matrix(nrow=0, ncol=n_stats)
  for (i in 1:n_trees_D) {
    result = create_tree_mat_fixed_N(N_max, pars[1:3], pars[4:n_pars], 0, max_attempts, max_t)
    if (length(result) == 0 || result$attempt >= max_attempts) {
      return(NULL)
    }
    tree_mat = result$tree_mat
    tree_extant = L2phylo(unname(tree_mat), dropextinct=TRUE)
    s = rbind(s, stats_PDM(tree_extant))
  }
  return(colMedians(s))
}

calc_D_stats <- function(n_stats, n_pars, pars_i, n_trees_D, max_attempts, max_t, print_info = TRUE) {
  eps <- pars_i / 10
  
  # Helper function to generate stats
  generate_stats <- function(pars) {
    s <- matrix(nrow = 0, ncol = n_stats)
    for (i in 1:n_trees_D) {
      result <- create_tree_mat_fixed_N(N_max, pars[1:3], pars[4:n_pars], 0, max_attempts, max_t)
      if (is.null(result) || result$attempt >= max_attempts) {
        return(NULL)
      }
      tree_mat <- result$tree_mat
      tree_extant <- L2phylo(unname(tree_mat), dropextinct = TRUE)
      s <- rbind(s, stats_PDM(tree_extant))
    }
    colMedians(s)
  }
  
  # Calculate original stats
  stats_orig <- generate_stats(pars_i)
  if (is.null(stats_orig)) {
    stop("Failed to generate original stats.")
  }
  
  # Initialize gradient matrix
  stats_grad <- matrix(nrow = n_pars, ncol = n_stats)
  
  # Calculate gradients
  for (p in 1:n_pars) {
    pars_eps <- pars_i
    pars_eps[p] <- pars_eps[p] + eps[p]
    stats_eps <- generate_stats(pars_eps)
    if (is.null(stats_eps)) {
      stop(paste("Failed to generate stats for parameter", p))
    }
    stats_grad[p, ] <- (stats_eps - stats_orig) / eps[p]
    if (print_info) {
      cat("Parameter", p, "done\n")
    }
  }
  
  stats_grad
}
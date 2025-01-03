# Ensure required packages are installed
#if (!requireNamespace("progress", quietly = TRUE)) {
#  install.packages("progress")
#}

emphasis_de_mod <- function(brts, num_iterations, num_points, max_missing, 
                        sd_vec, lower_bound, upper_bound, maxN = 10, 
                        max_lambda, disc_prop = 0.5, verbose = TRUE, 
                        num_threads = 1) {
  
  # Validate inputs
  validate_inputs(upper_bound, lower_bound, sd_vec)
  
  # Initialize
  init_time <- proc.time()
  alpha <- calculate_alpha(sd_vec, num_iterations)
  input <- initialize_input(brts, max_missing, lower_bound, upper_bound, maxN, max_lambda, disc_prop)
  
  # Get initial parameters
  pars <- emphasis:::get_random_grid(num_points, lower_bound = lower_bound, upper_bound = upper_bound)
  
  # Initialize storage vectors
  results_storage <- initialize_results_storage()
  
  # Main loop
  for (k in 1:num_iterations) {
    
    iter_res <- iteration_step(pars = pars,
                               input = input,
                               num_points =  num_points,
                               sd_vec = sd_vec,
                               alpha = alpha,
                               verbose = TRUE)
    # Update the results_storage with the results from iter_res
    #results_storage$pv <- rbind(results_storage$pv, iter_res$pv)
    results_storage$min_loglik <- c(results_storage$min_loglik, iter_res$min_loglik)
    results_storage$mean_loglik <- c(results_storage$mean_loglik, iter_res$mean_loglik)
    results_storage$min_pars <- rbind(results_storage$min_pars, iter_res$min_pars)
    results_storage$mean_pars <- rbind(results_storage$mean_pars, iter_res$mean_pars)
    results_storage$fhatdiff <- c(results_storage$fhatdiff, iter_res$fhatdiff)
    results_storage$dmval <- c(results_storage$dmval,iter_res$dmval)
    # Update pars and sd_vec for the next iteration
    pars <- iter_res$pars
    sd_vec <- iter_res$sd_vec
    
    if (verbose) {
      cat("Epoch:", k, "/", num_iterations, " - Likelihood:", iter_res$min_loglik, "\n")
    }
  }
  
  # Compile results
  out <- compile_results(results_storage, init_time)
  return(out)
}

# Additional helper functions

validate_inputs <- function(upper_bound, lower_bound, sd_vec) {
  if (length(upper_bound) != length(lower_bound)) {
    stop("lower bound and upper bound vectors need to be the same length")
  }
  
  if (length(upper_bound) != length(sd_vec)) {
    stop("sd vector is not the correct length")
  }
}

calculate_alpha <- function(sd_vec, num_iterations) {
  return(sd_vec / num_iterations)
}

initialize_input <- function(brts, max_missing, lower_bound, upper_bound, maxN, max_lambda, disc_prop) {
  return(list(brts = brts, max_missing = max_missing, lower_bound = lower_bound, 
              upper_bound = upper_bound, sample_size = 1, maxN = maxN, 
              max_lambda = max_lambda, disc_prop = disc_prop))
}

initialize_results_storage <- function() {
  return(list(
    pv = list(),
    rejl_count = NULL,
    rejo_count = NULL,
    min_loglik = c(),
    min_pars = c(),
    mean_loglik = c(),
    mean_pars = c(),
    fhatdiff = c(),
    times = Sys.time(),
    dmval = list()
  ))
}

iteration_step <- function(pars, 
                           input, 
                           num_points, 
                           sd_vec, 
                           alpha, 
                           verbose = FALSE,
                           num_threads=1) {
  # 1. Get Results
  dmval <- emphasis:::get_results(pars, input, num_threads, num_points)
  
  # 2. Store results
  to_store <- cbind(pars, dmval[, 5], dmval[, 6], dmval[, 1])
  colnames(to_store) <- c("par1", "par2", "par3", "par4", "logf", "logg", "fhat")
  
  # 3. Handling failures
  vals <- dmval[, 1]
  fails <- which(is.na(vals))
  rejl <- dmval[, 2]
  rejo <- dmval[, 3]
  
  # 4. Update pars
  wi <- which(!is.na(vals))
  pars <- pars[wi, ]
  vals <- -vals[wi]
  
  minval <- which.min(vals)
  min_loglik <- vals[minval]
  min_pars <- as.numeric(pars[minval, ])
  
  mean_loglik <- mean(vals)
  mean_pars <- colMeans(pars[1:4])
  
  fhatdiff <- stats::sd(vals)
  
  pars <- emphasis:::update_pars(pars, num_points, disc_prop, vals, lower_bound, upper_bound, sd_vec)
  
  # Decrease variation
  sd_vec <- sd_vec - alpha
  
  
  return(list(pv = to_store, 
              min_loglik = min_loglik, 
              mean_loglik = mean_loglik, 
              min_pars = min_pars, 
              mean_pars = mean_pars, 
              fhatdiff = fhatdiff, 
              pars = pars, 
              sd_vec = sd_vec,
              dmval=dmval))
}


compile_results <- function(results_storage, init_time) {
  obtained_estim <- colMeans(results_storage$min_pars)
  
  out <- list(
    "parameters" = results_storage$pv,
    "time" = proc.time() - init_time,
    "fhatdiff" = results_storage$fhatdiff,
    "minloglik" = results_storage$min_loglik,
    "meanloglik" = results_storage$mean_loglik,
    "min_pars" = results_storage$min_pars,
    "mean_pars" = results_storage$mean_pars,
    "obtained_estim" = obtained_estim,
    "times" = results_storage$times,
    "dmval" = results_storage$dmval
  )
  return(out)
}


process_and_plot <- function(file_path) {
  # Load necessary libraries
  library(ggplot2)
  
  # Load the data
  cat("Loading data...")
  load(file_path)

  computation_time <- time1 - time0
  cat("Data loaded")
  # Convert to hours, minutes, and seconds
  computation_hours <- computation_time['elapsed'] %/% 3600
  computation_minutes <- (computation_time['elapsed'] %% 3600) %/% 60
  computation_seconds <- computation_time['elapsed'] %% 60

  print(paste("Simulation took ",computation_hours, "hours,", computation_minutes, "minutes and", round(computation_seconds), "seconds"))


  # Let's assume 'parameters' is a named vector or list of your input parameters
  parameters <- c(betaN_interval = betaN_interval, betaP_interval =  betaP_interval, 
                lambda_interval = lambda_interval, max_lin = max_lin, max_tries = max_tries, 
                mu_interval = mu_interval, n_trees = n_trees, phylo_name = phylo_name)

  # Print each parameter with its name and value
  
  for (param_name in names(parameters)) {
    print(paste(param_name, ":", parameters[param_name]))
  }

  cat("Extracting trees information...")

  # Extract parameters
  parameters <- results$param

  # Print the number of simulations attempted and completed
  num_simulations <- length(results$trees)
  num_completed <- sum(!sapply(results$loglik_estimation, inherits, "try-error"))


  # Assuming each tree is a data frame with columns 'n' and 't_ext'
  count_extinct_species <- function(tree) {
    extinction_times = as.data.frame(tree)["t_ext"]
    sum(extinction_times < 1e+10)
  }


  completed_indices <- which(!sapply(results$loglik_estimation, inherits, "try-error"))
  
  trees = results$trees[completed_indices]

  num_extinct_per_tree <- sapply(trees, count_extinct_species)


  cat("Number of simulations attempted:", num_simulations, "\n")
  cat("Number of simulations completed successfully:", num_completed, "\n")


  # Print errors encountered during simulations
  #errors <- sapply(results$loglik_estimation, function(x) if(inherits(x, "try-error")) x else NA)
  #cat("Errors encountered during simulations:\n", paste(errors, collapse = "\n"))


  # Extract indices of successful simulations
  
  params_completed <- data.frame(mu = results$param$mu,
  lambda = results$param$lambda,
  betaN = results$param$betaN,
  betaP = results$param$betaP,
  loglik = unlist(results$loglik_estimation[completed_indices]),
  nspecies =  num_extinct_per_tree
  )

  cat("Creating plots, this might take a while...")

# Generate scatterplot for mu and lambda, colored by log likelihood estimation
plot1 = ggplot(params_completed, aes(x = mu, y = lambda, color = loglik)) +
  geom_point(alpha = 0.2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "log likelihood estimation", x = "mu", y = "lambda") +
  theme_minimal()

# Generate scatterplot for betaN and betaP, colored by log likelihood estimation
plot2 = ggplot(params_completed, aes(x = betaN, y = betaP, color = loglik)) +
  geom_point(alpha = 0.2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "log likelihood estimation", x = "betaN", y = "betaP") +
  theme_minimal()



# Generate scatterplot for mu and lambda, colored by log likelihood estimation
plot3 = ggplot(params_completed, aes(x = mu, y = lambda, color = nspecies)) +
  geom_point(alpha = 0.2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Number of extincted species", x = "mu", y = "lambda") +
  theme_minimal()

# Generate scatterplot for betaN and betaP, colored by log likelihood estimation
plot4 = ggplot(params_completed, aes(x = betaN, y = betaP, color = nspecies)) +
  geom_point(alpha = 0.2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Number of extincted species", x = "betaN", y = "betaP") +
  theme_minimal()


library(patchwork)


# Combine the plots
combined_plots <- (plot1 | plot2) / (plot3 | plot4)

# Print the combined plot
print(combined_plots)


# Generate scatterplot for betaN and betaP, colored by log likelihood estimation
plot5 = ggplot(params_completed, aes(x = nspecies, y =  loglik)) +
  geom_point(alpha = 0.2) +
  theme_minimal()

print(plot5)
#

}

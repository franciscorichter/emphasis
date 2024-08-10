
# Install and load the rotl package
#install.packages("rotl")
library(rotl)
library(emphasis)
tree <- sim_tree_pd_cpp(pars = c(0.1,0.6,-0.0175,0.001),max_t = 20)

plot(tree$tes)

# Perform the DE learning
result <- learn_de(phylo = tree$tes,num_points = 10,max_iterations = 1000)





library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Convert timestamps to seconds
times <- as.POSIXct(result$times, format="%Y-%m-%d %H:%M:%S")
start_time <- min(times)
time_seconds <- as.numeric(difftime(times, start_time, units = "secs"))

# Ensure lengths match for plotting
min_length <- min(length(time_seconds), length(result$minloglik), length(result$meanloglik))

# Adjust vectors to have the same length
time_seconds <- time_seconds[1:min_length]
minloglik <- result$minloglik[1:min_length]
meanloglik <- result$meanloglik[1:min_length]

# Initialize data frames for each parameter
param_names <- colnames(result$parameters[[1]]$pars)
param_stats <- list()
for (param in param_names) {
  param_stats[[param]] <- data.frame(
    time_seconds = numeric(),
    mean = numeric(),
    min = numeric(),
    max = numeric(),
    q25 = numeric(),
    q50 = numeric(),
    q75 = numeric()
  )
}

# Extract statistics using a for loop
for (i in 1:length(result$parameters)) {
  pars <- result$parameters[[i]]$pars
  for (param in param_names) {
    values <- pars[, param]
    stats <- data.frame(
      time_seconds = time_seconds[i],
      mean = mean(values),
      min = min(values),
      max = max(values),
      q25 = quantile(values, 0.25),
      q50 = quantile(values, 0.5),
      q75 = quantile(values, 0.75)
    )
    param_stats[[param]] <- rbind(param_stats[[param]], stats)
  }
}

# Plotting function for parameters
plot_stats <- function(param_stats, parameter) {
  ggplot(param_stats, aes(x = time_seconds)) +
    geom_ribbon(aes(ymin = min, ymax = max), fill = "grey80", alpha = 0.5) +
    geom_line(aes(y = q25), linetype = "dashed", size = 1) +
    geom_line(aes(y = q75), linetype = "dashed", size = 1) +
    geom_line(aes(y = mean), linetype = "solid", size = 1.2) +
    labs(title = paste("Evolution of Parameter", parameter), x = "Time (seconds)", y = "Value") +
    theme_minimal()
}

# Generate and display plots for each parameter
plots <- lapply(param_names, function(param) plot_stats(param_stats[[param]], param))

# Plot minloglik and meanloglik over time
plot_minloglik <- ggplot(data.frame(time_seconds, minloglik), aes(x = time_seconds, y = minloglik)) +
  geom_line() +
  labs(title = "Min Log Likelihood over Time", x = "Time (seconds)", y = "Min Log Likelihood") +
  theme_minimal()

plot_meanloglik <- ggplot(data.frame(time_seconds, meanloglik), aes(x = time_seconds, y = meanloglik)) +
  geom_line() +
  labs(title = "Mean Log Likelihood over Time", x = "Time (seconds)", y = "Mean Log Likelihood") +
  theme_minimal()

# Plot min_pars
plot_min_pars <- ggplot(data.frame(min_par1 = result$min_pars[,1], min_par2 = result$min_pars[,2]), aes(x = min_par1, y = min_par2)) +
  geom_point() +
  labs(title = "Min Parameters", x = "Parameter 1", y = "Parameter 2") +
  theme_minimal()

plot_min_pars2 <- ggplot(data.frame(min_par3 = result$min_pars[,3], min_par4 = result$min_pars[,4]), aes(x = min_par3, y = min_par4)) +
  geom_point() +
  labs(title = "Min Parameters", x = "Parameter 3", y = "Parameter 4") +
  theme_minimal()

# Arrange the plots in a grid and display them
grid.arrange(grobs = c(plots, list(plot_minloglik, plot_meanloglik, plot_min_pars,plot_min_pars2)), ncol = 2)



library(DDD)

cat("Estimating the intrinsic speciation rate lambda and the carrying capacity K")
cat("for a fixed extinction rate of 0.1, conditioning on clade survival and two missing species:")
tree = dd_sim(pars = c(0.8,0.1,40),age = 10)
dd_ML(brts = tree$brts)
plot(tree$tes)
result <- learn_de(phylo = tree$tes,num_points = 50,max_iterations = 1000,lower_bound = c(0,0,-0.1,0),upper_bound = c(1,5,0,0))

result$obtained_estim


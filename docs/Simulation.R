## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(emphasis) # Load the Emphasis package

## ----nhpp-sampling------------------------------------------------------------
# Define a rate function for the NHPP
rate_function_example <- function(t) { 0.5 + 0.7*sin(t) }

# Set parameters for sampling
n_samples = 5000   # Number of event times to sample
now = 0          # Start time
tMax = 20        # Maximum time limit

# Sample event times from the NHPP
nhpp_event_times <- nhExpRand(n_samples, rate_function_example, now, tMax)

# Print the sampled event times
#print(nhpp_event_times)

# Density histogram of sampled event times
hist(nhpp_event_times, breaks = 100, col = 'skyblue', freq = FALSE,
     main = 'Density Histogram of NHPP Sampled Event Times',
     xlab = 'Event Time', ylab = 'Density')


# Adding a line plot for the rate function
curve(rate_function_example, from = now, to = tMax, add = TRUE, col = 'red', lwd = 2)
legend('topright', legend = 'Rate Function', col = 'red', lwd = 2)


## ----rate_t-function----------------------------------------------------------
# Define example covariate functions
cov_func1 <- function(t) { sin(t) }
cov_func2 <- function(t) { cos(t) }

# Example parameters (baseline and coefficients for covariates)
params <- c(0.5, 1.2, -0.8)


# Generate a sequence of time points
time_points <- seq(0, 10, by = 0.1)

# Compute covariate values and rate at each time point
cov_values1 <- sapply(time_points, cov_func1)
cov_values2 <- sapply(time_points, cov_func2)
rate_values <- sapply(time_points, function(t) rate_t(t, params, 
                              list(cov_func1, cov_func2), use_exponential = TRUE))


# Plotting the covariate functions and rate function
plot(time_points, cov_values1, type = 'l', col = 'blue', ylim = range(c(cov_values1, cov_values2, rate_values)),
     main = 'Covariates and Rate Function', xlab = 'Time', ylab = 'Values / Rate', lty = 2)
lines(time_points, cov_values2, col = 'green', lty = 2)
lines(time_points, rate_values, col = 'red') # Using a different line type for rate

# Adding a legend
legend('topright', legend = c('Covariate Function 1', 'Covariate Function 2', 'Rate Function'),
       col = c('blue', 'green', 'red'), lty = c(1, 1, 2), cex = 0.8)



## -----------------------------------------------------------------------------

# Install and load the required packages
if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}
library(ape)

# Example parameter set
pars <- params <- c(mu = 0.1, lambda_0 = 0.5, beta_N = -0.02, beta_P = 0.03)
max_t <- 10
max_lin <- 1e6
max_tries <- 100

tree = emphasis:::sim_tree_pd_cpp(pars, max_t, max_lin, max_tries)$tes

time_points = seq(0,max_t,length.out = 100)
# Compute species richness (number of lineages) over time
richness_values <- sapply(time_points, function(t) sum(ape:::branching.times(tree) < t))

# Compute phylogenetic diversity over time
diversity_values <- sapply(time_points, function(t) sum(ape:::branching.times(tree)[ape:::branching.times(tree) < t]))

# Compute speciation rate
speciation_rate_values <- sapply(time_points, function(t) {
  Nt <- richness_values[t]
  Pt <- diversity_values[t]
  lambda_t <- max(0, params['lambda_0'] + params['beta_N'] * Nt + params['beta_P'] * Pt / max(Nt, 1))
  return(lambda_t)
})

# Plotting the results
plot(time_points, richness_values, type = 'l', col = 'blue', lty = 2, 
     ylim = range(c(richness_values, diversity_values, speciation_rate_values)),
     main = 'PD Model Visualization', xlab = 'Time', ylab = 'Values / Rate')
lines(time_points, diversity_values, col = 'green', lty = 2)
lines(time_points, speciation_rate_values, col = 'red', lty = 1) # Solid line for speciation rate

# Add a legend
legend('topright', legend = c('Species Richness (Covariate)', 'Phylogenetic Diversity (Covariate)', 'Speciation Rate'),
       col = c('blue', 'green', 'red'), lty = c(2, 2, 1), cex = 0.8)





## -----------------------------------------------------------------------------


  
plot(tree)

larger_tree = tree
# Define time points for visualization, ensuring they don't exceed the tree age
max_time <- max(branching.times(larger_tree))
time_points <- seq(0, max_time, length.out = 100)

# Compute species richness and phylogenetic diversity at each time point
richness_values <- sapply(time_points, function(t) sum(branching.times(larger_tree) < t))
diversity_values <- sapply(time_points, function(t) sum(branching.times(larger_tree)[branching.times(larger_tree) < t]))

# Compute the speciation rate at each time point
params <- c(lambda_0 = 0.5, beta_N = 1.2, beta_P = -0.8)
speciation_rate_values <- sapply(time_points, function(t) {
  Nt <- richness_values[t]
  Pt <- diversity_values[t]
  lambda_t <- max(0, params['lambda_0'] + params['beta_N'] * Nt + params['beta_P'] * Pt / max(Nt, 1))
  return(lambda_t)
})

# Plot species richness, phylogenetic diversity, and speciation rate
plot(time_points, richness_values, type = 'l', col = 'blue', lty = 2, 
     ylim = range(c(richness_values, diversity_values, speciation_rate_values)),
     main = 'PD Model Visualization on Larger Tree', xlab = 'Time', ylab = 'Values / Rate')
lines(time_points, diversity_values, col = 'green', lty = 2)
lines(time_points, speciation_rate_values, col = 'red', lty = 1)

# Add a legend
legend('topright', legend = c('Species Richness (Covariate)', 'Phylogenetic Diversity (Covariate)', 'Speciation Rate'),
       col = c('blue', 'green', 'red'), lty = c(2, 2, 1), cex = 0.8)

## -----------------------------------------------------------------------------

pars =c(0.21779984 ,1.16601633, -0.08090578 , 0.01261866)
tree = sim_tree_pd_cpp(pars = pars,max_t = 5,max_lin = 1e+7,max_tries = 1)
phylo = tree$tes
plot(phylo)
brts <- ape::branching.times(phylo)

# Define parameter bounds
lower_bound = c(0,0.5,-0.1,0)
upper_bound = c(1,3,0,0.1)

num_points <- 1000

pars <- matrix(nrow = num_points, ncol = length(lower_bound))
for (i in seq_along(lower_bound)) {
  pars[, i] <- runif(num_points, min = lower_bound[i], max = upper_bound[i])
}

dmval <- mcGrid(pars, 
                brts = brts, 
                sample_size = 1, 
                maxN = 1, 
                soc = 2, 
                max_missing = 1e+4, 
                max_lambda = 1e+4, 
                lower_bound = lower_bound, 
                upper_bound = upper_bound, 
                xtol_rel = 0.01, 
                num_threads = 4)
# Extract relevant columns
loglikelihood_est <- dmval[, 1]
sampled_x_loglikelihood <- dmval[, 5]
sampling_probability <- dmval[, 6]

# Remove rows with NA in loglikelihood_est
valid_indices <- !is.na(loglikelihood_est)
loglikelihood_est <- loglikelihood_est[valid_indices]
sampled_x_loglikelihood <- sampled_x_loglikelihood[valid_indices]
sampling_probability <- sampling_probability[valid_indices]

# Plotting the distribution of loglikelihood estimates
hist(sampling_probability/loglikelihood_est, main="Distribution of Loglikelihood Estimates", xlab="Loglikelihood Estimate")

# Scatter plot of sampled_x_loglikelihood vs sampling_probability
plot(sampled_x_loglikelihood, sampling_probability, main="Sampled Loglikelihood vs Sampling Probability",
     xlab="Sampled Loglikelihood", ylab="Sampling Probability", pch=19)

# Adding a line to visualize the ideal sampling distribution might require additional theoretical information
# For simplicity, this example focuses on the distributions and relationships in the sampled data


# Assuming dmval is your matrix
plot(dmval[,5], dmval[,6], main = "Loglikelihood vs. Sampling Probability",
     xlab = "Loglikelihood of Sampled x", ylab = "Sampling Probability",
     pch = 19, col = "blue")

# Adding a 1:1 line for perfect agreement
abline(a = 0, b = 1, col = "red")

# Calculate correlation
correlation <- cor(dmval[,5], dmval[,6], use = "complete.obs")
cat("Correlation between Loglikelihood of Sampled x and Sampling Probability: ", correlation, "\n")


## -----------------------------------------------------------------------------
# load simulation data
#path_data = "~/Library/CloudStorage/Dropbox/Pancho/78 - GAM paper/01 - Simulation data/Parulidae_20240205_233041.RData"
#process_and_plot(path_data)


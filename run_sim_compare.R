library(emphasis)
library(ggplot2)
library(microbenchmark)
library(ape)

# Example parameters (edit if needed)
pars <- c(0.1, 0.05, 0, 0, 0, 0) # adjust as needed
max_t <- 10
max_N <- 1000
max_tries <- 10
n_trees <- 1000

# R function wrapper (replace with your actual R simulation function if different)
simulate_tree_R <- function() {
  # Example: emphasis::simulate_tree(model = list(type = "birth-death", lambda = 0.1, mu = 0.05), n_taxa = 0, max_t = max_t)
  # Replace the line below with your actual R simulation call
  ape::rtree(n = sample(10:100, 1)) # placeholder for demonstration
}

# C++ function wrapper (using emphasis package)
simulate_tree_C <- function() {
  res <- emphasis::simulate_div_tree_cpp(pars, max_t, max_N, max_tries)
  tree_mat <- res$tree
  # Get number of unique tip labels (assuming column 3 is tip label)
  n_species <- length(unique(tree_mat[, 3]))
  return(n_species)
}

set.seed(123)
n_trees <- 5  # Reduce for debugging

cat("Simulating R trees...\n")
r_times <- system.time({
  r_trees <- replicate(n_trees, simulate_tree_R(), simplify = FALSE)
})
r_n_species <- sapply(r_trees, function(tr) ape::Ntip(tr))

cat("Simulating C++ trees...\n")
safe_simulate_tree_C <- function() {
  tryCatch({
    res <- emphasis::simulate_div_tree_cpp(pars, max_t, max_N, max_tries)
    tree_mat <- res$tree
    n_species <- length(unique(tree_mat[, 3]))
    return(n_species)
  }, error = function(e) {
    message("Error: ", e$message)
    return(NA)
  })
}
c_times <- system.time({
  c_n_species <- replicate(n_trees, safe_simulate_tree_C())
})

df <- data.frame(
  n_species = c(r_n_species, c_n_species),
  method = rep(c("R", "C++"), each = n_trees)
)

# Plot distribution of number of species
p1 <- ggplot(df, aes(x = n_species, fill = method)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Number of Species", x = "Number of Species", y = "Density")

timing_df <- data.frame(
  method = c("R", "C++"),
  time = c(r_times["elapsed"], c_times["elapsed"])
)
p2 <- ggplot(timing_df, aes(x = method, y = time, fill = method)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Simulation Time (1000 trees)", y = "Seconds")

print(p1)
print(p2)

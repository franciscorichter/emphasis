library(ape)
library(tidyverse)

# Function to compute the phylogenetic distance matrix
compute_phylo_distance_matrix <- function(tree) {
  # Get the terminal nodes (leaves)
  terminals <- tree$tip.label
  
  # Initialize a distance matrix with zeros
  n <- length(terminals)
  dist_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Populate the distance matrix
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Calculate the phylogenetic distance between the two terminals
      phylo_dist <- cophenetic(tree)[terminals[i], terminals[j]]
      
      # Update the distance matrix
      dist_matrix[i, j] <- phylo_dist
      dist_matrix[j, i] <- phylo_dist
    }
  }
  
  # Create a DataFrame for better visualization
  df <- as.data.frame(dist_matrix)
  colnames(df) <- terminals
  rownames(df) <- terminals
  
  return(df)
}

# Create an example tree
tree_newick <- "(((A:1,B:1):1,C:2):1,(D:1.5,(E:0.5,F:0.5):1):1.5);"
tree <- read.tree(text = tree_newick)
plot(tree)
# Compute the phylogenetic distance matrix
df <- compute_phylo_distance_matrix(tree)
print(df / 2)

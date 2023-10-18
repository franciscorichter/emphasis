from Bio import Phylo
from io import StringIO
import numpy as np
import pandas as pd

# Function to compute the phylogenetic distance matrix
def compute_phylo_distance_matrix(tree):
    # Get the terminal nodes (leaves)
    terminals = tree.get_terminals()
    
    # Initialize a distance matrix with zeros
    n = len(terminals)
    dist_matrix = np.zeros((n, n))
    
    # Populate the distance matrix
    for i in range(n):
        for j in range(i+1, n):
            # Get the most recent common ancestor of the two terminals
            mrca = tree.common_ancestor(terminals[i], terminals[j])
            
            # Calculate the phylogenetic distance between the two terminals
            dist_i_to_mrca = tree.distance(terminals[i], mrca)
            dist_j_to_mrca = tree.distance(terminals[j], mrca)
            phylo_dist = dist_i_to_mrca + dist_j_to_mrca
            
            # Update the distance matrix
            dist_matrix[i, j] = phylo_dist
            dist_matrix[j, i] = phylo_dist
    
    # Create a DataFrame for better visualization
    terminal_names = [t.name for t in terminals]
    df = pd.DataFrame(dist_matrix, index=terminal_names, columns=terminal_names)
    
    return df

# Create an example tree
tree_newick = "(((A:1,B:1):1,C:2):1,D:3);"
tree = Phylo.read(StringIO(tree_newick), "newick")

# Compute the phylogenetic distance matrix
df = compute_phylo_distance_matrix(tree)
print(df/2)

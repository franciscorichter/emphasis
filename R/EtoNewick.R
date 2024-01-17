toNewick <- function(E, present_time) {
  # Initialize tokens for each species
  species_names <- LETTERS[1:nrow(E)]
  tokens <- setNames(as.list(species_names), 1:lespecies_names)
  
  bt = sort(c(E[E!=0]),decreasing = T)
  # Iterate over each time point from present to past
  for (i in 1:length(bt)) {
    # Check for speciation events at this time point
    relationship <- which(E == bt[i], arr.ind = TRUE)
    token_parent = tokens[[relationship[1,1]]]
    token_child = tokens[[relationship[1,2]]]
    token_time = bt
    paste0("(", token_parent,":",present_time-bt[i],",",token_child,")")
    
    for (event in seq_len(nrow(events))) {
      parent_idx <- events[event, "row"]
      child_idx <- events[event, "col"]
      
      parent <- species_names[parent_idx]
      child <- species_names[child_idx]
      
      # Merge child into parent token
      child_time <- present_time - time
      parent_time <- ifelse(E[parent_idx, parent_idx] == 0, present_time, present_time - E[parent_idx, parent_idx])
      tokens[[parent]] <- sprintf("(%s:%s,%s:%s)", tokens[[parent]], parent_time, tokens[[child]], child_time)
      
      # Remove child token
      tokens[[child]] <- NULL
    }
  }
  
  # The final Newick string
  newick_string <- paste(tokens, collapse = "", sep = "")
  return(paste(newick_string, ";", sep = ""))
}



# Example usage
E <- matrix(c(0, -1, 0, 2, 0,
              0,  0, 1, 0, 0,
              0,  0, 0, 0, 3,
              0,  0, 0, 0, 0,
              0,  0, 0, 0, 0), byrow = TRUE, nrow = 5)

newick_string <- toNewick(E, present_time = 4)
print(newick_string)

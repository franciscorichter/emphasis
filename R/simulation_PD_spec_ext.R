simulate_evolution <- function(N_max, spec_pars, ext_pars, max_t, ...) {
  
  tree_mat <- matrix(NA, ncol = 4, nrow = 0)
  colnames(tree_mat) <- c("birth.time", "parent", "child", "death.time")
  
  tips_at_t = c(1, 2)  # Start with two species
  counter = 2  # Counter for the next species ID
  t = 0  # Simulation start time
  
  while (length(tips_at_t) <= N_max && t <= max_t) {
    D = diversity(E)

    spe_rates = ext_rates = NULL
    for(species in 1:nrow(E)){
      spec_rate = rate_diversity(pars = spe_pars,covariates = c(D$Richness,D$Phylogenetic_Diversity,mean(D$Distance_Matrix[i,]))) 
      spec_rates = c(spec_rates,spec_rate)
      
      ext_rate = rate_diversity(pars = ext_pars,covariates = c(D$Richness,D$Phylogenetic_Diversity,mean(D$Distance_Matrix[i,]))) 
      ext_rates = c(ext_rates,ext_rate)
    }
    total_rate <- function(t) sum(spe_rates + ext_rates)
    t_sample = nhExpRand(1, total_rate)
    event = sample(c(0,1),1,prob=c(sum(ext_rates)/sum(ext_rates+spe_rates),sum(spe_rates)/sum(ext_rates+spe_rates)))
    if(event==0){#extinction
      species = sample(1:nrow(E),prob=ext_rates)
      E = updateE(E=E,event="speciation",species=species)
    }else{ #speciation
      parent = sample(1:nrow(E),prob=ext_rates)
      E = updateE(E=E,event="extinction",species=parent)
    }
      
    n <- nrow(E)
    I_n <- diag(1, n)  # Identity matrix
    D <- t_p * (matrix(1, n, n) - I_n) - E
   
    
  }
  
  return(list(tree_mat = tree_mat, metrics = metrics))
}


sample_and_update_event <- function(tree_mat, tips_at_t, counter, spec_rates, ext_rates) {
  total_rate = sum(spec_rates) + sum(ext_rates)
  if (runif(1) < sum(spec_rates) / total_rate) {
    # Speciation event
    event_species = sample(tips_at_t, 1, prob = spec_rates)
    new_species_id = counter + 1
    tree_mat = update_tree_mat_speciation(tree_mat, event_species, new_species_id)
    tips_at_t = c(tips_at_t, new_species_id)
    counter = new_species_id
  } else {
    # Extinction event
    event_species = sample(tips_at_t, 1, prob = ext_rates)
    tree_mat = update_tree_mat_extinction(tree_mat, event_species)
    tips_at_t = tips_at_t[tips_at_t != event_species]
  }
  return(list(tree_mat = tree_mat, tips_at_t = tips_at_t, counter = counter))
}

sample_next_event_time <- function(current_time, spec_rates, ext_rates) {
  total_rate = sum(spec_rates) + sum(ext_rates)
  t_sample = nhExpRand(1, total_rate)
  return(current_time + t_sample)
}



diversity <- function(E,t_present) {
  # Extract present time (t_p) as the maximum value in E
  t_p <- max(E)
  
  # Species Richness (R)
  R <- nrow(E)
  
  # Phylogenetic Diversity (P)
  # Summing all positive elements of E
  P <- sum(E) + t_p
  
  # Phylogenetic Distance Matrix (D)
  n <- nrow(E)
  I_n <- diag(1, n)  # Identity matrix
  D <- t_p * (matrix(1, n, n) - I_n) - E
  
  return(list("Richness" = R, "Phylogenetic_Diversity" = P, "Distance_Matrix" = D))
}

# Function to handle speciation
handle_speciation <- function(E, parent_species, time) {
  # Add new species to the matrix
  n <- nrow(E)
  E <- rbind(E, rep(0, n))
  E <- cbind(E, rep(0, n + 1))
  
  # Record the speciation event
  E[parent_species, n + 1] <- time
  return(E)
}

# Function to handle extinction
handle_extinction <- function(E, extinct_species) {
  # Find daughters of the extinct species
  daughters_species <- which(E[extinct_species, ] > 0)
  
  # If the extinct species has daughters, reassign them to the MCA
  if (length(daughters_species) > 0) {
    parent_species <- which(E[, extinct_species] > 0)
    parent_time <- E[parent_species, extinct_species]
    daughters_times <- E[extinct_species, daughters_species]
    
    # Determine the most recent common ancestor
    recent_time <- max(c(parent_time, daughters_times))
    mca <- c(parent_species, daughters_species)[which(recent_time == c(parent_time, daughters_times))]
    
    E[mca, daughters_species] <- E[extinct_species, daughters_species]
  }
  
  # Remove extinct species from the matrix
  E <- E[-extinct_species, -extinct_species]
  return(E)
}


rate_diversity <- function(pars, covariates, g = function(x) x){
  rate =  g(sum(pars * c(1,covariates)))
  return(rate)
}


updateE_extant <- function(E,event,species,event_time){
  if(event == "speciation"){
    new_species = rep(0,nrow(E))
    new_species[species] = event_time
    E = cbind(E,new_species)
    E = rbind(E,E[,ncol(E)])
  }
  if(event == "extinction"){
    ext_species = E[species,]
    heritage = which(E[species,]>0)
    if(is.null(heritate)){
      E = E[-species,-species]
    }else{
      mca_time = max(E[species,])
      mca = which(E[species,]==mca_time)
      E[mca,heritage] = E[species,heritage]
      E = E[-species,-species]
    }
  }
  return(E)
}


updateE_extant <- function(E, event, species, event_time) {
  if (event == "speciation") {
    # Add a new column for the new species
    new_species = rep(0, nrow(E))
    new_species[species] = event_time
    E = cbind(E, new_species)
    # Add a new row representing the new species
    new_species_row = rep(0, ncol(E))
    new_species_row[ncol(E)] = event_time
    E = rbind(E, new_species_row)
  } else if (event == "extinction") {
    # Handling extinction
    heritage = which(E[species, ] > 0)
    if (length(heritage) == 0) {
      # If no heritage, remove the species row and column
      E = E[-species, -species]
    } else {
      # Find the most common ancestor and update its heritage
      mca_time = max(E[species, heritage])
      mca = which(E[, species] == mca_time)
      E[mca, heritage] = E[species, heritage]
      E = E[-species, -species]
    }
  }
  return(E)
}

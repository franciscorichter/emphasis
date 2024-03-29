#' @keywords internal
calc_p <- function(l_table, t) {
  l_table[, 1] <- t - l_table[, 1]
  return(sum(DDD::L2phylo(l_table, dropextinct = T)$edge.length))
}

#' simulation function to simulate a tree under the pd model, returns the full
#' tree.
#' @description super fast function to simulate the process of diversification
#' with diversity dependence and phylogenetic diversity dependence
#' @param pars parameter vector with c(mu, lambda_0, beta_N, beta_P)
#' @param max_t crown age
#' @return a list with the phylogeny, and a vector with population size,
#' time of extinction (equal to the crown #' age in the absence of extinction),
#' and the phylogenetic diversity at the time of extinction (or crown age).
#' @export
sim_tree_pd_R <- function(pars, max_t) {
  N1 <- 1
  N2 <- 1
  N <- N1 + N2
  t <- 0
  # birth date, parent label, ID, time of extinction, clade
  tree <- matrix(nrow = 2, ncol = 4)
  
  tree[1, ] <- c(0, 0, -1, -1)
  tree[2, ] <- c(0, -1, 2, -1)
  tree_ID <- 3
  
  mu <- pars[1]
  P <- 0
  while (t < max_t && 
         N1 >= 1 && N2 >= 1) {
    
    N <- N1 + N2
    spec_rate <- max(0, pars[2] + 
                       pars[3] * N  +  
                       ((P + N * (max_t - t) - t) / N) * pars[4])
    total_rate <- (spec_rate + mu) * N
    
    if (total_rate == 0) {
      t <- max_t
      break
    }
    
    next_event_time <- t + stats::rexp(n = 1, rate = total_rate)  

    P <- calc_p(tree, t)
    
    if (next_event_time < max_t) {
      focal_spec <- max(0, pars[2] + 
                          pars[3]*N  +  
                          ((P + N * (next_event_time - t) - t) / N) * pars[4])
      
      pt = ((focal_spec + mu) * N ) / total_rate
      
      if (stats::runif(1) < pt) {
        # event is accepted
        
        # pick event type:
        if (stats::runif(1) < focal_spec / (focal_spec + mu)) {
          # speciation
         # cat("speciation", "\n")
          parent <- sample(which(tree[, 4] == -1), 1)
          
          new_ID <- tree_ID
          tree_ID <- tree_ID + 1
          if (tree[parent, 3] < 0) new_ID <- new_ID * -1  # 'b' clade
          
          new_spec <- c(next_event_time, tree[parent, 3], new_ID, -1)
          
          tree <- rbind(tree, new_spec)
          
          if (new_ID < 0) {
            N2 <- N2 + 1
          } else {
            N1 <- N1 + 1
          }
        } else {
          # extinction
          # pick random species
         # cat("extinction", "\n")
          to_remove <- sample(which(tree[, 4] == -1), 1)
          
          tree[to_remove, 4] <- next_event_time
          
          if (tree[to_remove, 3] < 0) {
            N2 <- N2 - 1
          } else {
            N1 <- N1 - 1
          }
        }
      }
    }
    t <- next_event_time
  }
  
  if (N1 >= 1 && N2 >= 1) P <- calc_p(tree, t)
  N <- N1 + N2
  
  tree[, 1] <- max_t - tree[, 1]
  extinct <- which(tree[, 4] != -1)
  tree[extinct, 4] <- max_t - tree[extinct, 4]
  
   for_ddd <- as.matrix(tree[, 1:4])
  
   phy <- DDD::L2phylo(as.matrix(tree[, 1:4]), dropextinct = FALSE)
  
  t <- min(t, max_t)
  
  return(list("phy" = phy, 
              "result" = c(N, t, P)))
  # return(list("result" = c(N, t, P)))
}




#' simulation function to simulate a tree under the pd model, returning whether
#' the tree went extinct before max_t or not. 
#' This function does not return a phylogenetic tree to improve computation
#' speed. 
#' @description super fast function to simulate the process of diversification
#' with diversity dependence and phylogenetic diversity dependence
#' @param pars parameter vector with c(mu, lambda_0, beta_N, beta_P)
#' @param max_t crown age
#' @param num_repl number of replicates
#' @param max_lin number of lineages past which non-extinction is assumed.
#' @return a tibble with five columns: 1) whether or not the tree went extinct,
#' 2) the time of extinction (equal to the crown #' age in the absence of 
#' extinction), 3) the number of tips,4) the phylogenetic diversity at the 
# time of extinction (or crown age), and 5) the break condition indicating why 
# the simulation was stopped (options: no break (none), exceeded max_t, 
# extinction, or exceeded max_lin).
#' @export
sim_tree_is_extinct_pd <- function(pars, max_t, num_repl = 1, max_lin) {
    result <- simulate_pd_trees_cpp(pars, max_t, num_repl, max_lin)
    colnames(result) <- c("is_extinct", "t", "N", "P", "break_condition")
    result <- tibble::as_tibble(result)
    result$break_condition[result$break_condition == 3] <- "max_lin"
    result$break_condition[result$break_condition == 2] <- "extinction"
    result$break_condition[result$break_condition == 1] <- "max_t"
    result$break_condition[result$break_condition == 0] <- "none"
    
    return(result)
}

#' simulation function to simulate a tree under the pd model
#' @description super fast function to simulate the process of diversification
#' with diversity dependence and phylogenetic diversity dependence
#' @param pars parameter vector with c(mu, lambda_0, beta_N, beta_P)
#' @param max_t crown age
#' @param max_lin number of lineages past which non-extinction is assumed.
#' @param max_tries maximum number of tries to get a non-extinct tree.
#' @return list with: 1) tes - reconstructed tree, 2) tas - tree with extinct
#' lineages, 3) L = Ltable and 4) brts - branching times of the reconstructed
#' tree and 5) status of simulation, options 1) "extinct", 2) "too_large" or 3)
#' "done".
#' @export
sim_tree_pd_cpp <- function(pars,
                            max_t,
                            max_lin = 1e6,
                            max_tries = 100, 
                            useDDD=TRUE) {

  result <- simulate_single_pd_tree_cpp(pars,
                                        max_t,
                                        max_lin,
                                        max_tries)
  tes <- NULL
  tas <- NULL
  brts <- NULL
  
  if (result$status == "extinct") {
    warning("could not simulate tree, all trees went extinct, try increasing max_tries")
  }
  if (result$status == "too_large") {
    warning("could not simulate tree, all trees were too large, try increasing max_lin")
  }
  if (result$status == "done" & useDDD) {
    tes <- DDD:::L2phylo(result$Ltable, dropextinct = TRUE)
    tas <- DDD:::L2phylo(result$Ltable, dropextinct = FALSE)
    brts = DDD:::L2brts(result$Ltable, dropextinct = TRUE)
  }
  
  out = list(tes = tes, tas = tas, L = result$Ltable, brts = brts,
             status = result$status)
  
  return(out)
}


#' simulation function to simulate a tree under the pd model
#' @description super fast function to simulate the process of diversification
#' with diversity dependence and phylogenetic diversity dependence
#' @param mu_vec vector of extinction rates to explore
#' @param lambda_vec vector of lambda rates to explore
#' @param b_n_vec vector of B_n rates to explore
#' @param b_p_vec vector of B_p rates to explore
#' @param max_t crown age
#' @param max_N maximum number of lineages past which non-extinction is assumed.
#' @param num_repl number of replicates
#' @return a list with the used parameter combinations, 
#' extinction of the phylogeny, and a vector with population size, 
#' time of extinction (equal to the crown #' age in the absence of extinction), 
#' and the phylogenetic diversity at the time of extinction (or crown age).
#' @export
sim_tree_pd_grid <- function(mu_vec,
                             lambda_vec,
                             b_n_vec,
                             b_p_vec,
                             max_t,
                             num_repl,
                             max_N) {
  result <- explore_grid_cpp(mu_vec,
                             lambda_vec,
                             b_n_vec,
                             b_p_vec,
                             max_t,
                             num_repl,
                             max_N)
  colnames(result) <- c("is_extinct", "t", "N", "P")
  result <- tibble::as_tibble(result)
  return(result)
}



#' Simulate a single phylogenetic diversity tree using C++ backend
#'
#' This function interfaces with a C++ function to simulate a single
#' phylogenetic tree based on the provided parameters and constraints.
#'
#' @param pars A numeric vector of parameters required by the C++ simulation function.
#' @param max_t The maximum time or iterations to run the simulation for.
#' @param max_N The maximum number of nodes that the phylogenetic tree should have.
#' @param max_tries The maximum number of tries to attempt the simulation before stopping.
#'
#' @return A simulated phylogenetic tree object (specific format to be detailed).
#' @export
#' @examples
#' # Example usage:
#' # Define parameters (example parameters to be replaced with actual ones)
#' parameters <- c(0.1, 0.2, 0.3)
#' max_time <- 10000
#' max_nodes <- 100
#' max_attempts <- 10
#'
#' # Run the simulation
#' tree <- simulate_single_pd_tree_cpp(parameters, max_time, max_nodes, max_attempts)
#' 
simulate_single_pd_tree_cpp <- function(pars, max_t, max_N, max_tries) {
  .Call('_emphasis_simulate_single_pd_tree_cpp', PACKAGE = 'emphasis', pars, max_t, max_N, max_tries)
}

simulate_pd_trees_cpp <- function(pars, max_t, repl, max_N) {
  .Call('_emphasis_simulate_pd_trees_cpp', PACKAGE = 'emphasis', pars, max_t, repl, max_N)
}
  
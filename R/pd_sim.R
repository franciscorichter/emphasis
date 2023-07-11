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
#' This function does not return a phylogenetic tree to improve computation speed. 
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
sim_tree_pd_cpp <- function(pars, max_t, num_repl = 1, max_lin) {
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
#' tree
#' @export
sim_single_tree_pd_cpp <- function(pars,
                                   max_t,
                                   max_lin = 1e6,
                                   max_tries = 100) {
  result <- simulate_single_pd_tree_cpp(pars, max_t, max_lin,
                                        max_tries)
  if (nrow(result) < 2) {
    stop("could not simulate tree")
  }
  #phy_object <- DDD::L2phylo(result, dropextinct = dropextinct)
  tes <- DDD::L2phylo(result, dropextinct = TRUE)
  tas <- DDD::L2phylo(result, dropextinct = FALSE)
  brts = DDD::L2brts(result, dropextinct = TRUE)
  
  out = list(tes = tes, tas = tas, L = result, brts = brts)
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
  
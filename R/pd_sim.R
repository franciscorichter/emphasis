#' simulation function to simulate a minimal pd model
#' @description super fast function to simulate the process of diversification
#' with diversity dependence and phylogenetic diversity dependence
#' @param pars parameter vector with c(mu, lambda_0, beta_N, beta_P)
#' @param max_t crown age
#' @return a vector with population size, time of extinction (equal to the crown 
#' age in the absence of extinction), and the phylogenetic diversity at the time
#' of extinction (or crown age).
#' @export
sim_tree_pd_minimal <- function(pars, max_t) {
  N <- 2
  t <- 0
  tree <- c(t, t)
  mu <- pars[1]
  P <- 0
  while (t < max_t && N >= 2 && N < 1000) {
    spec_rate <- max(0, pars[2] + 
                        pars[3] * N  +  
                        ((P + N * (max_t - t) - t) / N) * pars[4])
    total_rate <- (spec_rate + mu) * N
    next_event_time <- t + stats::rexp(n = 1, rate = total_rate)  # max_t - log(x = u1) / total_rate
    P <- P + N * (next_event_time - t)
    
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
          parent <- sample(1:N, 1) # sample a parent
          tree[parent] <- next_event_time      # set the most recent daughter time for that parent
          tree <- c(tree, next_event_time)                  # add the offspring
          N <- N + 1
        } else {
          # extinction
          # pick random species
          to_remove <- sample(1:N, 1)
          P <- P - (t - tree[to_remove]) # we have to correct P for the extinct species 
          tree <- tree[-to_remove]
          N <- N - 1
        }
      }
    }
    t = next_event_time
  }
  
  t <- min(t, max_t)
  
  return(c(N, t, P))
}

#' simulation function to simulate a tree under the pd model
#' @description super fast function to simulate the process of diversification
#' with diversity dependence and phylogenetic diversity dependence
#' @param pars parameter vector with c(mu, lambda_0, beta_N, beta_P)
#' @param max_t crown age
#' @return a list with the phylogeny, and a vector with population size, 
#' time of extinction (equal to the crown #' age in the absence of extinction), 
#' and the phylogenetic diversity at the time of extinction (or crown age).
#' @export
sim_tree_pd_full <- function(pars, max_t) {
  N <- 2
  t <- 0
  tree <- matrix(nrow = 2, ncol = 5)  # birth date, parent label, ID, time of extinction, time of last daughter
  
  tree[1, ] <- c(0, 0, -1, -1, 0)
  tree[2, ] <- c(0, -1, 2, -1, 0)
  tree_ID <- 3
  
  mu <- pars[1]
  P <- 0
  while (t < max_t && N >= 2 && N < 1000) {
    spec_rate <- max(0, pars[2] + 
                        pars[3] * N  +  
                        ((P + N * (max_t - t) - t) / N) * pars[4])
    total_rate <- (spec_rate + mu) * N
    next_event_time <- t + stats::rexp(n = 1, rate = total_rate)  # max_t - log(x = u1) / total_rate
    P <- P + N * (next_event_time - t)
    
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
          parent <- sample(which(tree[, 4] == -1), 1)
          
          new_ID <- tree_ID
          tree_ID <- tree_ID + 1
          if (tree[parent, 3] < 0) new_ID <- new_ID * -1
          new_spec <- c(next_event_time, tree[parent, 3], new_ID, -1, next_event_time)
          tree[parent, 5] <- next_event_time  # update last daughter event
          tree <- rbind(tree, new_spec)
          N <- N + 1
        } else {
          # extinction
          # pick random species
          to_remove <- sample(which(tree[, 4] == -1), 1)
          P <- P - (t - tree[to_remove, 5]) # we have to correct P for the extinct species 
          tree[to_remove, 4] <- next_event_time
          N <- N - 1
        }
      }
    }
    t = next_event_time
  }
  
  tree[, 1] <- max_t - tree[, 1]
  
  for_ddd <- as.matrix(tree[, 1:4])
  
  phy <- DDD::L2phylo(as.matrix(tree[, 1:4]), dropextinct = FALSE)
  
  t <- min(t, max_t)
  
  return(list("phy" = phy, 
              "result" = c(N, t, P)))
}


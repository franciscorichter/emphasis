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
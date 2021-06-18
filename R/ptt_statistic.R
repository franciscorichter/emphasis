#' function to calculate the ptt statistic for a tree
#' @description this function calculates the ptt statistic for a tree
#' @param tree phylogenetic tree
#' @param dt time step of the discrete integration
#' @export
calculate_ptt <- function(tree, dt = 1e-4) {
  if (!ape::is.binary(tree)) {
    stop("phylogeny must be binary")
  }
  
  if (any(ape::branching.times(tree) < 0.0)) {
    stop("tree cannot have negative branching times")
  }
  
  if (!inherits(tree, "phylo")) {
    stop("phylogeny must be of type phylo")
  }
  
  # we have extinct lineages
  ltab <- DDD::phylo2L(tree)
  crown_age <- ltab[1, 1]
  ltab[, 1] <- crown_age - ltab[, 1] 
  
  ltab[which(ltab[, 4] == -1), 4 ] <- 0
  ltab[, 4]  <- crown_age - ltab[, 4]
  
  t <- seq(0, crown_age + dt, by = dt)
  Pt <- rep(0, length(t))
  Nt <- rep(0, length(t))
  for (i in 1:length(ltab[, 1])) {
    left <- ltab[i, 1]
    right <- ltab[i, 4]
    left_index <- min(which(t >= left))
    right_index <- min(which(t >= right))
    Nt[left_index:right_index] <- Nt[left_index:right_index] + 1
    
    to_add <- seq(from = 0, 
                  to = (right - left), 
                  length.out = length(left_index:right_index))
    Pt[left_index:right_index] <- Pt[left_index:right_index] + to_add
  }
  Ptt <- (Pt - t) / Nt
  
  return(Ptt)
}

sim_tree_pd <- function(pars, max_t) {
  N <- 2

  tree <- c(0, 0) # we only need the birth dates
  t <- 0
  mu <- pars[1]
  P <- 0
  while (t < max_t && N > 0) {
    spec_rate <- pars[2] + 
                 pars[3]*N  +  
                 ((P + N * (max_t - t) - t) / N) * pars[4]
    
    total_rate <- (spec_rate + mu) * N    
    u1 <- stats::runif(1, 0, 1)
    next_event_time = max_t - log(x = u1) / total_rate
    P = P + N * (next_event_time - max_t)
    if (next_event_time < max_t) {
      # determine type of event:
      focal_spec <- pars[2] + 
                    pars[3]*N  +  
                    ((P + N * (next_event_time - t) - t) / N) * pars[4]
      
      if (runif(1) < focal_spec / (focal_spec + mu)) {
        N <- N + 1
        tree <- c(tree, next_event_time)
      } else {
        N <- N - 1
        # pick random species
        to_remove <- sample(1:length(tree), 1)
        P <- P - (t - tree[to_remove])
        tree <- tree[-to_remove]
      }
    } else {
      # reached max time
      break
    }
    t = next_event_time
  }
     
  return(c(N, t))
}



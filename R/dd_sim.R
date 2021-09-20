#' simulation function to simulate a minimal dd model
#' @description super fast function to simulate the process of diversification
#' with diversity dependence
#' @param pars parameter vector with c(mu, lambda_0, beta_N)
#' @param max_t crown age
#' @return a vector with population size and time of extinction
#' (equal to the crown age in the absence of extinction)
#' @export
sim_tree_dd_minimal <- function(pars, max_t) {
  N <- 2
  mu <- max(0, pars[1])
  t <- 0
  while (t < max_t && N >= 2 && N < 1000) {

    lambda_ct <- max(0, pars[2] +
                        pars[3] * N)
    rate_max <- (lambda_ct + mu) * N
    next_t <- t + stats::rexp(n = 1, rate = rate_max)
    
    if (next_t < max_t) {
      if (stats::runif(1, 0, 1) < lambda_ct / (lambda_ct + mu)) {
        # speciation
        N <- N + 1
      } else {
        # extinction
        N <- N - 1
      }
    }
    t <- next_t
  }
  t <- min(t, max_t)
  return(c(N, t))
}

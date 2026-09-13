## H21 part E: does the range-scaled rule fire before the EM drift is finished?
## BDI cr fits from a FAR start (1.2, 0.9) with box [0,3]^2 vs [0,30]^2 vs [0,300]^2.
## Same tree as H21.R. Report iterations, stop, pars, distance to bd_ML, loglik deficit.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(TreeSim); library(DDD); library(ape)})
lam <- 0.5; mu <- 0.2
set.seed(21)
t <- TreeSim::sim.bd.taxa(25, 1, lam, mu, complete = FALSE)[[1]]
brts <- sort(as.numeric(ape::branching.times(t)), decreasing = TRUE)
ll <- function(p) DDD::bd_loglik(c(p[1], p[2], 0, 0), c(0, 0, 1, 0, 2), brts, 0)
mle <- c(0.435130, 0.225782); llmax <- ll(mle)   # from H21.R (DDD::bd_ML)
tol <- 1e-3
fit_once <- function(ub, init, max_time = 60) {
  t0 <- proc.time()[3]
  f <- estimate_rates(brts, model = "cr", method = "mcem", init_pars = init,
        control = list(lower_bound = c(0,0), upper_bound = c(ub,ub), sampling = "bdi",
                       num_threads = 1L, max_time = max_time, verbose = FALSE))
  m <- f$details$mcem
  data.frame(ub = ub, stop = f$details$stop_reason, iters = f$details$iterations,
             lambda = f$pars[1], mu = f$pars[2],
             d_lambda = f$pars[1] - mle[1], d_mu = f$pars[2] - mle[2],
             ll_deficit = llmax - ll(f$pars), last_delta = tail(m$delta_max, 1),
             last_step_abs = max(abs(unlist(m[nrow(m), c("par1","par5")]) - unlist(m[nrow(m)-1, c("par1","par5")]))),
             secs = proc.time()[3] - t0, row.names = NULL)
}
init <- c(1.2, 0.9)
cat(sprintf("start=(%.2f,%.2f)  bd_ML=(%.4f,%.4f)  loglik(start)=%.3f  llmax=%.3f\n", init[1], init[2], mle[1], mle[2], ll(init), llmax))
for (ub in c(3, 30, 300)) {
  res <- do.call(rbind, lapply(1:4, function(i) fit_once(ub, init)))
  cat(sprintf("\n--- box [0,%g]^2 : tol*range = |dtheta| < %.3g ---\n", ub, tol * ub))
  print(res, digits = 4, row.names = FALSE)
}
## one verbose trajectory for ub=300 to show the drift at the moment "converged" fires
cat("\n--- trajectory, ub=300 ---\n")
f <- estimate_rates(brts, model = "cr", method = "mcem", init_pars = init,
      control = list(lower_bound = c(0,0), upper_bound = c(300,300), sampling = "bdi",
                     num_threads = 1L, max_time = 60, verbose = FALSE))
m <- f$details$mcem; m$abs_step <- c(NA, pmax(abs(diff(m$par1)), abs(diff(m$par5))))
print(m[, c("par1","par5","fhat","delta_max","abs_step")], digits = 4)
cat("stop:", f$details$stop_reason, "\n")

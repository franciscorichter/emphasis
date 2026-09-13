## H21 probe: time one default MCEM cr fit (BDI and thinning) on a 25-tip tree.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(TreeSim); library(DDD); library(ape)})
lam <- 0.5; mu <- 0.2
set.seed(21)
t <- TreeSim::sim.bd.taxa(25, 1, lam, mu, complete = FALSE)[[1]]
brts <- sort(as.numeric(ape::branching.times(t)), decreasing = TRUE)
cat("n_tips", length(brts)+1, "crown", brts[1], "\n")
lb <- c(0, 0); ub <- c(3, 3)
for (samp in c("bdi", "dynamic_fresh")) {
  t0 <- proc.time()[3]
  f <- estimate_rates(brts, model = "cr", method = "mcem",
        control = list(lower_bound = lb, upper_bound = ub, sampling = samp,
                       num_threads = 1L, max_time = 120, verbose = FALSE))
  cat(sprintf("%s: %.1fs  iters=%d  stop=%s  pars=%s  loglik=%.3f\n", samp, proc.time()[3]-t0,
      f$details$iterations, f$details$stop_reason, paste(round(f$pars,4), collapse=","), f$loglik))
  print(utils::head(f$details$mcem[, c("par1","par5","fhat","delta_max","time")], 12))
}

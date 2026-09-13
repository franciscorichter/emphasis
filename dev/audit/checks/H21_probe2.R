.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(TreeSim); library(DDD); library(ape)})
lam <- 0.5; mu <- 0.2
set.seed(21)
t <- TreeSim::sim.bd.taxa(25, 1, lam, mu, complete = FALSE)[[1]]
brts <- sort(as.numeric(ape::branching.times(t)), decreasing = TRUE)
for (p in list(c(1.5,1.5), c(0.5,0.5), c(0.35,0.08))) {
  t0 <- proc.time()[3]
  r <- tryCatch(emphasis:::.augment_tree_bdi(brts, p, c(0L,0L,0L), 200L, 1e4, 0L, 1.0), error = function(e) e)
  if (inherits(r, "error")) cat("pars", p, "ERROR:", conditionMessage(r), "\n") else
    cat("pars", p, "ntrees", length(r$trees), "fhat", r$fhat, "time", proc.time()[3]-t0, "\n")
}
t0 <- proc.time()[3]
f <- estimate_rates(brts, model = "cr", method = "mcem",
      control = list(lower_bound = c(0,0), upper_bound = c(3,3), sampling = "bdi",
                     num_threads = 1L, max_time = 60, verbose = TRUE))
cat("stop:", f$details$stop_reason, " iters:", f$details$iterations, " t:", proc.time()[3]-t0, "\n")
t0 <- proc.time()[3]
f <- estimate_rates(brts, model = "cr", method = "mcem",
      control = list(lower_bound = c(0,0), upper_bound = c(1,1), sampling = "bdi",
                     num_threads = 1L, max_time = 100, verbose = FALSE))
cat("ub=1: stop:", f$details$stop_reason, " iters:", f$details$iterations, " t:", proc.time()[3]-t0, " pars:", f$pars, "\n")
print(utils::tail(f$details$mcem[, c("par1","par5","fhat","delta_max","time")], 8))

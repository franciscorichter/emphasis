.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(TreeSim); library(ape)})
set.seed(7)
t <- TreeSim::sim.bd.taxa(20, 1, 0.6, 0.3, complete = FALSE)[[1]]
brts <- sort(as.numeric(ape::branching.times(t)), decreasing = TRUE)
f4 <- estimate_rates(brts, model = "cr", method = "mcem",
        control = list(lower_bound = c(0,0), upper_bound = c(2,2), sampling = "bdi",
                       num_threads = 1L, max_time = 30, verbose = FALSE))
cat("stop:", f4$details$stop_reason, " iterations:", if (is.null(f4$details$iterations)) "NULL" else f4$details$iterations,
    " pars:", f4$pars, " loglik:", f4$loglik, "\n")
# asymmetric box: midpoint has lambda != mu
f5 <- estimate_rates(brts, model = "cr", method = "mcem",
        control = list(lower_bound = c(0,0), upper_bound = c(2,1), sampling = "bdi",
                       num_threads = 1L, max_time = 20, verbose = FALSE))
cat("asym box: stop:", f5$details$stop_reason, " iterations:", f5$details$iterations, " pars:", f5$pars, "\n")

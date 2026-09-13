## H12 (part 3) — symmetric default box: lower=c(0,0), upper=c(1,1) gives
## init = centre = (0.5, 0.5), i.e. lam == mu exactly.  What happens?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({library(emphasis); library(ape)})
set.seed(121)
tr <- ape::rphylo(15, 0.5, 0.2); tr$edge.length <- tr$edge.length / max(ape::branching.times(tr)) * 5
ns <- asNamespace("emphasis")
cat(".bdi_p_cr(1, 0.5, 0.5, 5) =", ns$.bdi_p_cr(1, 0.5, 0.5, 5), "\n")
msgs <- character(0)
fit <- withCallingHandlers(
  estimate_rates(tr, model = "cr", method = "mcem",
                 control = list(lower_bound = c(0, 0), upper_bound = c(1, 1),
                                sample_size = 20L, max_iter = 20L, num_threads = 1L, verbose = TRUE)),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
cat("pars:", fit$pars, " loglik:", fit$loglik, " stop:", fit$details$stop_reason, "\n")
cat("E-step failure messages:", sum(grepl("E-step failed", msgs)), "\n")
cat(tail(msgs, 3), sep = "\n")
e <- tryCatch(ns$.augment_tree_bdi(tr, pars = c(0.5, 0.5), sample_size = 2L), error = function(e) conditionMessage(e))
cat("direct .augment_tree_bdi at (0.5,0.5):", format(e), "\n")
## same box, thinning sampler
fit_t <- estimate_rates(tr, model = "cr", method = "mcem",
                 control = list(lower_bound = c(0, 0), upper_bound = c(1, 1), sampling = "dynamic_fresh",
                                sample_size = 20L, max_iter = 20L, max_time = 60, num_threads = 1L))
cat("thinning from (0.5,0.5): pars:", fit_t$pars, " stop:", fit_t$details$stop_reason, " iters:", fit_t$details$iterations, "\n")

## H12: can MCEM started at mu>lam PROCEED (with corrupt E-steps) on a
## low-rate tree instead of stopping with e_step_failure?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({library(emphasis); library(ape); library(DDD)})
set.seed(3); tr <- ape::rphylo(8, 0.5, 0.2)
tr$edge.length <- tr$edge.length / max(ape::branching.times(tr)) * 5
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
ml <- suppressWarnings(DDD::bd_ML(brts, cond = 0, btorph = 0, soc = 2, verbose = FALSE, initparsopt = c(0.3, 0.1)))
cat("bd_ML(cond=0):", ml$lambda0, ml$mu0, "\n")
for (r in 1:3) {
  msgs <- character(0)
  fit <- withCallingHandlers(suppressWarnings(
    estimate_rates(tr, model = "cr", method = "mcem",
                   control = list(lower_bound = c(0, 0), upper_bound = c(0.2, 0.4),  # init = (0.1, 0.2)
                                  sample_size = 20L, max_iter = 15L, max_time = 60, num_threads = 1L, verbose = TRUE))),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  cat(sprintf("rep %d: pars (%.4f, %.4f) loglik %.3f stop %s iters %s  E-step-failed msgs %d\n", r,
              fit$pars[1], fit$pars[2], fit$loglik, fit$details$stop_reason,
              format(fit$details$iterations), sum(grepl("E-step failed", msgs))))
  if (!is.null(fit$details$mcem)) print(round(head(fit$details$mcem[, 1:2], 6), 4))
}

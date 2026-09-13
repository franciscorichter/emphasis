source("/Users/pancho/Code/emphasis/dev/audit/review/1.6_common.R")
set.seed(7)
fit <- estimate_rates(brts22, model="cr", method="mcem", init_pars=c(1,0.3), control=list(sampling="dynamic_fresh", num_trees=20L, max_iter=3L, lower_bound=c(0,0), upper_bound=c(4,4), num_threads=1L, max_missing = 3))
dg <- tryCatch(diagnose_mcem(fit), error = function(e) e)
if (inherits(dg, "error")) cat("diagnose_mcem ERROR:", conditionMessage(dg), "\n") else {
  cat("diagnose_mcem convergence table rows:", nrow(dg$convergence), " fit$iterations:", fit$iterations, "\n"); print(dg$convergence)
  cat("trace rejected_overruns:", fit$details$mcem$rejected_overruns, "\n")
}
# default-control behaviour: how many iterations does a default thinning fit take now?
t0 <- proc.time()[3]
fit <- estimate_rates(brts22, model="cr", method="mcem", init_pars=c(1,0.3), control=list(sampling="dynamic_fresh", lower_bound=c(0,0), upper_bound=c(4,4), num_threads=1L, max_time=120))
cat(sprintf("POST default control thinning: stop=%s iterations=%d pars=(%s) loglik=%.3f %.0fs; delta tail %s\n", fit$stop_reason, fit$iterations, paste(signif(fit$pars,3), collapse=","), fit$loglik, proc.time()[3]-t0, paste(signif(tail(fit$details$mcem$delta_max, 6),2), collapse=" ")))

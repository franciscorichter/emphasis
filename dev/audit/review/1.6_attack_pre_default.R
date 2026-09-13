.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
brts22 <- c(2.984340, 2.456040, 2.041011, 1.919104, 1.703097, 1.313093, 0.939284, 0.715800, 0.553944, 0.455806, 0.427175, 0.308270, 0.174959, 0.141573, 0.112195, 0.069056, 0.058366, 0.044335, 0.035010)
t0 <- proc.time()[3]
fit <- estimate_rates(brts22, model="cr", method="mcem", init_pars=c(1,0.3), control=list(sampling="dynamic_fresh", lower_bound=c(0,0), upper_bound=c(4,4), num_threads=1L, max_time=120))
cat(sprintf("PRE default control thinning: stop=%s iterations=%d pars=(%s) loglik=%.3f %.0fs\n", fit$details$stop_reason, nrow(fit$details$mcem), paste(signif(fit$pars,3), collapse=","), fit$loglik, proc.time()[3]-t0))

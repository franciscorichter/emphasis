.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)
brts20 <- c(5,4.4,3.9,3.6,3.1,2.8,2.5,2.2,2.0,1.7,1.5,1.3,1.1,0.9,0.75,0.6,0.45,0.3,0.15)
cat("== C1 thinning run: iterations vs trace\n")
set.seed(21)
f <- suppressWarnings(estimate_rates(brts20, model="cr", method="mcem", init_pars=c(1,0.3),
  control=list(lower_bound=c(0,0), upper_bound=c(4,4), sampling="dynamic_fresh",
               sample_size=20L, max_iter=4L, max_time=60, num_threads=1L)))
cat("iterations:", f$iterations, " nrow(trace):", nrow(f$details$mcem),
    " sum(m_step):", sum(f$details$mcem$m_step), " stop:", f$stop_reason, "\n")
print(names(f$details$mcem))

cat("\n== C2 BDI all-fail (lambda fixed at 0)\n")
set.seed(22)
f2 <- suppressWarnings(estimate_rates(brts20, model="cr", method="mcem",
  control=list(lower_bound=c(0,0.1), upper_bound=c(0,0.1), sampling="bdi",
               sample_size=10L, max_iter=20L, max_time=60, num_threads=1L)))
cat("loglik:", f2$loglik, " stop:", f2$stop_reason, " iters:", f2$iterations,
    " n_pars:", f2$n_pars, " AIC:", f2$AIC, "\n")
print(f2)

cat("\n== C3 sample_size = 1 (loglik_var NA) + maxN boundary maxN == sample_size\n")
set.seed(23)
f3 <- tryCatch(suppressWarnings(estimate_rates(brts20, model="cr", method="mcem", init_pars=c(1,0.3),
  control=list(lower_bound=c(0,0), upper_bound=c(4,4), sampling="dynamic_fresh",
               num_trees=50L, maxN=50L, max_iter=1L, max_time=60, num_threads=1L))),
  error=function(e) conditionMessage(e))
if (is.character(f3)) cat("ERROR:", f3, "\n") else
  cat("loglik:", f3$loglik, "stop:", f3$stop_reason, "iters:", f3$iterations, "\n")

cat("\n== C4 num_trees=1\n")
set.seed(24)
f4 <- tryCatch(suppressWarnings(estimate_rates(brts20, model="cr", method="mcem", init_pars=c(1,0.3),
  control=list(lower_bound=c(0,0), upper_bound=c(4,4), sampling="bdi",
               num_trees=1L, max_iter=2L, num_threads=1L))),
  error=function(e) conditionMessage(e))
if (is.character(f4)) cat("ERROR:", f4, "\n") else {
  cat("loglik:", f4$loglik, "var:", f4$loglik_var, "stop:", f4$stop_reason, "iters:", f4$iterations, "\n")
  print(f4) }

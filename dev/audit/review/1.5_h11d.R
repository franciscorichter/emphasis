lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
brts_dd <- c(6, 4.848493, 4.401821, 3.108164, 3.073914, 2.835023, 1.838828, 0.50463)
lb <- c(0.01, -1, 0.001, 0); ub <- c(5, 0, 2, 0)
set.seed(11)
fit <- estimate_rates(brts_dd, method="mcem", model="dd", init_pars=c(1.5,-0.12,0.4,0),
  control=list(lower_bound=lb, upper_bound=ub, sample_size=200L, max_iter=4L,
               max_missing=30L, num_threads=1L, sampling="bdi"))
cat("stop_reason:", fit$details$stop_reason, " iters:", format(fit$details$iterations), "\n")
m <- fit$details$mcem
print(m)
cat("loglik:", fit$loglik, "\n")

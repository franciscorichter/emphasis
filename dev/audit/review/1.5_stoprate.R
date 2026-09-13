.libPaths(c(commandArgs(TRUE)[1], .libPaths())); library(emphasis)
brts12 <- c(6.385824,2.063997,1.19255,0.923743,0.884126,0.822976,0.725585,0.718801,0.539214,0.077287,0.068123)
set.seed(99)
f <- estimate_rates(brts12, model="cr", method="mcem", init_pars=c(1.2,0.9),
  control=list(lower_bound=c(0,0), upper_bound=c(3,3), sampling="bdi",
               sample_size=200L, num_threads=1L, max_iter=60L, tol=1e-2, patience=3L))
m <- f$details$mcem; em <- m[m$m_step,]
cat("stop:", f$details$stop_reason, " iters:", f$details$iterations, "\n")
cat("delta_max:", paste(signif(em$delta_max,2), collapse=" "), "\n")
cat("frac(delta<1e-2):", mean(em$delta_max < 1e-2), "  median:", signif(median(em$delta_max),3), "\n")
cat("pars:", round(unname(f$pars),4), "\n")

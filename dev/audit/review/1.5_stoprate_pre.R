.libPaths(c(commandArgs(TRUE)[1], .libPaths())); library(emphasis)
brts12 <- c(6.385824,2.063997,1.19255,0.923743,0.884126,0.822976,0.725585,0.718801,0.539214,0.077287,0.068123)
for (tl in c(1e-3, 1e-2)) { set.seed(99)
f <- estimate_rates(brts12, model="cr", method="mcem", init_pars=c(1.2,0.9),
  control=list(lower_bound=c(0,0), upper_bound=c(3,3), sampling="bdi",
               sample_size=200L, num_threads=1L, max_iter=60L, tol=tl, patience=3L))
m <- f$details$mcem
cat("tol",tl,"stop:", f$details$stop_reason, " nrow(mcem):", nrow(m),
    " pars:", round(unname(f$pars),4), " loglik:", round(f$loglik,3), "\n")
cat("  delta:", paste(signif(m$delta_max,2), collapse=" "), "\n") }

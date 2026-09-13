lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths()))
suppressMessages(library(emphasis))
set.seed(42)
brts11 <- c(5,4.745201,4.530461,4.067871,3.622029,2.929002,1.468698,1.386875,1.302139,0.044729)
cat("build:", basename(lib), "\n")
ab <- tryCatch(auto_bounds(brts11, model="cr", link="exponential"), error=function(e) conditionMessage(e))
cat("auto_bounds cr/exponential:\n"); print(ab)
ab0 <- tryCatch(auto_bounds(brts11, model="cr", link="linear"), error=function(e) conditionMessage(e))
cat("auto_bounds cr/linear:\n"); print(ab0)

t0 <- proc.time()[3]
f <- tryCatch(estimate_rates(brts11, model="cr", link="linear", method="mcem", sampling="bdi",
            lower_bound=c(0,0), upper_bound=c(1,1), num_threads=1,
            control=list(sample_size=20, max_iter=6, patience=2, verbose=FALSE)),
       error=function(e) paste("ERROR:", conditionMessage(e)))
cat("\nsymmetric box [0,1]^2 (init = (0.5,0.5), lam == mu):\n")
if (is.character(f)) cat(f,"\n") else
  cat(sprintf("  pars=(%s) loglik=%s stop=%s iters=%s  elapsed=%.1fs\n",
      paste(signif(f$pars,5),collapse=","), format(f$loglik),
      paste(f$stop_reason, collapse=""), paste(f$iterations,collapse=""), proc.time()[3]-t0))

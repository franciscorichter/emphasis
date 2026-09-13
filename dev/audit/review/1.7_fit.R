lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
b12 <- c(10,8.3,7.1,6.2,5.5,4.4,3.9,3.1,2.2,1.5,0.9,0.3)
go <- function(samp, nt) {
  t0 <- proc.time()[3]
  f <- try(estimate_rates(tree = b12, method = "mcem", model = "cr", init_pars = c(0.5,0.1),
        control = list(num_trees = 30L, max_iter = 5L, sampling = samp, num_threads = nt,
                       lower_bound = c(0.01,0.001), upper_bound = c(5,5))), silent=TRUE)
  if (inherits(f,"try-error")) { cat(sprintf("%-14s nt=%d ERROR %s\n", samp, nt, conditionMessage(attr(f,"condition")))); return(invisible()) }
  d <- f$details
  cat(sprintf("%-14s nt=%d pars=(%.3f,%.3f) loglik=%.3f var=%s iters=%s stop=%s ess=%s  %.1fs\n",
    samp, nt, f$pars[1], f$pars[2], f$loglik,
    format(f$loglik_var, digits=3),
    format(if (!is.null(f$iterations)) f$iterations else d$iterations),
    format(if (!is.null(f$stop_reason)) f$stop_reason else d$stop_reason),
    format(if (!is.null(d$final_IS$ess)) d$final_IS$ess else NA, digits=4),
    proc.time()[3]-t0))
}
go("dynamic_fresh", 1L); go("dynamic_fresh", 8L); go("bdi", 1L)

lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths()))
suppressMessages(library(emphasis))
brts11 <- c(5,4.745201,4.530461,4.067871,3.622029,2.929002,1.468698,1.386875,1.302139,0.044729)
cat("build:", basename(lib), "\n")
go <- function(lab, init, link="linear", lb=c(0,0), ub=c(1,1)) {
  set.seed(7); t0 <- proc.time()[3]
  f <- tryCatch(estimate_rates(brts11, method="mcem", model="cr", link=link, init_pars=init,
        control=list(sampling="bdi", sample_size=20, max_iter=5, patience=2,
                     num_threads=1, lower_bound=lb, upper_bound=ub, verbose=FALSE)),
        error=function(e) paste("ERROR:", conditionMessage(e)),
        warning=function(w) paste("WARN:", conditionMessage(w)))
  if (is.character(f)) { cat(sprintf("%-28s %s\n", lab, f)); return(invisible()) }
  cat(sprintf("%-28s pars=(%s) loglik=%s stop=%s iters=%s [%.0fs]\n", lab,
     paste(signif(f$pars,5),collapse=","), format(signif(f$loglik,8)),
     paste(c(f$stop_reason, f$details$stop_reason),collapse="/"),
     paste(c(f$iterations, f$details$iterations),collapse="/"), proc.time()[3]-t0))
}
go("linear init lam==mu (.5,.5)", c(0.5,0.5))
go("linear init mu>lam (.3,.8)",  c(0.3,0.8))
go("linear init mu<lam (.8,.3)",  c(0.8,0.3))
go("exp link init (0,5) mu>>lam", c(0,5), link="exponential", lb=c(-5,-5), ub=c(6,6))

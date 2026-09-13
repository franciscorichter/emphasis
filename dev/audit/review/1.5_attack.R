.libPaths(c(commandArgs(TRUE)[1], .libPaths())); library(emphasis)
brts12 <- c(6.385824,2.063997,1.19255,0.923743,0.884126,0.822976,0.725585,0.718801,0.539214,0.077287,0.068123)
brts3  <- c(4.0, 1.2)
go <- function(tag, ...) {
  r <- tryCatch(estimate_rates(...), error=function(e) paste("ERROR:", conditionMessage(e)),
                warning=function(w) paste("WARN:", conditionMessage(w)))
  if (is.character(r)) { cat(sprintf("%-28s %s\n", tag, substr(r,1,110))); return(invisible()) }
  cat(sprintf("%-28s stop=%-14s it=%-3s ll=%-10s var=%-9s pars=%s ESS=%s\n", tag,
      r$stop_reason, r$iterations, signif(r$loglik,6), signif(r$loglik_var,3),
      paste(round(unname(r$pars),3), collapse=","),
      signif(r$details$final_IS$ESS,4)))
}
ctl <- function(...) utils::modifyList(list(lower_bound=c(0,0), upper_bound=c(3,3),
  sampling="bdi", sample_size=20L, num_threads=1L, max_iter=5L), list(...))
set.seed(1); go("cr linear",      brts12, model="cr", method="mcem", init_pars=c(1.2,0.9), control=ctl())
set.seed(1); go("cr sample_size=1", brts12, model="cr", method="mcem", init_pars=c(1.2,0.9), control=ctl(sample_size=1L))
set.seed(1); go("cr mu>lambda init", brts12, model="cr", method="mcem", init_pars=c(0.5,1.2), control=ctl())
set.seed(1); go("cr mu==lambda init", brts12, model="cr", method="mcem", init_pars=c(0.9,0.9), control=ctl())
set.seed(1); go("cr default init symbox", brts12, model="cr", method="mcem", control=ctl())
set.seed(1); go("cr exp link", brts12, model="cr", method="mcem", link="exponential",
                init_pars=c(log(1.2),log(0.9)), control=ctl(lower_bound=c(-5,-5), upper_bound=c(2,2)))
set.seed(1); go("cr 3-tip tree", brts3, model="cr", method="mcem", init_pars=c(1.2,0.9), control=ctl())
set.seed(1); go("cr rates near 0", brts12, model="cr", method="mcem", init_pars=c(0.02,0.005), control=ctl())
set.seed(1); go("cr maxN<num_trees", brts12, model="cr", method="mcem", init_pars=c(1.2,0.9), control=ctl(maxN=10L))

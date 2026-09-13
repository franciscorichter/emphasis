lib <- commandArgs(TRUE)[1]; .libPaths(c(lib,.libPaths())); library(emphasis)
set.seed(11)
sim <- DDD::dd_sim(pars = c(1.5, 0.4, 1.1/0.12), age = 6, ddmodel = 1)
brts <- as.numeric(ape::branching.times(sim$tes))
pars <- c(1.5, -0.12, 0.4, 0)
lb <- c(0.01, -1, 0.001, 0); ub <- c(5, 0, 3, 0)
t0 <- Sys.time()
fit <- try(estimate_rates(brts, method="mcem", model="dd", init_pars=pars,
  control=list(lower_bound=lb, upper_bound=ub, sample_size=100L, max_iter=4L,
               max_missing=30L, num_threads=1L, sampling="bdi", verbose=FALSE)), silent=TRUE)
if (inherits(fit,"try-error")) cat("ERROR:", conditionMessage(attr(fit,"condition")), "\n") else {
  cat("pars:", paste(round(as.numeric(fit$pars),4), collapse=" "), "\n")
  cat("loglik:", format(fit$loglik), "\n")
  cat("stop_reason:", format(fit$details$stop_reason), " iterations:", format(fit$details$iterations), "\n")
  fi <- fit$details$final_IS
  if (!is.null(fi)) cat(sprintf("final_IS fhat=%s ESS=%s n+Inf=%d n-Inf=%d nNaN=%d\n",
    format(fi$fhat), format(fi$ESS), sum(fi$logf==Inf,na.rm=TRUE), sum(fi$logf==-Inf,na.rm=TRUE), sum(is.nan(fi$logf))))
  cat("moved:", !isTRUE(all.equal(as.numeric(fit$pars), pars)), "\n")
}
cat("secs:", round(as.numeric(difftime(Sys.time(),t0,units="secs")),1), "\n")

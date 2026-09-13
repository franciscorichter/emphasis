lib <- Sys.getenv("EMPH_LIB"); .libPaths(c(lib, .libPaths()))
suppressMessages(library(emphasis)); ns <- asNamespace("emphasis")
brts1 <- c(5.000000, 3.822147, 3.059192, 1.966150, 1.604249, 1.358484,
           0.622046, 0.467531, 0.393533, 0.285495, 0.276148, 0.225482,
           0.201295, 0.154060, 0.057677, 0.044974, 0.041380, 0.038051, 0.013895)
cat("build:", lib, "\n")
ctrl <- list(sampling="dynamic_fresh", num_trees=30L, maxN=3000L, max_iter=25L,
             lower_bound=c(0,0), upper_bound=c(2,1), num_threads=1L)
for (s in 1:2) {
  set.seed(100+s)
  f <- estimate_rates(brts1, model="cr", method="mcem", init_pars=c(0.5,0.1), control=ctrl)
  cat(sprintf("seed %d: pars=%s loglik=%.4f se=%.4f iters=%s stop=%s rows=%d\n", 100+s,
      paste(sprintf("%.5f", f$pars), collapse=","), f$loglik, sqrt(f$loglik_var),
      if (is.null(f$details$iterations)) NA else f$details$iterations,
      f$details$stop_reason, nrow(f$details$mcem)))
}
# wide box: does the estimate still depend on the box?
for (ub in list(c(2,1), c(30,30), c(300,300))) {
  set.seed(55)
  c2 <- ctrl; c2$upper_bound <- ub; c2$max_iter <- 25L
  f <- estimate_rates(brts1, model="cr", method="mcem", init_pars=c(0.5,0.1), control=c2)
  cat(sprintf("box ub=%-10s pars=%s stop=%s iters=%s\n", paste(ub, collapse=","),
      paste(sprintf("%.5f", f$pars), collapse=","), f$details$stop_reason,
      if (is.null(f$details$iterations)) nrow(f$details$mcem) else f$details$iterations))
}
# all parameters fixed: does it report "converged" with zero movement?
set.seed(9)
c3 <- ctrl; c3$lower_bound <- c(0.5,0.1); c3$upper_bound <- c(0.5,0.1); c3$max_iter <- 6L
f <- tryCatch(estimate_rates(brts1, model="cr", method="mcem", init_pars=c(0.5,0.1), control=c3),
              error=function(e) conditionMessage(e))
if (is.list(f)) cat(sprintf("fixed box: stop=%s iters=%s delta=%s\n", f$details$stop_reason,
    if (is.null(f$details$iterations)) NA else f$details$iterations,
    paste(signif(f$details$mcem$delta_max,3), collapse=","))) else cat("fixed box ERROR:", f, "\n")

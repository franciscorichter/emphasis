lib <- Sys.getenv("EMPH_LIB"); .libPaths(c(lib, .libPaths()))
suppressMessages(library(emphasis)); ns <- asNamespace("emphasis")
ex <- function(p, m) ns$.expand_pars(p, m)
brts1 <- c(5.000000, 3.822147, 3.059192, 1.966150, 1.604249, 1.358484,
           0.622046, 0.467531, 0.393533, 0.285495, 0.276148, 0.225482,
           0.201295, 0.154060, 0.057677, 0.044974, 0.041380, 0.038051, 0.013895)
cr <- c(0L,0L,0L); dd <- c(1L,0L,0L); dmod <- c(0L,0L,1L)
run <- function(tag, pars, model, link, lb, ub, ss=20L, maxN=2000L, it=4L,
                max_missing=1e4, seed=7, nt=1L) {
  set.seed(seed)
  r <- tryCatch(ns$.mcem_dynamic_fresh(brts1, ex(pars,model), sample_size=ss, maxN=maxN,
        max_missing=max_missing, lower_bound=ex(lb,model), upper_bound=ex(ub,model),
        max_iter=it, xtol=1e-3, tol=1e-2, patience=3L, num_threads=nt, verbose=FALSE,
        model=model, link=link, max_time=60), error=function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) { cat(sprintf("%-24s %s\n", tag, r)); return(invisible(NULL)) }
  nr <- if (is.null(r$mcem)) 0L else nrow(r$mcem)
  cat(sprintf("%-24s stop=%-12s it=%s nrow=%s fin=%s nf=%s maxN=%s pars=%s\n", tag,
      r$stop_reason, r$iterations, nr, r$final_estep, r$n_failed, r$maxN,
      paste(sprintf("%.4g", ns$.contract_pars(r$pars, model)), collapse=",")))
  invisible(r)
}
cat("---- dd / d models (4 compact pars) ----\n")
run("dd/linear", c(0.6,-0.01,0.1,0), dd, 0L, c(0,-1,0,0), c(2,0,1,0))
run("dd/exponential", c(-0.5,-0.01,-2.3,0), dd, 1L, c(-5,-1,-5,0), c(2,0,2,0))
run("dd/gaussian", c(0.6,-0.01,0.1,0), dd, 2L, c(0,-1,0,0), c(2,0,1,0))
run("d/linear", c(0.6,0.01,0.1,0), dmod, 0L, c(0,-1,0,0), c(2,1,1,0))
cat("---- threads ----\n")
run("cr nt=2", c(0.5,0.1), cr, 0L, c(0,0), c(2,1), nt=2L)
cat("---- rejections: max_missing = 2 ----\n")
r <- run("cr max_missing=2", c(0.6,0.5), cr, 0L, c(0,0), c(2,1), ss=20L, maxN=5000L,
         it=2L, max_missing=2, seed=4)
if (!is.null(r)) {
  print(r$mcem[, c("fhat","rejected","rejected_overruns","rejected_lambda",
                   "rejected_zero_weights","n_rejected","final_estep")])
  cat("final_IS n_rejected:", r$final_IS$n_rejected,
      " zero_w:", r$final_IS$rejected_zero_weights,
      " fhat:", r$final_IS$fhat, " last row fhat:", r$mcem$fhat[nrow(r$mcem)], "\n")
  cat("identical fhat?", isTRUE(all.equal(r$final_IS$fhat, r$mcem$fhat[nrow(r$mcem)])), "\n")
}
cat("---- mcem_diagnostics iteration count ----\n")
set.seed(21)
fit <- estimate_rates(brts1, model="cr", method="mcem", init_pars=c(0.5,0.1),
  control=list(sampling="dynamic_fresh", num_trees=30L, maxN=3000L, max_iter=3L,
               tol=1e-2, patience=3L, lower_bound=c(0,0), upper_bound=c(2,1), num_threads=1L))
cat("fit$iterations:", fit$iterations, " nrow(details$mcem):", nrow(fit$details$mcem), "\n")
dg <- tryCatch(diagnose_mcem(fit, plot=FALSE), error=function(e) tryCatch(mcem_diagnostics(fit, plot=FALSE), error=function(e2) conditionMessage(e2)))
if (is.list(dg)) {
  cat("diag rows:", nrow(dg$convergence), " last delta:", dg$convergence$delta_max[nrow(dg$convergence)], "\n")
  print(utils::tail(dg$convergence, 2))
  cat("diag rejected col:", paste(dg$convergence$rejected, collapse=","),
      " IS n_rejected:", dg$IS_quality$n_rejected, "\n")
} else cat("diag:", dg, "\n")
print(fit)

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
suppressMessages(library(emphasis)); ns <- asNamespace("emphasis")
brts1 <- c(5,3.822147,3.059192,1.966150,1.604249,1.358484,0.622046,0.467531,
           0.393533,0.285495,0.276148,0.225482,0.201295,0.154060,0.057677,
           0.044974,0.041380,0.038051,0.013895)
cr <- c(0L,0L,0L); ex <- function(p) ns$.expand_pars(p, cr)

cat("== A. successes then permanent failure ==\n")
k <- new.env(); k$n <- 0L
cond <- function(p) { k$n <- k$n + 1L; if (k$n > 60L) stop("dead") ; 0 }
r <- suppressWarnings(ns$.mcem_dynamic_fresh(brts1, ex(c(0.5,0.1)), sample_size=5L,
   maxN=2000L, max_missing=1e4, lower_bound=ex(c(0,0)), upper_bound=ex(c(2,1)),
   max_iter=30L, xtol=1e-3, tol=1e-8, patience=3L, num_threads=1L, verbose=FALSE,
   conditional=cond, model=cr, link=0L, max_time=120))
cat("stop:", r$stop_reason, " iters:", r$iterations, " nfail:", r$n_failed,
    " final_estep:", r$final_estep, " rows:", nrow(r$mcem), "\n")
cat("pars       :", sprintf("%.6f", r$pars[c(1,5)]), "\n")
cat("last row   :", sprintf("%.6f", as.numeric(r$mcem[nrow(r$mcem), c("par1","par5")])), "\n")
cat("loglik (.run_mcem rule) =", utils::tail(r$mcem$fhat[is.finite(r$mcem$fhat)],1),
    " final_IS$fhat =", r$final_IS$fhat, "\n")
cat("pars == last row pars? ", isTRUE(all.equal(as.numeric(r$pars[c(1,5)]),
     as.numeric(r$mcem[nrow(r$mcem), c("par1","par5")]))), "\n")

cat("== B. time budget overshoot ==\n")
t0 <- proc.time()[3]
r2 <- ns$.mcem_dynamic_fresh(brts1, ex(c(0.5,0.1)), sample_size=200L, maxN=20000L,
   max_missing=1e4, lower_bound=ex(c(0,0)), upper_bound=ex(c(2,1)), max_iter=30L,
   xtol=1e-3, tol=1e-8, patience=3L, num_threads=1L, verbose=FALSE, model=cr,
   link=0L, max_time=0.001)
cat("stop:", r2$stop_reason, " iters:", r2$iterations, " rows:", nrow(r2$mcem),
    " elapsed:", round(proc.time()[3]-t0, 2), "s (budget 0.001s)\n")

cat("== C. NA/NaN robustness of rel_delta (documentation probe) ==\n")
rel <- function(new, old) max(abs(new-old)/pmax(abs(old), 1e-2))
cat("rel(c(NA,1), c(1,1)) =", rel(c(NA,1), c(1,1)), "\n")

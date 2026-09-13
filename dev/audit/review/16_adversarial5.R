.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
suppressMessages(library(emphasis)); ns <- asNamespace("emphasis")
brts1 <- c(5,3.822147,3.059192,1.966150,1.604249,1.358484,0.622046,0.467531,
           0.393533,0.285495,0.276148,0.225482,0.201295,0.154060,0.057677,
           0.044974,0.041380,0.038051,0.013895)
cr <- c(0L,0L,0L); ex <- function(p) ns$.expand_pars(p, cr)
fl <- new.env(); fl$dead <- FALSE
cond <- function(p) { if (fl$dead) stop("dead"); 0 }
r <- suppressWarnings(withCallingHandlers(
  ns$.mcem_dynamic_fresh(brts1, ex(c(0.5,0.1)), sample_size=20L, maxN=2000L,
    max_missing=1e4, lower_bound=ex(c(0,0)), upper_bound=ex(c(2,1)),
    max_iter=30L, xtol=1e-3, tol=1e-8, patience=3L, num_threads=1L, verbose=TRUE,
    conditional=cond, model=cr, link=0L, max_time=120),
  message=function(m){ if (grepl("^Iteration 3:", conditionMessage(m))) fl$dead <- TRUE
                       invokeRestart("muffleMessage") }))
cat("stop:", r$stop_reason, "iters:", r$iterations, "nfail:", r$n_failed,
    "final_estep:", r$final_estep, "rows:", nrow(r$mcem), "\n")
cat("returned pars :", sprintf("%.6f", r$pars[c(1,5)]), "\n")
cat("last trace row:", sprintf("%.6f", as.numeric(r$mcem[nrow(r$mcem), c("par1","par5")])), "\n")
cat("loglik (.run_mcem rule) =", utils::tail(r$mcem$fhat[is.finite(r$mcem$fhat)],1),
    "  final_IS$fhat =", r$final_IS$fhat, "\n")
cat("pars equals last trace row? ",
    isTRUE(all.equal(as.numeric(r$pars[c(1,5)]),
      as.numeric(r$mcem[nrow(r$mcem), c("par1","par5")]))), "\n")
cat("mcem trace:\n"); print(r$mcem[, c("par1","par5","fhat","delta_max","final_estep")])

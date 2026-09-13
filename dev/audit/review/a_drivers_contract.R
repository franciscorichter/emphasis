# (a) Seam: the two MCEM drivers -- one stopping rule, one default, one
# trace/return contract.  cr fit with each sampler on the POST-FIX build.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)
set.seed(1)

ctrl0 <- emphasis:::estimate_rates_control("mcem")
cat("## control defaults\n")
cat("tol      =", ctrl0$tol, "\n")
cat("patience =", ctrl0$patience, "\n")
cat("sampling =", ctrl0$sampling, "\n")
cat("maxN     =", if (is.null(ctrl0$maxN)) "NULL" else ctrl0$maxN, "\n")

# default eps / rel_floor used inside each driver
b <- deparse(body(emphasis:::.mcem_bdi))
d <- deparse(body(emphasis:::.mcem_dynamic_fresh))
cat("bdi eps line      :", grep("eps <- ", b, value = TRUE), "\n")
cat("thinning floor    :", grep("rel_floor <- ", d, value = TRUE), "\n")
cat("bdi tol default   :", deparse(formals(emphasis:::.mcem_bdi)$tol), "\n")
cat("thin tol default  :", deparse(formals(emphasis:::.mcem_dynamic_fresh)$tol), "\n")

# ---- a cr tree both samplers can fit --------------------------------------
set.seed(42)
brts <- sort(ape::branching.times(ape::rcoal(30)), decreasing = TRUE)
lb <- c(0, 0); ub <- c(2, 2)

fit <- function(samp, ...) {
  set.seed(7)
  estimate_rates(brts, method = "mcem", model = "cr",
                 control = list(lower_bound = lb, upper_bound = ub,
                                sampling = samp, num_trees = 50L,
                                max_iter = 12L, num_threads = 1L, ...))
}

f_bdi  <- fit("bdi")
f_thin <- fit("dynamic_fresh")

cat("\n## top-level fields\n")
show <- function(f, nm) {
  cat(nm, ": names =", paste(names(f), collapse = ","), "\n")
  cat("   pars=", paste(round(f$pars, 4), collapse = ","),
      " loglik=", f$loglik, " loglik_var=", f$loglik_var,
      " iters=", f$iterations, " stop=", f$stop_reason,
      " n_pars=", f$n_pars, " AIC=", round(f$AIC, 3), "\n", sep = "")
  d <- f$details
  cat("   details names =", paste(names(d), collapse = ","), "\n")
  cat("   trace cols    =", paste(colnames(d$mcem), collapse = ","), "\n")
  cat("   nrow(trace)   =", nrow(d$mcem), "  driver iterations =", d$iterations, "\n")
  cat("   final_IS names=", paste(names(d$final_IS), collapse = ","), "\n")
  cat("   driver loglik field present:", !is.null(d$loglik), "\n")
  cat("   tail fhat     =", utils::tail(d$mcem$fhat, 2), "\n")
  cat("   final_IS$fhat =", d$final_IS$fhat, "  fit$loglik =", f$loglik, "\n")
  invisible(NULL)
}
show(f_bdi, "bdi")
show(f_thin, "thinning")

cat("\n## contract checks\n")
req <- c("iterations", "stop_reason", "loglik_var", "final_IS", "n_failed")
for (nm in req)
  cat(sprintf("  %-12s bdi=%s thin=%s\n", nm,
              !is.null(f_bdi$details[[nm]]), !is.null(f_thin$details[[nm]])))
cat("  loglik in driver return: bdi=", !is.null(f_bdi$details$loglik),
    " thin=", !is.null(f_thin$details$loglik), "\n", sep = "")

cat("  trace cols only in bdi :",
    paste(setdiff(colnames(f_bdi$details$mcem), colnames(f_thin$details$mcem)),
          collapse = ","), "\n")
cat("  trace cols only in thin:",
    paste(setdiff(colnames(f_thin$details$mcem), colnames(f_bdi$details$mcem)),
          collapse = ","), "\n")

cat("\n## loglik == fhat at pars?\n")
cat("  bdi : last trace row fhat =", utils::tail(f_bdi$details$mcem$fhat, 1),
    " final_IS$fhat =", f_bdi$details$final_IS$fhat, "\n")
cat("  thin: last trace row fhat =", utils::tail(f_thin$details$mcem$fhat, 1),
    " final_IS$fhat =", f_thin$details$final_IS$fhat, "\n")

cat("\n## diagnose_mcem on each\n")
for (nm in c("bdi", "thin")) {
  f <- if (nm == "bdi") f_bdi else f_thin
  r <- tryCatch({ d <- emphasis::diagnose_mcem(f); "ok" },
                error = function(e) paste("ERROR:", conditionMessage(e)),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(" ", nm, "->", r, "\n")
}

cat("\n## print.emphasis_fit\n")
print(f_bdi); print(f_thin)

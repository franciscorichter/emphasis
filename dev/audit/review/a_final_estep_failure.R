# (a) Seam: what the two drivers report when the FINAL E-step (the one that is
# supposed to give loglik = fhat at the returned pars) fails, and what
# .run_mcem then reports.  The sampler is made to fail only on that last call
# by wrapping the namespace binding in this process (no files are touched).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)
ns <- asNamespace("emphasis")

set.seed(4)
brts <- sort(ape::branching.times(ape::rcoal(15)), decreasing = TRUE)
lb <- c(0, 0); ub <- c(3, 3)

patch <- function(name, n_ok) {
  orig <- get(name, envir = ns)
  i <- 0L
  new <- function(...) {
    i <<- i + 1L
    if (i > n_ok) stop("injected E-step failure (call ", i, ")")
    orig(...)
  }
  unlockBinding(name, ns); assign(name, new, envir = ns); lockBinding(name, ns)
  invisible(function() { unlockBinding(name, ns); assign(name, orig, envir = ns); lockBinding(name, ns) })
}

report <- function(tag, f) {
  d <- f$details
  cat("\n--", tag, "--\n")
  cat("  fit$loglik      =", f$loglik, "\n")
  cat("  fit$iterations  =", f$iterations, "  stop =", f$stop_reason, "\n")
  cat("  fit$AIC         =", f$AIC, "\n")
  cat("  driver 'loglik' (exact match):",
      if ("loglik" %in% names(d)) d[["loglik"]] else "FIELD ABSENT", "\n")
  cat("  details$loglik  (partial match as a user would write it):", d$loglik, "\n")
  cat("  final_IS NULL?  ", is.null(d$final_IS),
      if (!is.null(d$final_IS)) paste0("  final_IS$fhat=", d$final_IS$fhat) else "", "\n")
  cat("  trace rows      =", nrow(d$mcem), "  last fhat =", utils::tail(d$mcem$fhat, 1), "\n")
  flag <- intersect(c("m_step", "final_estep"), colnames(d$mcem))
  cat("  final-row flag column:", flag, "=", utils::tail(d$mcem[[flag]], 1), "\n")
}

cat("#### BDI: final E-step fails (3 iterations succeed, 4th call = final) ####\n")
undo <- patch(".augment_tree_bdi", 3L)
f <- estimate_rates(brts, method = "mcem", model = "cr",
                    control = list(lower_bound = lb, upper_bound = ub,
                                   sampling = "bdi", num_trees = 30L,
                                   max_iter = 3L, num_threads = 1L))
undo()
report("bdi, final E-step failed", f)

cat("\n#### thinning: final E-step fails (3 iterations succeed, 4th = final) ####\n")
undo <- patch("em_cpp", 3L)
f2 <- estimate_rates(brts, method = "mcem", model = "cr",
                     control = list(lower_bound = lb, upper_bound = ub,
                                    sampling = "dynamic_fresh", num_trees = 30L,
                                    max_iter = 3L, num_threads = 1L))
undo()
report("thinning, final E-step failed", f2)

cat("\n#### healthy runs, for the same fields ####\n")
f3 <- estimate_rates(brts, method = "mcem", model = "cr",
                     control = list(lower_bound = lb, upper_bound = ub,
                                    sampling = "bdi", num_trees = 30L,
                                    max_iter = 3L, num_threads = 1L))
report("bdi ok", f3)
f4 <- estimate_rates(brts, method = "mcem", model = "cr",
                     control = list(lower_bound = lb, upper_bound = ub,
                                    sampling = "dynamic_fresh", num_trees = 30L,
                                    max_iter = 3L, num_threads = 1L))
report("thinning ok", f4)

cat("\n#### diagnose_mcem: does it count the final E-step row as an iteration?\n")
for (nm in c("bdi ok", "thinning ok")) {
  f <- if (nm == "bdi ok") f3 else f4
  dg <- diagnose_mcem(f, plot = FALSE)
  cat(" ", nm, ": nrow(convergence) =", nrow(dg$convergence),
      " fit$iterations =", f$iterations,
      " last delta_max =", utils::tail(dg$convergence$delta_max, 1),
      " rejected col =", paste(utils::tail(dg$convergence$rejected, 2), collapse = ","), "\n")
}

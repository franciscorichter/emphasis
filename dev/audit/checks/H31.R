# H31/H32: all-fail MCEM path -> fit object, print, and the .mcem_warn_estep diagnostic
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(2)
tr <- ape::rcoal(15); brts <- sort(ape::branching.times(tr), decreasing = TRUE); brts <- brts/max(brts)*5
run <- function(rho) {
  w <- NULL
  fit <- withCallingHandlers(
    estimate_rates(brts, model = "cr", link = "linear", method = "mcem",
                   control = list(lower_bound = c(0, 0.1), upper_bound = c(0, 0.1),   # lambda fixed at 0 -> every tree has zero weight
                                  sampling = "dynamic_fresh", sample_size = 10L, maxN = 100L,
                                  max_iter = 20L, max_time = 60, rho = rho)),
    warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
  cat("---- rho =", rho, "----\n"); cat("WARNINGS:\n", paste(w, collapse = "\n"), "\n")
  cat("loglik:", fit$loglik, " AIC:", fit$AIC, "\n")
  cat("details$iterations: "); str(fit$details$iterations)
  cat("details$stop_reason:", fit$details$stop_reason, "\n")
  cat("names(fit):", paste(names(fit), collapse=","), "\n"); print(fit)
  invisible(fit)
}
f1 <- run(1.0)
f2 <- try(run(0.5))
# what the failing configuration actually reports (compare with the diagnostic's rho=1 re-run)
lb8 <- emphasis:::.expand_pars(c(0,0.1), c(0L,0L,0L))
e <- tryCatch(emphasis:::em_cpp(brts, lb8, 10L, 100L, 1e4, 1e6, lb8, lb8, 1e-3, 1L, FALSE, c(0L,0L,0L), 0L, 0.5), error = function(e) conditionMessage(e))
cat("direct em_cpp error (rho=0.5):", e, "\n")

# (c) the thinning maxN ratchet: does it grow on failure and never reset?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)
ns <- asNamespace("emphasis")
set.seed(4)
brts <- sort(ape::branching.times(ape::rcoal(12)), decreasing = TRUE)

fail_on <- function(which_calls) {
  orig <- get("em_cpp", envir = ns); i <- 0L
  new <- function(...) { i <<- i + 1L
    if (i %in% which_calls) stop("injected E-step failure") else orig(...) }
  unlockBinding("em_cpp", ns); assign("em_cpp", new, envir = ns); lockBinding("em_cpp", ns)
  invisible(function() { unlockBinding("em_cpp", ns); assign("em_cpp", orig, envir = ns); lockBinding("em_cpp", ns) })
}

undo <- fail_on(c(1L, 2L, 5L))
r <- emphasis:::.mcem_dynamic_fresh(
  brts = brts, pars = c(1, 0, 0, 0, 0.5, 0, 0, 0),
  sample_size = 20L, maxN = 500L, max_missing = 1e4,
  lower_bound = rep(0, 8), upper_bound = c(3, 0, 0, 0, 3, 0, 0, 0),
  max_iter = 8L, xtol = 1e-3, patience = 3L, num_threads = 1L,
  model = c(0L, 0L, 0L), link = 0L)
undo()
cat("start maxN = 500, end maxN =", r$maxN, " (expect 2000 after 3 failures)\n")
cat("n_failed =", r$n_failed, " iterations =", r$iterations, " stop =", r$stop_reason, "\n")
cat("trace maxN column =", paste(r$mcem$maxN, collapse = ","), "\n")
cat("trace n_rejected  =", paste(r$mcem$n_rejected, collapse = ","), "\n")
cat("final_IS$n_rejected =", r$final_IS$n_rejected, "\n")
cat("does n_rejected include rejected_nonfinite? columns present:",
    paste(setdiff(colnames(r$mcem), grep("^par", colnames(r$mcem), value = TRUE)),
          collapse = ","), "\n")

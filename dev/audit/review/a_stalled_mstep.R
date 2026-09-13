# (a) Seam: an M-step that returns its starting point unchanged.
# .mcem_bdi refuses to count such an iteration toward patience (item 1.5,
# H11(d)).  .mcem_dynamic_fresh has no equivalent guard, so delta = 0 counts
# and the fit reports "converged" at the point where the M-step stalled.
# The stall is injected by wrapping the C++ entry in this process only.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)
ns <- asNamespace("emphasis")
set.seed(4)
brts <- sort(ape::branching.times(ape::rcoal(15)), decreasing = TRUE)
lb <- c(0, 0); ub <- c(3, 3)

stall <- function(name, arg) {
  orig <- get(name, envir = ns)
  new <- function(...) {
    a <- list(...)
    r <- orig(...)
    r$estimates <- as.numeric(a[[arg]])   # M-step returns its start
    r
  }
  unlockBinding(name, ns); assign(name, new, envir = ns); lockBinding(name, ns)
  invisible(function() { unlockBinding(name, ns); assign(name, orig, envir = ns); lockBinding(name, ns) })
}

fit <- function(samp) estimate_rates(brts, method = "mcem", model = "cr",
  control = list(lower_bound = lb, upper_bound = ub, sampling = samp,
                 num_trees = 30L, max_iter = 8L, num_threads = 1L))

undo <- stall("m_cpp", "init_pars")
f <- fit("bdi"); undo()
cat("BDI, stalled M-step      : stop =", f$stop_reason, " iterations =", f$iterations,
    " delta_max =", paste(round(f$details$mcem$delta_max, 6), collapse = ","),
    " m_moved =", paste(f$details$mcem$m_moved, collapse = ","), "\n")

undo <- stall("em_cpp", "init_pars")
f2 <- fit("dynamic_fresh"); undo()
cat("thinning, stalled M-step : stop =", f2$stop_reason, " iterations =", f2$iterations,
    " delta_max =", paste(round(f2$details$mcem$delta_max, 6), collapse = ","), "\n")
cat("  pars returned =", round(f2$pars, 6), " (init was the box midpoint",
    round((lb + ub) / 2, 6), "with mu halved)\n")

# H30: does a convergence streak survive an interleaved E-step failure?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(1)
tr <- ape::rcoal(20); brts <- sort(ape::branching.times(tr), decreasing = TRUE)
brts <- brts / max(brts) * 5
lb8 <- emphasis:::.expand_pars(c(0, 0),   c(0L,0L,0L))
ub8 <- emphasis:::.expand_pars(c(2, 1),   c(0L,0L,0L))
ip8 <- emphasis:::.expand_pars(c(0.5, 0.1), c(0L,0L,0L))
flag <- new.env(); flag$fail <- FALSE; flag$ncall <- 0L
cond <- function(p) { flag$ncall <- flag$ncall + 1L
  if (flag$fail) { flag$fail <- FALSE; stop("injected E-step failure") }
  0 }
msgs <- character()
res <- withCallingHandlers(
  emphasis:::.mcem_dynamic_fresh(brts, ip8, sample_size = 20L, maxN = 2000L, max_missing = 1e4,
      lower_bound = lb8, upper_bound = ub8, max_iter = 20L, xtol = 1e-3,
      tol = 1,            # every successful iteration counts as "stable"
      patience = 3L, num_threads = 1L, verbose = TRUE, conditional = cond,
      model = c(0L,0L,0L), link = 0L, max_time = 120),
  message = function(m) { txt <- conditionMessage(m); msgs <<- c(msgs, txt)
    if (grepl("^Iteration 2:", txt)) flag$fail <- TRUE
    invokeRestart("muffleMessage") })
cat(msgs, sep = "")
cat("stop_reason:", res$stop_reason, " iterations recorded:", res$iterations, "\n")
m <- res$mcem
center <- (lb8 + ub8)/2; range_vec <- ub8 - lb8; range_vec[range_vec == 0] <- 1
p2 <- unlist(m[2, grep("^par", names(m))]); p3 <- unlist(m[3, grep("^par", names(m))])
pert <- pmin(pmax(0.8*p2 + 0.2*center, lb8), ub8)
cat(sprintf("delta_max reported at 3rd success: %.3e\n", m$delta_max[3]))
cat(sprintf("delta vs previous success (p2):    %.3e\n", max(abs(p3 - p2)/range_vec)))
cat(sprintf("delta vs perturbed iterate:        %.3e\n", max(abs(p3 - pert)/range_vec)))
cat("Streak survived failure? ", res$stop_reason == "converged" && nrow(m) == 3L, "\n")

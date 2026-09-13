# H23 (companion): the 120 s check runs only BETWEEN attempts (E_step.cpp:71-78),
# so a single augmentation that inserts up to max_missing lineages cannot be
# interrupted.  Measure one attempt (sample_size = 1, maxN = 1) under
# lambda = mu = 40 on a 30-tip tree with increasing max_missing; every attempt
# overruns, so the E-step throws after exactly one augmentation and its
# elapsed time is the time of that one augmentation.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
set.seed(30)
tr <- TreeSim::sim.bd.taxa(30, 1, 0.5, 0.1, complete = FALSE)[[1]]
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
pars <- c(40, 0, 0, 0, 40, 0, 0, 0)
for (mm in c(5000L, 10000L, 20000L, 40000L)) {
  t1 <- system.time(res <- tryCatch(
    emphasis:::augment_trees(brts, pars, sample_size = 1L, maxN = 1L,
                              max_missing = mm, max_lambda = 1e9, num_threads = 1L),
    error = function(e) e))[3]
  msg <- if (inherits(res, "error")) conditionMessage(res) else "ACCEPTED"
  cat(sprintf("max_missing=%6d  one attempt: %7.2f s   %s\n", mm, t1, msg))
}

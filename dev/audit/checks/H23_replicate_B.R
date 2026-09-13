# H23 (B only): forced E-step failure (sample_size > maxN) via .mcem_dynamic_fresh
# with max_time = 1 s. Claim: the failure branch `next`s before the max_time
# check, doubles maxN, perturbs pars 20% toward the box centre, and after 8
# failures returns the perturbed vector.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
set.seed(7)
tr <- TreeSim::sim.bd.taxa(25, 1, 0.5, 0.2, complete = FALSE)[[1]]
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
lb <- c(0.01, 0, 0, 0, 0.0, 0, 0, 0); ub <- c(2, 0, 0, 0, 1, 0, 0, 0)
p0 <- c(0.5, 0, 0, 0, 0.2, 0, 0, 0)
tB <- system.time(rB <- suppressWarnings(emphasis:::.mcem_dynamic_fresh(
  brts, p0, sample_size = 100000L, maxN = 50L, max_missing = 10000L,
  lower_bound = lb, upper_bound = ub, max_iter = 200L, xtol = 1e-3, tol = 1e-3,
  patience = 3L, num_threads = 1L, verbose = FALSE, model = c(0L, 0L, 0L),
  link = 0L, max_time = 1)))[3]
ctr <- (lb + ub) / 2
expected <- p0; for (k in 1:8) expected <- pmin(pmax(0.8 * expected + 0.2 * ctr, lb), ub)
cat(sprintf("[B] elapsed %.1f s with max_time = 1 s: stop_reason = %s (time_budget never triggered: %s)\n",
            tB, rB$stop_reason, rB$stop_reason != "time_budget"))
cat(sprintf("    returned (la, mu) = (%.4f, %.4f); 8x centre-perturbed prediction = (%.4f, %.4f); start = (%.2f, %.2f); match: %s\n",
            rB$pars[1], rB$pars[5], expected[1], expected[5], p0[1], p0[5],
            isTRUE(all.equal(rB$pars, expected))))

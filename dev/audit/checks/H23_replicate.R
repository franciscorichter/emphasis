# H23 replication (independent). Varies what the verifier did not:
#  (A) different tree (20 tips, la=0.8, mu=0.4), num_threads = 2, N = 3e6,
#      maxN = 6e6: does the C++ E-step still throw "maxN exceeded" at ~120 s
#      with attempts << maxN and zero rejections?
#  (B) fast forced failure (sample_size > maxN, no timeout needed) through
#      .mcem_dynamic_fresh with max_time = 1: is the R loop's max_time check
#      skipped on the failure branch, are pars perturbed toward the centre,
#      and is maxN doubled?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))

# ---- (B) first: cheap ---------------------------------------------------------
set.seed(7)
tr <- TreeSim::sim.bd.taxa(25, 1, 0.5, 0.2, complete = FALSE)[[1]]
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
lb <- c(0.01, 0, 0, 0, 0.0, 0, 0, 0); ub <- c(2, 0, 0, 0, 1, 0, 0, 0)
p0 <- c(0.5, 0, 0, 0, 0.2, 0, 0, 0)
tB <- system.time(rB <- emphasis:::.mcem_dynamic_fresh(
  brts, p0, sample_size = 100000L, maxN = 50L, max_missing = 10000L,
  lower_bound = lb, upper_bound = ub, max_iter = 200L, xtol = 1e-3, tol = 1e-3,
  patience = 3L, num_threads = 1L, verbose = TRUE, model = c(0L, 0L, 0L),
  link = 0L, max_time = 1))[3]
ctr <- (lb + ub) / 2
expected <- p0; for (k in 1:8) expected <- pmin(pmax(0.8 * expected + 0.2 * ctr, lb), ub)
cat(sprintf("[B] elapsed %.1f s (max_time = 1): stop_reason = %s, iterations = %d\n",
            tB, rB$stop_reason, rB$iterations))
cat(sprintf("    returned pars (la, mu) = (%.4f, %.4f); 8x centre-perturbed prediction = (%.4f, %.4f); start = (%.2f, %.2f)\n",
            rB$pars[1], rB$pars[5], expected[1], expected[5], p0[1], p0[5]))

# ---- (A) the cap, different tree / threads ----------------------------------
set.seed(20)
tr <- TreeSim::sim.bd.taxa(20, 1, 0.8, 0.4, complete = FALSE)[[1]]
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
pars <- c(0.8, 0, 0, 0, 0.4, 0, 0, 0)
N <- 3000000L; maxN <- 6000000L
tA <- system.time(res <- tryCatch(
  emphasis:::augment_trees(brts, pars, sample_size = N, maxN = maxN,
                            max_missing = 10000L, max_lambda = 1e6, num_threads = 2L),
  error = function(e) e))[3]
if (inherits(res, "error")) {
  msg <- conditionMessage(res)
  trees <- as.integer(sub(".*Trees so far: (\\d+).*", "\\1", msg))
  rej <- sum(as.integer(regmatches(msg, gregexpr("\\d+(?= (lambda|overruns|zero weights|unhandled))", msg, perl = TRUE))[[1]]))
  cat(sprintf("[A] 20 tips, 2 threads, N=%d, maxN=%d: THREW after %.1f s\n    %s\n    accepted=%d rejected=%d attempts=%d (%.1f%% of maxN); mentions time: %s\n",
              N, maxN, tA, msg, trees, rej, trees + rej, 100 * (trees + rej) / maxN,
              grepl("time", msg, ignore.case = TRUE)))
} else {
  cat(sprintf("[A] did NOT throw: %d trees in %.1f s\n", length(res$logf), tA))
}

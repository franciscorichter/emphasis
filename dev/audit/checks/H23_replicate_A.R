# H23 (A only): the 120 s cap on a tree the verifier did not use, with 2 threads.
# 100 tips, la = 0.6, mu = 0.45 (mu/la = 0.75, slow attempts), N = 1e6, maxN = 2e6.
# Expected if the claim holds: throw at ~120 s, "maxN exceeded ..." with
# attempts << maxN and no time word in the message.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
set.seed(100)
tr <- TreeSim::sim.bd.taxa(100, 1, 0.6, 0.45, complete = FALSE)[[1]]
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
pars <- c(0.6, 0, 0, 0, 0.45, 0, 0, 0)
# per-attempt cost first
t1 <- system.time(r1 <- tryCatch(emphasis:::augment_trees(brts, pars, sample_size = 20L, maxN = 200L,
        max_missing = 10000L, max_lambda = 1e6, num_threads = 1L), error = function(e) conditionMessage(e)))[3]
cat(sprintf("[A0] 100 tips: 20 accepted trees in %.2f s (rej0=%s)\n", t1,
            if (is.character(r1)) r1 else r1$rejected_zero_weights))
N <- 1000000L; maxN <- 2000000L
tA <- system.time(res <- tryCatch(
  emphasis:::augment_trees(brts, pars, sample_size = N, maxN = maxN,
                            max_missing = 10000L, max_lambda = 1e6, num_threads = 2L),
  error = function(e) e))[3]
if (inherits(res, "error")) {
  msg <- conditionMessage(res)
  trees <- as.integer(sub(".*Trees so far: (\\d+).*", "\\1", msg))
  rej <- sum(as.integer(regmatches(msg, gregexpr("\\d+(?= (lambda|overruns|zero weights|unhandled))", msg, perl = TRUE))[[1]]))
  cat(sprintf("[A] 100 tips, 2 threads, N=%d, maxN=%d: THREW after %.1f s\n    %s\n    accepted=%d rejected=%d attempts=%d (%.2f%% of maxN); mentions time: %s\n",
              N, maxN, tA, msg, trees, rej, trees + rej, 100 * (trees + rej) / maxN,
              grepl("time", msg, ignore.case = TRUE)))
} else {
  cat(sprintf("[A] did NOT throw: %d trees in %.1f s\n", length(res$logf), tA))
}

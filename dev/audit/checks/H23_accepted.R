# H23 (companion): same cap when trees ARE being accepted (no rejection at all).
# 30-tip tree, ~0.5 ms per accepted tree, sample_size = 1e6, maxN = 2e6.
# At that rate 1e6 trees need ~500 s; the cap should fire at ~120 s with
# ~2e5 accepted trees, attempts << maxN, message "maxN exceeded ... 0 zero weights".
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
set.seed(30)
tr <- TreeSim::sim.bd.taxa(30, 1, 0.5, 0.1, complete = FALSE)[[1]]
brts <- sort(ape::branching.times(tr), decreasing = TRUE)
pars <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)
N <- 1000000L; maxN <- 2000000L
t1 <- system.time(res <- tryCatch(
  emphasis:::augment_trees(brts, pars, sample_size = N, maxN = maxN,
                            max_missing = 10000L, max_lambda = 1e6, num_threads = 1L),
  error = function(e) e))[3]
stopifnot(inherits(res, "error"))
msg <- conditionMessage(res)
trees <- as.integer(sub(".*Trees so far: (\\d+).*", "\\1", msg))
rej <- sum(as.integer(regmatches(msg, gregexpr("\\d+(?= (lambda|overruns|zero weights|unhandled))", msg, perl = TRUE))[[1]]))
cat(sprintf("30 tips, N=%d, maxN=%d: threw after %.1f s\n  message: %s\n  accepted=%d rejected=%d attempts=%d (%.1f%% of maxN) -> maxN NOT exceeded\n",
            N, maxN, t1, msg, trees, rej, trees + rej, 100 * (trees + rej) / maxN))

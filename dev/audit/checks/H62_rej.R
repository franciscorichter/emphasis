## H62 supplement: in overflow runs, are rejections zero (so bias == log(M/N))
## and does M < maxN come from iterations skipped while `stop` was briefly true?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
brts  <- c(4, 2.5, 1.2, 0.6); pars8 <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)
one <- function(nt) {
  r <- emphasis:::augment_trees(brts, pars8, sample_size = 50L, maxN = 5000L,
                                max_missing = 1000L, max_lambda = 1e6,
                                num_threads = nt, model = c(0L,0L,0L), link = 0L)
  c(M = length(r$logf), rzw = r$rejected_zero_weights, rej = r$rejected,
    ovr = r$rejected_overruns, lam = r$rejected_lambda,
    attempts_accounted = length(r$logf) + r$rejected_zero_weights + r$rejected +
                         r$rejected_overruns + r$rejected_lambda)
}
res <- t(replicate(25, one(8L)))
print(res)
cat(sprintf("all rejection counters zero: %s; M < maxN in %d/25 runs => %d..%d of 5000 iterations were skipped (not rejected)\n",
            all(res[, c("rzw","rej","ovr","lam")] == 0), sum(res[, "M"] < 5000),
            min(5000 - res[, "M"]), max(5000 - res[, "M"])))

# H23: the C++ E-step (src/E_step.cpp:73) enforces a 120 s wall-clock cap when
# max_time_seconds == 0 (the header says 0 = "no time limit"); no caller passes
# the argument; on timeout with < N trees the throw reads "maxN exceeded ...".
# Test: a 300-tip CR tree where every thinning attempt takes ~1.6 ms, with
# maxN = 200000 (~316 s of work at that rate).  If the cap exists, the call
# throws after ~120 s with total attempts << maxN and a message that says
# "maxN exceeded".
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))

mk <- function(n, la, mu) {
  tr <- TreeSim::sim.bd.taxa(n, 1, la, mu, complete = FALSE)[[1]]
  sort(ape::branching.times(tr), decreasing = TRUE)
}
parse_msg <- function(msg) {
  g <- function(p) as.integer(sub(paste0(".*?(\\d+) ", p, ".*"), "\\1", msg))
  c(lambda = g("lambda"), overruns = g("overruns"), zero = g("zero weights"),
    unhandled = g("unhandled"), trees = as.integer(sub(".*Trees so far: (\\d+).*", "\\1", msg)))
}

# --- (0) validation-sized E-step, for scale ---------------------------------
set.seed(30); brts30 <- mk(30, 0.5, 0.1)
p30 <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)
t30 <- system.time(r30 <- emphasis:::augment_trees(brts30, p30, sample_size = 200L, maxN = 2000L,
                     max_missing = 10000L, max_lambda = 1e6, num_threads = 1L))[3]
cat(sprintf("[0] 30 tips, N=200: %.3f s per E-step (rej0=%d) -> cap is %.0fx away\n",
            t30, r30$rejected_zero_weights, 120 / t30))

# --- (1) the timeout ---------------------------------------------------------
set.seed(300); brts <- mk(300, 0.6, 0.45)
pars <- c(0.6, 0, 0, 0, 0.45, 0, 0, 0)
maxN <- 200000L
t1 <- system.time(res <- tryCatch(
  emphasis:::augment_trees(brts, pars, sample_size = 200L, maxN = maxN,
                            max_missing = 10000L, max_lambda = 1e6, num_threads = 1L),
  error = function(e) e))[3]
stopifnot(inherits(res, "error"))
msg <- conditionMessage(res)
cnt <- parse_msg(msg)
attempts <- sum(cnt)
cat(sprintf("[1] 300 tips, maxN=%d: threw after %.1f s\n    message: %s\n    attempts made = %d (%.1f%% of maxN)\n",
            maxN, t1, msg, attempts, 100 * attempts / maxN))
cat(sprintf("    says 'maxN exceeded': %s | attempts < maxN: %s | elapsed in [110,140] s: %s | mentions time/timeout: %s\n",
            grepl("maxN exceeded", msg), attempts < maxN, t1 >= 110 && t1 <= 140,
            grepl("time", msg, ignore.case = TRUE)))

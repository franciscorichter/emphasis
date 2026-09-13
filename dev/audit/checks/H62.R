## H62 — E_step stop-flag race: `if (!stop)` is read outside the mutex and
## `stop = (trees.size() == N)` uses `==`, so a second thread that has already
## passed the check can push tree N+1 and flip `stop` back to false, after which
## every thread runs to maxN (or the 120 s cap).  `S_completed = N + rzw` then
## under-counts the sample and fhat is biased by ~log((M+rzw)/(N+rzw)).
##
## Test: repeated augment_trees / em_cpp calls with num_threads = 8; assert
## length(logf) == sample_size.  num_threads = 1 is the control.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))

brts  <- c(4, 2.5, 1.2, 0.6)                 # 4-tip CR tree (fast augmentation)
pars8 <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)       # lambda = 0.5, mu = 0.1
cr    <- c(0L, 0L, 0L)

aug <- function(N, maxN, nt) {
  t0 <- proc.time()[["elapsed"]]
  r  <- emphasis:::augment_trees(brts, pars8, sample_size = N, maxN = maxN,
                      max_missing = 1000L, max_lambda = 1e6,
                      num_threads = nt, model = cr, link = 0L)
  M   <- length(r$logf)
  rzw <- r$rejected_zero_weights
  lw  <- r$logf - r$logg
  m   <- max(lw)
  fhat_code    <- log(sum(exp(lw - m)) / (N + rzw)) + m   # what E_step.cpp:143-144 computes
  fhat_correct <- log(sum(exp(lw - m)) / (M + rzw)) + m   # denominator = trees actually summed
  c(M = M, rzw = rzw, rej = r$rejected, ovr = r$rejected_overruns,
    lam = r$rejected_lambda, fhat_code = fhat_code, fhat_correct = fhat_correct,
    bias = fhat_code - fhat_correct, sec = proc.time()[["elapsed"]] - t0)
}

run_block <- function(label, N, maxN, nt, reps) {
  res <- t(replicate(reps, aug(N, maxN, nt)))
  over <- res[, "M"] != N
  cat(sprintf("\n[%s] N=%d maxN=%d threads=%d reps=%d\n", label, N, maxN, nt, reps))
  cat(sprintf("  runs with length(logf) != N : %d / %d  (%.1f%%)\n",
              sum(over), reps, 100 * mean(over)))
  cat("  table of M = length(logf):\n")
  print(table(res[, "M"]))
  if (any(over)) {
    cat(sprintf("  overflow runs: M range %d..%d; fhat bias (code - correct) range %.4f..%.4f; predicted log(M/N) range %.4f..%.4f\n",
                min(res[over, "M"]), max(res[over, "M"]),
                min(res[over, "bias"]), max(res[over, "bias"]),
                min(log(res[over, "M"] / N)), max(log(res[over, "M"] / N))))
    cat(sprintf("  elapsed: overflow runs mean %.4f s vs clean runs mean %.4f s\n",
                mean(res[over, "sec"]), if (any(!over)) mean(res[!over, "sec"]) else NA))
  }
  invisible(res)
}

cat("hardware threads:", parallel::detectCores(), "\n")

## (a) the hypothesis' own recipe
r_ctrl <- run_block("control", N = 50L, maxN = 5000L, nt = 1L, reps = 40L)
r_8    <- run_block("race",    N = 50L, maxN = 5000L, nt = 8L, reps = 150L)

## (b) package defaults (sample_size 200, maxN 2000) with the default num_threads = 1
##     and with 8 threads
r_def1 <- run_block("defaults-1thr", N = 200L, maxN = 2000L, nt = 1L, reps = 30L)
r_def8 <- run_block("defaults-8thr", N = 200L, maxN = 2000L, nt = 8L, reps = 60L)

## (c) full em_cpp (E-step + M-step): does `trees` exceed sample_size, and what
##     happens to fhat?
lb <- c(0.01, 0, 0, 0, 0.001, 0, 0, 0)
ub <- c(5,    0, 0, 0, 5,     0, 0, 0)
em1 <- function(nt) {
  r <- emphasis:::em_cpp(brts = brts, init_pars = pars8, sample_size = 50L, maxN = 5000L,
              max_missing = 1000L, max_lambda = 1e6,
              lower_bound = lb, upper_bound = ub, xtol_rel = 1e-3,
              num_threads = nt, copy_trees = FALSE, model = cr, link = 0L,
              rho = 1, rconditional = NULL)
  c(trees = r$trees, fhat = r$fhat, rzw = r$rejected_zero_weights,
    lam = r$estimates[1], mu = r$estimates[5])
}
e1 <- t(replicate(20, em1(1L)))
e8 <- t(replicate(60, em1(8L)))
cat("\n[em_cpp] num_threads=1: trees table\n"); print(table(e1[, "trees"]))
cat(sprintf("  fhat mean %.4f sd %.4f\n", mean(e1[, "fhat"]), sd(e1[, "fhat"])))
cat("[em_cpp] num_threads=8: trees table\n"); print(table(e8[, "trees"]))
ov <- e8[, "trees"] != 50
cat(sprintf("  overflow runs %d/%d; fhat in overflow runs: mean %.4f (min %.4f); fhat in clean runs: mean %.4f\n",
            sum(ov), nrow(e8),
            if (any(ov)) mean(e8[ov, "fhat"]) else NA,
            if (any(ov)) min(e8[ov, "fhat"]) else NA,
            if (any(!ov)) mean(e8[!ov, "fhat"]) else NA))
if (any(ov)) {
  cat(sprintf("  predicted bias log(trees/50) in overflow runs: %.4f..%.4f\n",
              min(log(e8[ov, "trees"] / 50)), max(log(e8[ov, "trees"] / 50))))
  cat(sprintf("  observed fhat(overflow) - mean fhat(1 thread): %.4f..%.4f\n",
              min(e8[ov, "fhat"]) - mean(e1[, "fhat"]),
              max(e8[ov, "fhat"]) - mean(e1[, "fhat"])))
}

## verdict line
cat(sprintf("\nVERDICT: augment_trees 8-thr overflow rate %.1f%%; defaults 8-thr %.1f%%; 1-thr %.1f%%; em_cpp 8-thr %.1f%%\n",
            100 * mean(r_8[, "M"] != 50), 100 * mean(r_def8[, "M"] != 200),
            100 * mean(c(r_ctrl[, "M"] != 50, r_def1[, "M"] != 200)), 100 * mean(ov)))

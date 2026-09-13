## H22 replication part 2: intermediate box, exp-link replicates, 8-failure return value.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
options(width = 150)
fit1 <- function(tr, init, lb, ub, num_trees, maxN, max_iter, link = "linear") {
  msgs <- character(0)
  fit <- withCallingHandlers(
    suppressWarnings(estimate_rates(tr, model = "cr", method = "mcem", init_pars = init, link = link,
      control = list(sampling = "dynamic_fresh", num_trees = num_trees, maxN = maxN,
                     max_iter = max_iter, tol = 1e-3, patience = 3L, num_threads = 1L,
                     lower_bound = lb, upper_bound = ub, verbose = TRUE))),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  d <- fit$details
  c(lambda = unname(fit$pars[1]), mu = unname(fit$pars[2]), loglik = fit$loglik,
    n_success = if (is.null(d$iterations)) 0 else d$iterations,
    n_fail = sum(grepl("E-step failed", msgs)), stop = d$stop_reason)
}
set.seed(22); tr20 <- rphylo(20, 1, 0.3)

cat("=== 4. README-style box (0,0)-(2,1), centre (1,0.5), linear, tr20, 3 reps each ===\n")
t0 <- proc.time()[3]
res4 <- rbind(
  clean1 = fit1(tr20, c(1, 0.3), c(0,0), c(2,1), 200L, 4000L, 30L),
  clean2 = fit1(tr20, c(1, 0.3), c(0,0), c(2,1), 200L, 4000L, 30L),
  clean3 = fit1(tr20, c(1, 0.3), c(0,0), c(2,1), 200L, 4000L, 30L),
  alt1   = fit1(tr20, c(1, 0.3), c(0,0), c(2,1), 200L, 150L,  60L),
  alt2   = fit1(tr20, c(1, 0.3), c(0,0), c(2,1), 200L, 150L,  60L),
  alt3   = fit1(tr20, c(1, 0.3), c(0,0), c(2,1), 200L, 150L,  60L))
print(res4); cat(sprintf("(elapsed %.0fs)\n", proc.time()[3] - t0))

cat("\n=== 5. exponential link, box (-3,-3)-(1,1), centre (-1,-1) [rates 0.368], 3 reps each ===\n")
t0 <- proc.time()[3]
res5 <- rbind(
  clean1 = fit1(tr20, c(0, -1.2), c(-3,-3), c(1,1), 200L, 4000L, 30L, "exponential"),
  clean2 = fit1(tr20, c(0, -1.2), c(-3,-3), c(1,1), 200L, 4000L, 30L, "exponential"),
  clean3 = fit1(tr20, c(0, -1.2), c(-3,-3), c(1,1), 200L, 4000L, 30L, "exponential"),
  alt1   = fit1(tr20, c(0, -1.2), c(-3,-3), c(1,1), 200L, 150L,  60L, "exponential"),
  alt2   = fit1(tr20, c(0, -1.2), c(-3,-3), c(1,1), 200L, 150L,  60L, "exponential"),
  alt3   = fit1(tr20, c(0, -1.2), c(-3,-3), c(1,1), 200L, 150L,  60L, "exponential"))
print(res5); cat(sprintf("(elapsed %.0fs)\n", proc.time()[3] - t0))

cat("\n=== 6. exponential link, WIDE box (-6,-6)-(3,3), centre (-1.5,-1.5) [rates 0.22], 2 reps each ===\n")
t0 <- proc.time()[3]
res6 <- rbind(
  clean1 = fit1(tr20, c(0, -1.2), c(-6,-6), c(3,3), 200L, 4000L, 30L, "exponential"),
  clean2 = fit1(tr20, c(0, -1.2), c(-6,-6), c(3,3), 200L, 4000L, 30L, "exponential"),
  alt1   = fit1(tr20, c(0, -1.2), c(-6,-6), c(3,3), 200L, 150L,  60L, "exponential"),
  alt2   = fit1(tr20, c(0, -1.2), c(-6,-6), c(3,3), 200L, 150L,  60L, "exponential"))
print(res6); cat(sprintf("(elapsed %.0fs)\n", proc.time()[3] - t0))

cat("\n=== 7. eight consecutive failures (6-tip tree, num_trees=50001 > maxN cap 50000) ===\n")
set.seed(3); tr6 <- rphylo(6, 1, 0.3)
t0 <- proc.time()[3]
warnC <- character(0)
fitC <- withCallingHandlers(
  suppressMessages(estimate_rates(tr6, model = "cr", method = "mcem", init_pars = c(1, 0.3),
    control = list(sampling = "dynamic_fresh", num_trees = 50001L, maxN = 50000L, max_iter = 20L,
                   num_threads = 1L, lower_bound = c(0,0), upper_bound = c(4,4)))),
  warning = function(w) { warnC <<- c(warnC, conditionMessage(w)); invokeRestart("muffleWarning") })
expected <- 0.8^8 * c(1, 0.3) + (1 - 0.8^8) * c(2, 2)
cat(sprintf("elapsed %.0fs; stop=%s; is.null(details$iterations)=%s; pars=(%.6f, %.6f); expected=(%.6f, %.6f); maxdiff=%.1e\n",
    proc.time()[3] - t0, fitC$details$stop_reason, is.null(fitC$details$iterations),
    fitC$pars[1], fitC$pars[2], expected[1], expected[2], max(abs(fitC$pars - expected))))
cat(sprintf("loglik=%s AIC=%s; warning: %s\n", fitC$loglik, fitC$AIC, substr(paste(warnC, collapse=" | "), 1, 200)))

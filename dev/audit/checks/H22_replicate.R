## H22 replication: vary tree, box centre direction, and link.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(ape)})
options(width = 150)

fit1 <- function(tr, init, lb, ub, num_trees, maxN, max_iter, link = "linear", verbose = FALSE) {
  msgs <- character(0)
  fit <- withCallingHandlers(
    suppressWarnings(estimate_rates(tr, model = "cr", method = "mcem", init_pars = init, link = link,
      control = list(sampling = "dynamic_fresh", num_trees = num_trees, maxN = maxN,
                     max_iter = max_iter, tol = 1e-3, patience = 3L, num_threads = 1L,
                     lower_bound = lb, upper_bound = ub, verbose = verbose))),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  d <- fit$details
  c(lambda = fit$pars[1], mu = fit$pars[2], loglik = fit$loglik,
    n_success = if (is.null(d$iterations)) 0 else d$iterations,
    n_fail = sum(grepl("E-step failed", msgs)), stop = d$stop_reason,
    max_rej = if (!is.null(d$mcem)) max(d$mcem$rejected) else NA)
}

cat("=== 0. guards: any check of maxN vs sample_size? ===\n")
ctrl <- emphasis:::.resolve_control_aliases(
  modifyList(emphasis:::estimate_rates_control("mcem"), list(num_trees = 3000L)), "mcem", list(num_trees = 3000L))
cat(sprintf("num_trees=3000 -> sample_size=%d maxN=%d\n", ctrl$sample_size, ctrl$maxN))
r <- tryCatch(emphasis:::em_cpp(brts = c(2, 1, 0.5), init_pars = c(1, 0.3), sample_size = 30L, maxN = 20L,
  max_missing = 1e4, max_lambda = 1e6, lower_bound = c(0,0), upper_bound = c(4,4), xtol_rel = 1e-3,
  num_threads = 1L, copy_trees = FALSE), error = function(e) conditionMessage(e))
cat("em_cpp sample_size=30 maxN=20 ->", if (is.character(r)) paste("ERROR:", r) else "returned", "\n")

cat("\n=== 1. different tree: 30 tips, rphylo(30, 0.8, 0.4) seed 7; box (0,0)-(3,3), centre (1.5,1.5) ===\n")
set.seed(7); tr30 <- rphylo(30, 0.8, 0.4)
b30 <- sort(branching.times(tr30), decreasing = TRUE)
cat(sprintf("crown age %.3f\n", b30[1]))
ref <- DDD::bd_ML(brts = b30, cond = 0, btorph = 0, soc = 2, initparsopt = c(0.8, 0.4),
                  idparsopt = 1:2, parsfix = c(0,0), idparsfix = 3:4, verbose = FALSE)
cat(sprintf("DDD::bd_ML cond=0: lambda=%.4f mu=%.4f\n", ref$lambda0, ref$mu0))
t0 <- proc.time()[3]
res1 <- rbind(
  clean1 = fit1(tr30, c(0.8, 0.4), c(0,0), c(3,3), 200L, 4000L, 20L, verbose = TRUE),
  clean2 = fit1(tr30, c(0.8, 0.4), c(0,0), c(3,3), 200L, 4000L, 20L, verbose = TRUE),
  alt1   = fit1(tr30, c(0.8, 0.4), c(0,0), c(3,3), 200L, 150L,  40L, verbose = TRUE),
  alt2   = fit1(tr30, c(0.8, 0.4), c(0,0), c(3,3), 200L, 150L,  40L, verbose = TRUE))
print(res1); cat(sprintf("(elapsed %.0fs)\n", proc.time()[3] - t0))

cat("\n=== 2. centre BELOW the MLE: 20-tip tree (seed 22, as verifier); box (0,0)-(1.2,0.2), centre (0.6,0.1) ===\n")
set.seed(22); tr20 <- rphylo(20, 1, 0.3)
t0 <- proc.time()[3]
res2 <- rbind(
  clean = fit1(tr20, c(1, 0.15), c(0,0), c(1.2,0.2), 200L, 4000L, 20L, verbose = TRUE),
  alt   = fit1(tr20, c(1, 0.15), c(0,0), c(1.2,0.2), 200L, 150L,  40L, verbose = TRUE))
print(res2); cat(sprintf("(elapsed %.0fs)\n", proc.time()[3] - t0))

cat("\n=== 3. exponential link (log-rate parameters): 20-tip tree; box (-3,-3)-(1,1), centre (-1,-1) ===\n")
t0 <- proc.time()[3]
res3 <- rbind(
  clean = fit1(tr20, c(0, -1.2), c(-3,-3), c(1,1), 200L, 4000L, 20L, link = "exponential", verbose = TRUE),
  alt   = fit1(tr20, c(0, -1.2), c(-3,-3), c(1,1), 200L, 150L,  40L, link = "exponential", verbose = TRUE))
print(res3); cat(sprintf("(elapsed %.0fs)\n", proc.time()[3] - t0))
cat("log-scale centre (-1,-1) = rates (0.368, 0.368); clean MLE on log scale approx (0, -1.9)\n")

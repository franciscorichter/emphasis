## H62 independent replication: vary what the verifier did not.
##   (1) bigger trees (12 and 25 tips) -- slower augmentations, less mutex contention
##   (2) num_threads = 2 and 4 (not only 8)
##   (3) N close to maxN (N=100, maxN=200, grainsize 25)
##   (4) dd model on the exponential link
##   (5) is bias == log((M+rzw)/(N+rzw)) exactly, and does em_cpp hand all M trees to the M-step?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
cat("hardware threads:", parallel::detectCores(), "\n")

set.seed(11)
tr12 <- ape::rphylo(12, 0.5, 0.1); b12 <- sort(as.numeric(ape::branching.times(tr12)), decreasing = TRUE)
set.seed(12)
tr25 <- ape::rphylo(25, 0.5, 0.1); b25 <- sort(as.numeric(ape::branching.times(tr25)), decreasing = TRUE)
b4   <- c(4, 2.5, 1.2, 0.6)
cr8  <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)
dd8  <- c(log(0.6), -0.03, 0, 0, log(0.1), 0, 0, 0)   # exponential link: lambda = exp(b0 + bN*n)

aug <- function(brts, pars, N, maxN, nt, model = c(0L,0L,0L), link = 0L) {
  t0 <- proc.time()[["elapsed"]]
  r  <- emphasis:::augment_trees(brts, pars, sample_size = N, maxN = maxN,
                                 max_missing = 1000L, max_lambda = 1e6,
                                 num_threads = nt, model = model, link = link)
  M   <- length(r$logf); rzw <- r$rejected_zero_weights
  lw  <- r$logf - r$logg; m <- max(lw)
  fc  <- log(sum(exp(lw - m)) / (N + rzw)) + m
  fx  <- log(sum(exp(lw - m)) / (M + rzw)) + m
  c(M = M, rzw = rzw, other_rej = r$rejected + r$rejected_overruns + r$rejected_lambda,
    bias = fc - fx, pred = log((M + rzw) / (N + rzw)), sec = proc.time()[["elapsed"]] - t0)
}
block <- function(label, brts, pars, N, maxN, nt, reps, model = c(0L,0L,0L), link = 0L) {
  res <- t(replicate(reps, aug(brts, pars, N, maxN, nt, model, link)))
  ov  <- res[, "M"] != N
  cat(sprintf("\n[%s] tips=%d N=%d maxN=%d threads=%d reps=%d: overflow %d/%d (%.0f%%)",
              label, length(brts) + 1, N, maxN, nt, reps, sum(ov), reps, 100 * mean(ov)))
  if (any(ov)) cat(sprintf("\n   M range %d..%d; max|bias - log((M+rzw)/(N+rzw))| = %.2e; sec overflow %.3f vs clean %.3f",
                           min(res[ov, "M"]), max(res[ov, "M"]), max(abs(res[ov, "bias"] - res[ov, "pred"])),
                           mean(res[ov, "sec"]), if (any(!ov)) mean(res[!ov, "sec"]) else NA))
  if (any(res[, "M"] < N)) cat("\n   !! some run returned FEWER than N trees")
  cat("\n"); invisible(res)
}

r1  <- block("4tip-ctrl",   b4,  cr8, 50L, 5000L, 1L, 20L)
r2  <- block("4tip-2thr",   b4,  cr8, 50L, 5000L, 2L, 40L)
r3  <- block("4tip-4thr",   b4,  cr8, 50L, 5000L, 4L, 40L)
r4  <- block("12tip-4thr",  b12, cr8, 50L, 5000L, 4L, 40L)
r5  <- block("12tip-8thr",  b12, cr8, 50L, 5000L, 8L, 40L)
r6  <- block("25tip-8thr",  b25, cr8, 50L, 2000L, 8L, 30L)
r7  <- block("25tip-1thr",  b25, cr8, 50L, 2000L, 1L, 10L)
r8  <- block("4tip-N100-maxN200-8thr", b4, cr8, 100L, 200L, 8L, 40L)
r9  <- block("4tip-dd-exp-8thr", b4, dd8, 50L, 5000L, 8L, 40L, model = c(1L,0L,0L), link = 1L)
r10 <- block("12tip-dd-exp-4thr", b12, dd8, 50L, 5000L, 4L, 30L, model = c(1L,0L,0L), link = 1L)

## em_cpp: does the M-step see M trees, and is fhat off by log(M/N)?
lb <- c(0.01, 0, 0, 0, 0.001, 0, 0, 0); ub <- c(5, 0, 0, 0, 5, 0, 0, 0)
em1 <- function(nt, brts) {
  r <- emphasis:::em_cpp(brts = brts, init_pars = cr8, sample_size = 50L, maxN = 5000L,
                         max_missing = 1000L, max_lambda = 1e6, lower_bound = lb, upper_bound = ub,
                         xtol_rel = 1e-3, num_threads = nt, copy_trees = TRUE, model = c(0L,0L,0L), link = 0L,
                         rho = 1, rconditional = NULL)
  c(trees = r$trees, ntrees_list = length(r$trees_list %||% r$augmented_trees %||% list()), fhat = r$fhat,
    lam = r$estimates[1], mu = r$estimates[5])
}
`%||%` <- function(a, b) if (is.null(a)) b else a
e1 <- t(replicate(15, em1(1L, b12))); e4 <- t(replicate(30, em1(4L, b12)))
cat(sprintf("\n[em_cpp 12tip] 1thr: trees range %d..%d, fhat mean %.3f sd %.3f, lambda mean %.3f, mu mean %.3f\n",
            min(e1[,"trees"]), max(e1[,"trees"]), mean(e1[,"fhat"]), sd(e1[,"fhat"]), mean(e1[,"lam"]), mean(e1[,"mu"])))
ov <- e4[, "trees"] != 50
cat(sprintf("[em_cpp 12tip] 4thr: overflow %d/%d, trees range %d..%d; fhat(overflow) - mean fhat(1thr) range %.3f..%.3f vs log(trees/50) %.3f..%.3f; lambda mean %.3f, mu mean %.3f\n",
            sum(ov), nrow(e4), min(e4[,"trees"]), max(e4[,"trees"]),
            if (any(ov)) min(e4[ov,"fhat"]) - mean(e1[,"fhat"]) else NA,
            if (any(ov)) max(e4[ov,"fhat"]) - mean(e1[,"fhat"]) else NA,
            if (any(ov)) min(log(e4[ov,"trees"]/50)) else NA, if (any(ov)) max(log(e4[ov,"trees"]/50)) else NA,
            mean(e4[,"lam"]), mean(e4[,"mu"])))
str(names(emphasis:::em_cpp(brts = b4, init_pars = cr8, sample_size = 5L, maxN = 50L, max_missing = 1000L, max_lambda = 1e6,
                            lower_bound = lb, upper_bound = ub, xtol_rel = 1e-3, num_threads = 1L, copy_trees = TRUE,
                            model = c(0L,0L,0L), link = 0L, rho = 1, rconditional = NULL)))

## pipeline level: default sampler (bdi) with num_threads = 4 must be untouched; dynamic_fresh must show inflated loglik
fit_one <- function(sampling, nt) {
  f <- estimate_rates(tr12, method = "mcem", model = "cr", init_pars = c(0.5, 0.1),
                      control = list(sampling = sampling, num_threads = nt, num_trees = 50L, maxN = 5000L,
                                     max_iter = 4L, patience = 1L, tol = 1e-6,
                                     lower_bound = c(0.01, 0.001), upper_bound = c(5, 5), max_time = 120))
  c(loglik = f$loglik, tab_num_trees = max(f$details$mcem$num_trees))
}
p_bdi1 <- t(replicate(3, fit_one("bdi", 1L)));  p_bdi4 <- t(replicate(3, fit_one("bdi", 4L)))
p_df1  <- t(replicate(3, fit_one("dynamic_fresh", 1L))); p_df4 <- t(replicate(6, fit_one("dynamic_fresh", 4L)))
cat(sprintf("\n[estimate_rates 12tip] bdi 1thr loglik %s | bdi 4thr loglik %s\n",
            paste(round(p_bdi1[,1],2), collapse=","), paste(round(p_bdi4[,1],2), collapse=",")))
cat(sprintf("[estimate_rates 12tip] dynamic_fresh 1thr loglik %s | dynamic_fresh 4thr loglik %s (table num_trees %s)\n",
            paste(round(p_df1[,1],2), collapse=","), paste(round(p_df4[,1],2), collapse=","),
            paste(p_df4[,2], collapse=",")))

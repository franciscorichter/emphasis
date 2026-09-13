## H1: E_step.cpp:94 acceptance `isfinite(log_w) && 0 < exp(log_w)` is an absolute
## threshold at log_w ~ -745.13 (double underflow). Tests:
##  A. time-unit rescaling of one 60-tip CR tree: lw shifts by -(n-2)*log(s),
##     so the same tree in different units crosses the threshold.
##  B. hypothesis as stated: CR trees of 150/250/400 tips at lambda=.2, mu=.05.
##  C. R-level consequence: .mcem_dynamic_fresh (thinning) vs default BDI on the big tree.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(TreeSim); library(DDD); library(ape)})
cat("exp(-745) =", exp(-745), " exp(-746) =", exp(-746), " log(.Machine$double.xmin*2^-52) =", log(.Machine$double.xmin*2^-52), "\n")
lam <- 0.2; mu <- 0.05
sim_brts <- function(n, seed = n) { set.seed(seed)
  t <- TreeSim::sim.bd.taxa(n, 1, lam, mu, complete = FALSE)[[1]]
  sort(as.numeric(ape::branching.times(t)), decreasing = TRUE) }
pars8 <- function(l, m) c(l, 0, 0, 0, m, 0, 0, 0)
aug <- function(brts, l, m, N, maxN) tryCatch(
  emphasis:::augment_trees(brts, pars8(l, m), as.integer(N), as.integer(maxN), 10000L, 1e6, 1L, c(0L,0L,0L), 0L, 1.0),
  error = function(e) e)
summ <- function(r, label) {
  if (inherits(r, "error")) { cat(sprintf("%-28s ERROR: %s\n", label, conditionMessage(r))); return(invisible(NULL)) }
  lw <- r$logf - r$logg
  cat(sprintf("%-28s kept=%d zero_w=%d overrun=%d lambda=%d  lw: min=%.1f med=%.1f max=%.1f  logf med=%.1f logg med=%.1f  t=%.0fms\n",
              label, length(lw), r$rejected_zero_weights, r$rejected_overruns, r$rejected_lambda,
              min(lw), median(lw), max(lw), median(r$logf), median(r$logg), r$time))
  invisible(lw)
}
fhat_of <- function(lw, n_zero) { m <- max(lw); log(sum(exp(lw - m))) + m - log(length(lw) + n_zero) }

cat("\n=== A. one 60-tip CR tree, time units rescaled by s (brts*s, pars/s) ===\n")
b60 <- sim_brts(60); n <- length(b60) + 1
cat("n_tips =", n, " crown age =", b60[1], "\n")
r1 <- aug(b60, lam, mu, 40, 400); lw1 <- summ(r1, "s=1")
f1 <- fhat_of(lw1, r1$rejected_zero_weights)
cat("fhat(s=1) =", round(f1, 2), "  DDD bd_loglik(cond=0,btorph=1,soc=2) =",
    round(DDD::bd_loglik(c(lam, mu, 0, 0), c(0, 0, 1, 0, 2), b60, 0), 2), "\n")
res <- list()
for (s in c(1e2, 1e4, 1e6, 1e7)) {
  r <- aug(b60 * s, lam / s, mu / s, 40, 400)
  lw <- summ(r, sprintf("s=%g (pred shift %.0f)", s, -(n - 2) * log(s)))
  if (!is.null(lw)) cat(sprintf("   observed shift of median lw: %.1f ; predicted -(n-2)log s = %.1f ; T*s = %.2g (< 5e10 sentinel)\n",
                                median(lw) - median(lw1), -(n - 2) * log(s), b60[1] * s))
}
## straddle regime: pick s so that median lw lands at ~ -745
s_str <- exp((-745 - median(lw1)) / (-(n - 2)))
cat(sprintf("\n-- straddle: s* = %.4g puts median lw at -745 --\n", s_str))
rs <- aug(b60 * s_str, lam / s_str, mu / s_str, 40, 4000); lws <- summ(rs, "s=s*")
if (!is.null(lws)) {
  f_emph <- fhat_of(lws, rs$rejected_zero_weights)
  f_pred <- f1 - (n - 2) * log(s_str)
  cat(sprintf("   fraction of completed augmentations counted zero-weight: %d/%d = %.2f\n",
              rs$rejected_zero_weights, rs$rejected_zero_weights + length(lws),
              rs$rejected_zero_weights / (rs$rejected_zero_weights + length(lws))))
  cat(sprintf("   fhat as E_step computes it (S = N + zero_w): %.2f ; fhat(s=1) shifted by -(n-2)log s*: %.2f ; gap = %.2f\n",
              f_emph, f_pred, f_emph - f_pred))
  cat(sprintf("   survivors' min lw = %.2f (all > -745.13); DDD reference shift for same rescale = %.2f\n", min(lws),
              DDD::bd_loglik(c(lam/s_str, mu/s_str, 0, 0), c(0, 0, 1, 0, 2), b60 * s_str, 0) - DDD::bd_loglik(c(lam, mu, 0, 0), c(0, 0, 1, 0, 2), b60, 0)))
  ## IS noise at s=1 for comparison
  fr <- replicate(5, { r <- aug(b60, lam, mu, 40, 400); fhat_of(r$logf - r$logg, r$rejected_zero_weights) })
  cat(sprintf("   IS noise of fhat at s=1 over 5 replicate runs: sd = %.2f (range %.2f..%.2f)\n", sd(fr), min(fr), max(fr)))
}

cat("\n=== B. hypothesis as stated: natural units, lambda=0.2 mu=0.05 ===\n")
for (nt in c(150, 250, 400)) {
  b <- sim_brts(nt); ll <- DDD::bd_loglik(c(lam, mu, 0, 0), c(0, 0, 1, 0, 2), b, 0)
  t0 <- Sys.time(); r <- aug(b, lam, mu, 10, 60); dt <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("n=%d crown=%.1f  DDD bd_loglik=%.1f  (%.1fs)\n", nt, b[1], ll, dt))
  summ(r, sprintf("  augment n=%d", nt))
  if (inherits(r, "error")) {
    ## show the surviving magnitudes via a single attempt with maxN=1? not possible; instead show that the
    ## per-tree lw at these parameters is below -745 using the smallest kept sample at lower lambda*T? -> use eval_logf on a BDI draw
    bd <- tryCatch(emphasis:::.augment_tree_bdi(b, pars8(lam, mu), c(0L,0L,0L), 3L, 10000L, 0L, 1.0), error = function(e) e)
    if (!inherits(bd, "error")) cat(sprintf("  BDI draws on same tree: lw = %s (finite, all < -745 => thinning E_step would call them zero-weight)\n",
                                            paste(round(bd$weights, 1), collapse = ", ")))
  }
  assign(sprintf("b%d", nt), b)
}

cat("\n=== C. R-level consequence on the 400-tip tree ===\n")
lb <- c(0.01, 0.0); ub <- c(1, 0.5)
t0 <- Sys.time()
fit_thin <- withCallingHandlers(
  tryCatch(estimate_rates(b400, method = "mcem", model = "cr", init_pars = c(lam, mu),
                          control = list(sampling = "dynamic_fresh", sample_size = 10L, maxN = 60L, max_iter = 2L,
                                         lower_bound = lb, upper_bound = ub, max_time = 120)),
           error = function(e) { cat("  thinning MCEM ERROR:", conditionMessage(e), "\n"); NULL }),
  warning = function(w) { cat("  thinning MCEM WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
cat(sprintf("  thinning MCEM (%.0fs): pars = %s ; loglik = %s\n", as.numeric(Sys.time() - t0, units = "secs"),
            if (is.null(fit_thin)) "NULL" else paste(round(fit_thin$pars, 4), collapse = ","),
            if (is.null(fit_thin)) "NULL" else format(fit_thin$loglik)))
t0 <- Sys.time()
fit_bdi <- tryCatch(estimate_rates(b400, method = "mcem", model = "cr", init_pars = c(lam, mu),
                                   control = list(sample_size = 10L, max_iter = 2L, lower_bound = lb, upper_bound = ub, max_time = 120)),
                    error = function(e) { cat("  BDI MCEM ERROR:", conditionMessage(e), "\n"); NULL })
cat(sprintf("  default BDI MCEM (%.0fs): pars = %s ; loglik = %s\n", as.numeric(Sys.time() - t0, units = "secs"),
            if (is.null(fit_bdi)) "NULL" else paste(round(fit_bdi$pars, 4), collapse = ","),
            if (is.null(fit_bdi)) "NULL" else format(fit_bdi$loglik)))

cat("\n=== D. what .mcem_warn_estep says about the 400-tip tree at the TRUE parameters ===\n")
lb8 <- c(0.01, 0, 0, 0, 0, 0, 0, 0); ub8 <- c(1, 0, 0, 0, 0.5, 0, 0, 0)
withCallingHandlers(emphasis:::.mcem_warn_estep(b400, pars8(lam, mu), lb8, ub8, c(0L,0L,0L), 0L),
                    warning = function(w) { cat("  WARNING TEXT:\n", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
cat("  (thinning MCEM above returned pars 0.3098,0.122 = two 'perturb toward center' steps from (0.2,0.05) with stop_reason max_iter and no warning,\n   because max_iter=2 < the 8-failure threshold)\n")

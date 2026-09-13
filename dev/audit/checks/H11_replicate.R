## H11 replication — independent variation: different tree (20 tips, seed 3), thinning sampler contrast,
## exponential link contrast, and a per-E-step "drop the -Inf trees" counterfactual.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(DDD)})
options(width = 120)
set.seed(3)
# ddmodel 1: lambda = la0 - (la0-mu0) N/K. la0=0.8, mu0=0.1, K=30 -> slope -0.7/30, lambda=0 at N=34.3
sim  <- DDD::dd_sim(pars = c(0.8, 0.1, 30), age = 10, ddmodel = 1)
brts <- sort(as.numeric(sim$brts), decreasing = TRUE)
cat(sprintf("tree: %d tips, crown age %.3f\n", length(brts) + 1L, brts[1]))
pars <- c(0.8, -0.7 / 30, 0.1, 0)
lb <- c(0.01, -1, 0.001, 0); ub <- c(5, 0, 2, 0)
mb <- c(1L, 0L, 0L)

cat("\n== R1: BDI E-steps at generating pars, max_missing 30 (zero of lambda at N=34.3, reachable up to N=50) ==\n")
for (r in 1:4) {
  e <- emphasis:::.augment_tree_bdi(brts, pars, model_bin = mb, sample_size = 200L, max_missing = 30L, link = 0L, rho = 1)
  maxN <- max(sapply(e$trees, function(t) max(t$n)))
  cat(sprintf("  rep %d: n=%d  -Inf=%d  +Inf=%d  fhat=%s  maxN=%d\n", r, length(e$trees),
              sum(e$logf == -Inf), sum(e$logf == Inf), format(round(e$fhat, 3)), maxN))
}

cat("\n== R2: same tree, THINNING sampler at the same pars: does it filter non-finite trees? ==\n")
for (r in 1:2) {
  et <- tryCatch(emphasis:::.augment_tree_internal(brts, pars = pars, model_bin = mb, sample_size = 200L,
                                                   max_missing = 30L, link = 0L, rho = 1), error = function(e) e)
  if (inherits(et, "error")) { cat("  error:", conditionMessage(et), "\n"); next }
  cat("  names:", paste(names(et), collapse = " "), "\n")
  lf <- if (!is.null(et$logf)) et$logf else NA
  cat(sprintf("  rep %d: n_trees=%d  non-finite logf=%d  fhat=%s  rejected_zero_weights=%s\n", r, length(et$trees),
              sum(!is.finite(lf)), format(et$fhat), format(et$rejected_zero_weights)))
}

cat("\n== R3: BDI dd fit from a steep init so that the zero of lambda is reachable: pars (1.2, -0.05, 0.2, 0) zero at N=24 ==\n")
fit <- estimate_rates(brts, method = "mcem", model = "dd", init_pars = c(1.2, -0.05, 0.2, 0),
                      control = list(lower_bound = lb, upper_bound = ub, sample_size = 200L, max_iter = 10L,
                                     max_missing = 30L, num_threads = 1L, sampling = "bdi", verbose = FALSE))
tr <- fit$details$mcem; tr$zero_at_N <- -tr$par1 / tr$par2
print(round(tr[, c("par1", "par2", "par5", "fhat", "delta_max", "zero_at_N")], 4))
cat(sprintf("stop_reason: %s  iterations: %d  loglik: %s\n", fit$details$stop_reason, nrow(tr), format(fit$loglik)))
if (!is.null(fit$details$final_IS)) cat(sprintf("final_IS: fhat=%s +Inf=%d -Inf=%d\n", format(fit$details$final_IS$fhat),
   sum(fit$details$final_IS$logf == Inf), sum(fit$details$final_IS$logf == -Inf)))

cat("\n== R4: freeze counterfactual at the R3 final pars: same E-step, m_cpp with (a) w_norm as bdi.R does, (b) -Inf trees dropped ==\n")
th <- as.numeric(fit$details$pars)
run_m <- function(trees, w, init) {
  es <- list(trees = trees, weights = w, rejected = 0L, rejected_overruns = 0L, rejected_lambda = 0L,
             rejected_zero_weights = 0L, time = 0, fhat = 0)
  m <- emphasis:::m_cpp(e_step = es, init_pars = init, plugin = "rpd1",
                        lower_bound = c(0.01, -1, 0, 0, 0.001, 0, 0, 0), upper_bound = c(5, 0, 0, 0, 2, 0, 0, 0),
                        xtol_rel = 1e-3, num_threads = 1L, model = mb, link = 0L, rho = 1)
  c(max(abs(as.numeric(m$estimates) - init)), m$nlopt)
}
for (r in 1:6) {
  e <- emphasis:::.augment_tree_bdi(brts, th, model_bin = mb, sample_size = 200L, max_missing = 30L, link = 0L, rho = 1)
  lw <- e$weights; w <- exp(lw - max(lw)); w <- w / sum(w) * length(w)
  a <- run_m(e$trees, w, th)
  keep <- is.finite(e$logf)
  wk <- exp(lw[keep] - max(lw[keep])); wk <- wk / sum(wk) * length(wk)
  b <- run_m(e$trees[keep], wk, th)
  cat(sprintf("  rep %d: -Inf=%2d +Inf=%d | as-is: max|delta|=%.2e (nlopt %d) | dropped: max|delta|=%.2e (nlopt %d)\n",
              r, sum(e$logf == -Inf), sum(e$logf == Inf), a[1], a[2], b[1], b[2]))
}

cat("\n== R5: exponential link (rate = exp(eta) > 0, no zero reachable), same tree, 4 iterations ==\n")
fit_e <- tryCatch(estimate_rates(brts, method = "mcem", model = "dd", link = "exponential",
                        init_pars = c(log(0.8), -0.03, log(0.1), 0),
                        control = list(sample_size = 200L, max_iter = 4L, max_missing = 30L, num_threads = 1L,
                                       sampling = "bdi", verbose = FALSE)), error = function(e) e)
if (inherits(fit_e, "error")) cat("  error:", conditionMessage(fit_e), "\n") else {
  print(round(fit_e$details$mcem[, c("par1", "par2", "par3", "fhat", "delta_max")], 4))
  cat("  stop_reason:", fit_e$details$stop_reason, " final_IS non-finite:", sum(!is.finite(fit_e$details$final_IS$logf)), "\n")
}

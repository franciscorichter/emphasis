## H20 — Is the reported MCEM loglik fhat(theta_{K-1}) (E-step before the last
##       M-step) while pars = theta_K, and how large is the lag under the
##       different stop reasons?
##
## Self-contained. num_threads = 1. rho = 1, cond = NULL throughout.
##
## Parts
##  A  CR, BDI sampler (the default path). Converged fit: read theta_{K-1},
##     theta_K from the trace; recompute fhat at both with N = 500 BDI trees;
##     compare with fit$loglik and with DDD::bd_loglik at both points.
##     (BDI under CR is zero-variance IS, so fhat(theta) is the exact
##     Nee likelihood; every difference is lag, not MC noise.)
##  B  CR, BDI, forced early stops: max_iter = 1, 2, 5 from the midpoint init.
##     Lag = fhat(theta_K) - fit$loglik; also AIC consequences.
##  C  DD (linear, N-only), BDI: same with DDD::dd_loglik.
##  D  CR, thinning sampler (dynamic_fresh): one extra em_cpp at theta_K with
##     N = 2000 (the hypothesis' literal test) and at theta_{K-1} with N = 2000
##     to separate lag from MC noise.
##  E  Does the pipeline's stage selection depend on loglik at all? (code read)

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
ns <- asNamespace("emphasis")
set.seed(20)

## ---- data: one fixed 30-tip CR tree ----------------------------------------
tr   <- ape::rphylo(30, 0.8, 0.3)
brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
cat(sprintf("tree: %d tips, crown age %.3f\n", length(brts) + 1, brts[1]))

bd_ll <- function(p) DDD::bd_loglik(c(p[1], p[2], 0, 0), c(0, 0, 0, 0, 2), brts, 0)
ml <- DDD::bd_ML(brts, initparsopt = c(0.8, 0.3), idparsopt = 1:2, cond = 0, soc = 2,
                 btorph = 0, verbose = FALSE)
cat(sprintf("DDD::bd_ML (cond = 0): lam = %.4f mu = %.4f loglik = %.4f\n",
            ml$lambda0, ml$mu0, ml$loglik))

lb <- c(0.01, 0); ub <- c(3, 3)
fhat_bdi <- function(p8, model_bin = c(0L, 0L, 0L), N = 500L)
  ns$.augment_tree_bdi(brts, pars = p8, model_bin = model_bin,
                       sample_size = N, max_missing = 1e4, link = 0L, rho = 1)$fhat

report_fit <- function(fit, label, model_bin = c(0L, 0L, 0L), exact = NULL, N = 500L) {
  m  <- fit$details$mcem
  K  <- nrow(m)
  th_K   <- as.numeric(fit$details$pars)
  th_Km1 <- if (K >= 2) as.numeric(m[K - 1, seq_along(th_K)]) else NA
  f_K    <- fhat_bdi(th_K, model_bin, N)
  f_Km1  <- if (K >= 2) fhat_bdi(th_Km1, model_bin, N) else NA
  cat(sprintf("[%s] stop=%s K=%d  pars(theta_K)=(%s)\n", label, fit$details$stop_reason, K,
              paste(sprintf("%.4f", fit$pars), collapse = ", ")))
  cat(sprintf("   fit$loglik            = %.4f  (trace fhat[K] = %.4f)\n", fit$loglik, m$fhat[K]))
  cat(sprintf("   fhat(theta_K)   N=%d  = %.4f   lag = fhat(theta_K) - fit$loglik = %+.4f\n",
              N, f_K, f_K - fit$loglik))
  if (K >= 2)
    cat(sprintf("   fhat(theta_K-1) N=%d  = %.4f   (fit$loglik - fhat(theta_K-1) = %+.4f: MC noise only)\n",
                N, f_Km1, fit$loglik - f_Km1))
  if (!is.null(exact)) {
    cat(sprintf("   exact(theta_K) = %.4f  exact(theta_K-1) = %s\n", exact(th_K),
                if (K >= 2) sprintf("%.4f", exact(th_Km1)) else "NA"))
  }
  cat(sprintf("   AIC reported = %.4f  AIC at fhat(theta_K) = %.4f\n",
              fit$AIC, -2 * f_K + 2 * fit$n_pars))
  invisible(list(lag = f_K - fit$loglik, K = K, stop = fit$details$stop_reason,
                 f_K = f_K, f_Km1 = f_Km1, th_K = th_K, th_Km1 = th_Km1))
}

## ---- Part A: CR, BDI, converged -------------------------------------------
cat("\n==== Part A: CR, BDI, default stopping (converged) ====\n")
lags <- numeric(0)
for (r in 1:3) {
  fit <- estimate_rates(brts, model = "cr", method = "mcem",
                        control = list(lower_bound = lb, upper_bound = ub,
                                       sample_size = 50L, max_iter = 200L,
                                       max_time = 60, num_threads = 1L))
  rp <- report_fit(fit, sprintf("A rep %d", r), exact = bd_ll)
  lags <- c(lags, rp$lag)
}
cat(sprintf("converged-stop lags (fhat(theta_K) - fit$loglik): %s\n",
            paste(sprintf("%+.4f", lags), collapse = " ")))

## ---- Part B: CR, BDI, forced early stops -----------------------------------
cat("\n==== Part B: CR, BDI, max_iter stops from the midpoint init ====\n")
for (mi in c(1L, 2L, 5L)) {
  fit <- estimate_rates(brts, model = "cr", method = "mcem",
                        control = list(lower_bound = lb, upper_bound = ub,
                                       sample_size = 50L, max_iter = mi,
                                       max_time = 60, num_threads = 1L))
  report_fit(fit, sprintf("B max_iter=%d", mi), exact = bd_ll)
}
cat("init (midpoint) =", (lb + ub) / 2, " exact loglik there =", round(bd_ll((lb + ub) / 2), 4), "\n")

## ---- Part C: DD linear, BDI ------------------------------------------------
cat("\n==== Part C: DD (linear, N-only), BDI ====\n")
dd_ll <- function(p) {
  # compact p = (b0, bN, g0, gN); DDD ddmodel=1: lambda = lam0 - (lam0-mu0)N/K
  # emphasis linear dd is lambda = b0 + bN*N, mu = g0 + gN*N. Map only when gN = 0:
  # lam0 = b0, mu0 = g0, K = -b0/bN... DDD model 1 has mu constant; use it if gN == 0.
  if (abs(p[4]) > 1e-12) return(NA_real_)
  K <- -(p[1] - p[3]) / p[2]
  if (!is.finite(K) || K <= 0) return(NA_real_)
  DDD::dd_loglik(c(p[1], p[3], K), c(200, 1, 0, 0, 0, 2), brts, 0)
}
lb_dd <- c(0.01, -0.5, 0, 0); ub_dd <- c(3, 0, 3, 0)   # gammaN fixed at 0
for (mi in c(200L, 1L, 2L)) {
  fit <- estimate_rates(brts, model = "dd", method = "mcem",
                        control = list(lower_bound = lb_dd, upper_bound = ub_dd,
                                       sample_size = 50L, max_iter = mi,
                                       max_time = 90, num_threads = 1L))
  report_fit(fit, sprintf("C dd max_iter=%d", mi), model_bin = c(1L, 0L, 0L),
             exact = dd_ll, N = 200L)
}

## ---- Part D: CR, thinning sampler, extra em_cpp at theta_K with N = 2000 ---
cat("\n==== Part D: CR, thinning (dynamic_fresh), em_cpp N=2000 at theta_K ====\n")
fhat_thin <- function(p8, N = 2000L) {
  r <- ns$em_cpp(brts = brts, init_pars = p8, sample_size = N, maxN = 20000L,
                 max_missing = 1e4, max_lambda = 1e6, lower_bound = c(lb[1], 0, 0, 0, lb[2], 0, 0, 0),
                 upper_bound = c(ub[1], 0, 0, 0, ub[2], 0, 0, 0), xtol_rel = 1e-3,
                 num_threads = 1L, copy_trees = FALSE, model = c(0L, 0L, 0L), link = 0L,
                 rho = 1, rconditional = NULL)
  c(fhat = r$fhat, ess = ns$.ess_from_lw(r$logf - r$logg))
}
for (mi in c(200L, 2L)) {
  fit <- estimate_rates(brts, model = "cr", method = "mcem",
                        control = list(lower_bound = lb, upper_bound = ub,
                                       sampling = "dynamic_fresh",
                                       sample_size = 200L, maxN = 5000L, max_iter = mi,
                                       max_time = 90, num_threads = 1L))
  m <- fit$details$mcem; K <- nrow(m)
  th_K <- as.numeric(fit$details$pars); th_Km1 <- as.numeric(m[K - 1, 1:8])
  a <- replicate(3, fhat_thin(th_K)); b <- replicate(3, fhat_thin(th_Km1))
  cat(sprintf("[D max_iter=%d] stop=%s K=%d pars=(%.4f, %.4f) fit$loglik=%.4f (N=200)\n",
              mi, fit$details$stop_reason, K, fit$pars[1], fit$pars[2], fit$loglik))
  cat(sprintf("   fhat(theta_K)   N=2000 x3: %s  (ESS %s)   exact bd_loglik(theta_K)   = %.4f\n",
              paste(sprintf("%.4f", a["fhat", ]), collapse = " "),
              paste(sprintf("%.0f", a["ess", ]), collapse = "/"), bd_ll(th_K[c(1, 5)])))
  cat(sprintf("   fhat(theta_K-1) N=2000 x3: %s  (ESS %s)   exact bd_loglik(theta_K-1) = %.4f\n",
              paste(sprintf("%.4f", b["fhat", ]), collapse = " "),
              paste(sprintf("%.0f", b["ess", ]), collapse = "/"), bd_ll(th_Km1[c(1, 5)])))
  cat(sprintf("   lag (mean fhat(theta_K) - fit$loglik) = %+.4f ; exact lag = %+.4f\n",
              mean(a["fhat", ]) - fit$loglik, bd_ll(th_K[c(1, 5)]) - bd_ll(th_Km1[c(1, 5)])))
}

## ---- Part E: pipeline stage selection --------------------------------------
cat("\n==== Part E: does the pipeline rank stages by loglik? ====\n")
cat("see R/pipeline.R:289-302 — first finite among mcem > cem > gam; no loglik comparison.\n")

## H20 replication — independent re-test of "fit$loglik = fhat(theta_{K-1})
## next to pars = theta_K".
##
## Varies from the verifier's H20.R: different tree (25 tips, lam=1.0 mu=0.6,
## seed 7 -> MLE with mu > 0), whole-trace identity check (every row k, not
## only K-1), the time_budget stop path, a log-link CR fit (forces the thinning
## sampler), and the exposure of stop_reason in the fit object / print method.
## num_threads = 1, rho = 1, cond = NULL.

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
ns <- asNamespace("emphasis")
set.seed(7)

tr   <- ape::rphylo(25, 1.0, 0.6)
brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
cat(sprintf("tree: %d tips, crown age %.3f\n", length(brts) + 1, brts[1]))
bd_ll <- function(p) DDD::bd_loglik(c(p[1], p[2], 0, 0), c(0, 0, 0, 0, 2), brts, 0)
ml <- DDD::bd_ML(brts, initparsopt = c(1, 0.5), idparsopt = 1:2, cond = 0, soc = 2,
                 btorph = 0, verbose = FALSE)
cat(sprintf("DDD::bd_ML: lam = %.4f mu = %.4f loglik = %.4f\n", ml$lambda0, ml$mu0, ml$loglik))

lb <- c(0.01, 0); ub <- c(4, 3)
fhat_bdi <- function(p8, N = 500L)
  ns$.augment_tree_bdi(brts, pars = p8, model_bin = c(0L, 0L, 0L),
                       sample_size = N, max_missing = 1e4, link = 0L, rho = 1)$fhat
init8 <- ns$.expand_pars((lb + ub) / 2, c(0L, 0L, 0L))
## ---- 4. log link CR (thinning sampler forced) --------------------------------
cat("\n==== 4. CR, link = 'log' (thinning sampler), converged + max_iter=2 ====\n")
lb_log <- c(-3, -3); ub_log <- c(1.5, 1.5)
fhat_thin <- function(p8, N = 1500L) {
  r <- ns$em_cpp(brts = brts, init_pars = p8, sample_size = N, maxN = 20000L,
                 max_missing = 1e4, max_lambda = 1e6,
                 lower_bound = ns$.expand_pars(lb_log, c(0L, 0L, 0L)),
                 upper_bound = ns$.expand_pars(ub_log, c(0L, 0L, 0L)), xtol_rel = 1e-3,
                 num_threads = 1L, copy_trees = FALSE, model = c(0L, 0L, 0L), link = 1L,
                 rho = 1, rconditional = NULL)
  c(fhat = r$fhat, ess = ns$.ess_from_lw(r$logf - r$logg))
}
for (mi in c(200L, 2L)) {
  fit <- tryCatch(estimate_rates(brts, model = "cr", method = "mcem", link = "exponential",
                        control = list(lower_bound = lb_log, upper_bound = ub_log, sampling = "dynamic_fresh",
                                       sample_size = 200L, maxN = 5000L, max_iter = mi,
                                       max_time = 90, num_threads = 1L)),
                  error = function(e) { cat("  error:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) next
  m <- fit$details$mcem; K <- nrow(m)
  cat(sprintf("  sampling used: %s\n", if (is.null(fit$details$sampling)) "(not recorded)" else fit$details$sampling))
  th_K <- as.numeric(fit$details$pars)
  th_Km1 <- if (K >= 2) as.numeric(m[K - 1, 1:8]) else ns$.expand_pars((lb_log + ub_log) / 2, c(0L, 0L, 0L))
  a <- replicate(3, fhat_thin(th_K)); b <- replicate(3, fhat_thin(th_Km1))
  cat(sprintf("  max_iter=%d: stop=%s K=%d pars(log)=(%.4f, %.4f) fit$loglik=%.4f\n",
              mi, fit$details$stop_reason, K, fit$pars[1], fit$pars[2], fit$loglik))
  cat(sprintf("    fhat(theta_K)   N=1500 x3: %s (ESS %s)  exact=%.4f\n",
              paste(sprintf("%.4f", a["fhat", ]), collapse = " "),
              paste(sprintf("%.0f", a["ess", ]), collapse = "/"), bd_ll(exp(th_K[c(1, 5)]))))
  cat(sprintf("    fhat(theta_K-1) N=1500 x3: %s (ESS %s)  exact=%.4f\n",
              paste(sprintf("%.4f", b["fhat", ]), collapse = " "),
              paste(sprintf("%.0f", b["ess", ]), collapse = "/"), bd_ll(exp(th_Km1[c(1, 5)]))))
  cat(sprintf("    lag (mean fhat(theta_K) - fit$loglik) = %+.4f ; fit$loglik - mean fhat(theta_K-1) = %+.4f ; exact lag = %+.4f\n",
              mean(a["fhat", ]) - fit$loglik, fit$loglik - mean(b["fhat", ]),
              bd_ll(exp(th_K[c(1, 5)])) - bd_ll(exp(th_Km1[c(1, 5)]))))
}

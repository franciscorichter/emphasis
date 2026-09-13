## H12 (part 2) — can BDI-MCEM reach an MLE with mu > lam?
## Tree: crown age 10, all other branching times in [0.2, 1.0] (extreme
## pull of the present) so that DDD::bd_ML (cond = 0) returns mu0 > lam0.
## Run estimate_rates(method = "mcem") (default sampling = "bdi") from an
## init with mu < lam, bounds spanning mu > lam, verbose = TRUE to count
## E-step failures; then the same with the one-line fix applied in-session;
## then the thinning sampler for comparison.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({library(emphasis); library(DDD)})
set.seed(1212)
ns <- asNamespace("emphasis")

## search a few pull-of-the-present trees until bd_ML gives mu0 > lam0
for (w in c(0.3, 0.2, 0.1, 0.05)) {
  brts <- c(10, sort(runif(13, 0.02, w), decreasing = TRUE))
  ml <- suppressMessages(DDD::bd_ML(brts, cond = 0, btorph = 0, soc = 2, verbose = FALSE))
  cat(sprintf("window [0.02,%.2f]: bd_ML lambda0=%.4f mu0=%.4f\n", w, ml$lambda0, ml$mu0))
  if (ml$mu0 > ml$lambda0) break
}
cat("brts:", round(brts, 3), "\n")
cat(sprintf("DDD::bd_ML cond=0: lambda0=%.4f mu0=%.4f loglik=%.4f\n", ml$lambda0, ml$mu0, ml$loglik))
lb <- c(0, 0); ub <- c(50, 50); init <- c(2, 1)
cat("bd_loglik at (2,1):", DDD::bd_loglik(c(2, 1, 0, 0), c(0, 0, 0, 0, 2), brts, 0),
    " at MLE:", DDD::bd_loglik(c(ml$lambda0, ml$mu0, 0, 0), c(0, 0, 0, 0, 2), brts, 0), "\n")

run <- function(label, sampling = "bdi", ss = 20L, iters = 40L) {
  msgs <- character(0)
  t0 <- proc.time()[3]
  fit <- withCallingHandlers(
    estimate_rates(brts, model = "cr", method = "mcem", init_pars = init,
                   control = list(lower_bound = lb, upper_bound = ub, sampling = sampling,
                                  sample_size = ss, max_iter = iters, max_time = 120,
                                  num_threads = 1L, verbose = TRUE)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  el <- proc.time()[3] - t0
  n_fail <- sum(grepl("E-step failed", msgs))
  cat(sprintf("[%s] %.1fs  pars=(%.4f, %.4f)  stop=%s  iters=%s  E-step failures=%d  loglik=%.4f\n",
              label, el, fit$pars[1], fit$pars[2], fit$details$stop_reason,
              format(fit$details$iterations), n_fail, fit$loglik))
  m <- fit$details$mcem
  if (!is.null(m)) {
    cat("  trajectory (lam, mu) every 4th recorded iteration:\n")
    idx <- unique(c(seq(1, nrow(m), by = 4), nrow(m)))
    print(round(cbind(lam = m$par1[idx], mu = m$par5[idx], fhat = m$fhat[idx]), 4))
    cat(sprintf("  max mu/lam visited among recorded iterations: %.4f\n", max(m$par5 / m$par1)))
  }
  fit
}

cat("\n--- unfixed BDI ---\n")
f_bdi <- run("bdi-unfixed")

cat("\n--- apply fix (abs guard) in-session ---\n")
fixed <- function(t1, t2, n, k, lam0, mu0, tp) {
  if (t2 - t1 < 1e-15) return(0)
  d <- lam0 - mu0
  if (abs(d) < 1e-15) {
    a1 <- 1 + lam0 * (tp - t1); a2 <- 1 + lam0 * (tp - t2)
    I_lam <- log(a1 / a2)
    I_mu  <- mu0 * (t2 - t1) + mu0 * lam0 * ((tp - t1)^2 - (tp - t2)^2) / 2  # corrected critical I_mu
    return((n + 2L * k) * I_lam + n * I_mu)
  }
  E1 <- exp(-d * (tp - t1)); E2 <- exp(-d * (tp - t2))
  I_lam <- mu0 * (t2 - t1) + log((lam0 - mu0 * E2) / (lam0 - mu0 * E1))
  if (n > 0L) {
    oE1 <- 1 - E1; oE2 <- 1 - E2
    if (abs(oE2) < 1e-300) return(Inf)
    I_mu <- lam0 * (t2 - t1) + log(oE1 / oE2)
  } else I_mu <- 0
  (n + 2L * k) * I_lam + n * I_mu
}
environment(fixed) <- ns
assignInNamespace(".bdi_integral_cr", fixed, ns)
f_fix <- run("bdi-fixed")


cat("\n--- fixed BDI: fhat at the bd_ML MLE and at (2,1) vs DDD::bd_loglik (offset check) ---\n")
for (p in list(c(ml$lambda0, ml$mu0), c(2, 1), c(1, 2))) {
  a <- ns$.augment_tree_bdi(brts, pars = p, sample_size = 20L)
  cat(sprintf("pars=(%.3f,%.3f): fhat=%.4f  bd_loglik=%.4f  diff=%.4f  sd(lw)=%.1e\n",
              p[1], p[2], a$fhat, DDD::bd_loglik(c(p, 0, 0), c(0, 0, 0, 0, 2), brts, 0),
              a$fhat - DDD::bd_loglik(c(p, 0, 0), c(0, 0, 0, 0, 2), brts, 0), sd(a$weights)))
}

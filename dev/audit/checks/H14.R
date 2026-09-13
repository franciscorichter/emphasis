## H14: loglik_var bootstraps the K=2 moment-corrected estimator
##        ell_BC = mean(lw) + log(1 + m2/2)          (.is_fhat bias_correct=TRUE, no n_zero)
## while the reported loglik is the log-mean-exp
##        ell_LME = log(sum exp(lw)) - log(N + n_zero) (E_step.cpp:143, .is_fhat default)
## and print.emphasis_fit / compare_models present sqrt(loglik_var) as the MC s.e. of loglik.
##
## Tests
##  A. Synthetic lw (seeded, pure R): for increasing weight-tail heaviness compare
##     V_pkg  = .bootstrap_fhat_var(logf, logg)      [what the package reports]
##     V_lme  = bootstrap variance of ell_LME on the same lw (what it should report)
##     V_true = variance of ell_LME over fresh independent draws (ground truth for the MC s.e.)
##     and the point values ell_BC vs ell_LME.
##  B. Real E-steps: augment_trees on a 30-tip CR tree, model cr and dd, N = 200,
##     several parameter values (some far from the MLE -> low ESS); same three quantities,
##     plus a case with rejected_zero_weights > 0 to show the bootstrap drops the n_zero term.
##  C. compare_models: two fits sharing one tree; pairwise p from V_pkg vs V_lme.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(TreeSim); library(ape); library(DDD)})
options(width = 150)

lme   <- function(lw, n_zero = 0L) { m <- max(lw); log(sum(exp(lw - m))) + m - log(length(lw) + n_zero) }
bcK2  <- function(lw) { e <- mean(lw); e + log(1 + mean((lw - e)^2) / 2) }
ess   <- function(lw) { w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }
boot_var <- function(lw, f, B = 200L) {
  n <- length(lw)
  stats::var(vapply(seq_len(B), function(.) f(lw[sample.int(n, n, replace = TRUE)]), 0))
}
pkg_var <- function(lw) emphasis:::.bootstrap_fhat_var(lw, rep(0, length(lw)), K = 2L, B = 200L)

## ---------------------------------------------------------------------------
cat("=== A. synthetic lw = mu + sigma*Z  (lognormal weights), N = 200, B = 200 ===\n")
cat("   ell_true = log E[w] = mu + sigma^2/2.  V_true from 2000 fresh draws of ell_LME.\n")
set.seed(14)
N <- 200L
res_A <- do.call(rbind, lapply(c(0.5, 1, 2, 3, 4), function(sig) {
  draw <- function() rnorm(N, 0, sig)
  lw <- draw()
  v_pkg <- pkg_var(lw)
  v_lme <- boot_var(lw, lme)
  v_true_lme <- stats::var(replicate(2000, lme(draw())))
  v_true_bc  <- stats::var(replicate(2000, bcK2(draw())))
  data.frame(sigma = sig, ESS = round(ess(lw), 1),
             ell_true = sig^2 / 2, ell_LME = round(lme(lw), 3), ell_BC = round(bcK2(lw), 3),
             se_pkg = round(sqrt(v_pkg), 4), se_boot_LME = round(sqrt(v_lme), 4),
             se_true_LME = round(sqrt(v_true_lme), 4), se_true_BC = round(sqrt(v_true_bc), 4),
             ratio_pkg_over_trueLME = round(sqrt(v_pkg / v_true_lme), 2))
}))
print(res_A, row.names = FALSE)
cat("   (se_pkg is what print.emphasis_fit shows as 'MC se'; se_true_LME is the MC se of the reported loglik)\n")
lw5 <- rnorm(5); cat(sprintf("   .is_fhat(bias_correct=TRUE) = %.6f ; mean+log(1+m2/2) = %.6f ; lme = %.6f\n",
    emphasis:::.is_fhat(lw5, rep(0,5), bias_correct = TRUE, K = 2L), bcK2(lw5), lme(lw5)))

## ---------------------------------------------------------------------------
cat("\n=== B. real E-steps on a 30-tip CR tree (TreeSim, lambda=.5, mu=.2), N = 200 trees ===\n")
set.seed(30)
phy  <- TreeSim::sim.bd.taxa(30, 1, 0.5, 0.2, complete = FALSE)[[1]]
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
cat(sprintf("crown age %.3f, %d tips\n", brts[1], ape::Ntip(phy)))
# closed-form reference for cr: DDD::bd_loglik (cond 0, soc 2, labelled-density convention differs
# by a constant; only used to sanity-check the level, not the variance)
aug <- function(p8, model_bin, N = 200L, maxN = 20000L) tryCatch(
  emphasis:::augment_trees(brts, as.numeric(p8), as.integer(N), as.integer(maxN), 10000L, 1e6, 1L,
                           as.integer(model_bin), 0L, 1.0),
  error = function(e) e)
cases <- list(
  list(lab = "cr  lam=.5 mu=.2  (near truth)",  p8 = c(0.5, 0, 0, 0, 0.2, 0, 0, 0), mb = c(0,0,0)),
  list(lab = "cr  lam=1.2 mu=1.0 (very high)",   p8 = c(1.2, 0, 0, 0, 1.0, 0, 0, 0), mb = c(0,0,0)),
  list(lab = "dd  b0=.8 bN=-.01 mu=.3",          p8 = c(0.8, -0.01, 0, 0, 0.3, 0, 0, 0), mb = c(1,0,0)),
  list(lab = "dd  b0=1.5 bN=-.045 mu=1.0",       p8 = c(1.5, -0.045, 0, 0, 1.0, 0, 0, 0), mb = c(1,0,0))
)
res_B <- do.call(rbind, lapply(cases, function(cs) {
  r <- aug(cs$p8, cs$mb)
  if (inherits(r, "error")) { cat(sprintf("%-32s ERROR: %s\n", cs$lab, conditionMessage(r))); return(NULL) }
  lw <- r$logf - r$logg; nz <- r$rejected_zero_weights
  # replicate E-steps for the ground-truth MC se of the *reported* estimator (fresh trees, same theta)
  reps <- replicate(8, { rr <- aug(cs$p8, cs$mb); if (inherits(rr, "error")) NA else lme(rr$logf - rr$logg, rr$rejected_zero_weights) })
  data.frame(case = cs$lab, ESS = round(ess(lw), 1), n_zero = nz, overrun = r$rejected_overruns, lam_rej = r$rejected_lambda,
             ell_LME = round(lme(lw, nz), 3), ell_BC = round(bcK2(lw), 3),
             se_pkg = round(sqrt(pkg_var(lw)), 3), se_boot_LME = round(sqrt(boot_var(lw, function(x) lme(x, nz))), 3),
             se_8reps_LME = round(stats::sd(reps, na.rm = TRUE), 3), stringsAsFactors = FALSE)
}))
print(res_B, row.names = FALSE)

## ---------------------------------------------------------------------------
cat("\n=== B2. zero-weight completions: does the bootstrap see n_zero? ===\n")
cat("   .bootstrap_fhat_var has no n_zero argument; with n_zero>0 .is_fhat(bias_correct=TRUE) would fall back to LME,\n")
cat("   but the bootstrap never passes it. Simulate: N=200 valid lw plus n_zero zero-weight completions.\n")
set.seed(141)
lw <- rnorm(200, 0, 1.5)
for (nz in c(0L, 50L, 200L, 800L)) {
  # true MC variance of the reported estimator under the inverse-binomial design:
  # keep drawing until 200 accepted; each attempt is zero-weight w.p. p = nz/(200+nz)
  p <- nz / (200 + nz)
  v_true <- stats::var(replicate(2000, { k <- if (p > 0) rnbinom(1, size = 200, prob = 1 - p) else 0L
                                          lme(rnorm(200, 0, 1.5), k) }))
  cat(sprintf("   n_zero=%4d  ell_LME=%8.3f  ell_BC(bootstrapped)=%8.3f  se_pkg=%.4f  se_boot_LME(with nz)=%.4f  se_true=%.4f\n",
              nz, lme(lw, nz), bcK2(lw), sqrt(pkg_var(lw)), sqrt(boot_var(lw, function(x) lme(x, nz))), sqrt(v_true)))
}

## ---------------------------------------------------------------------------
cat("\n=== C. compare_models pairwise p-value with V_pkg vs V_lme ===\n")
mk_fit <- function(ll, v, np, model) structure(list(pars = rep(0, np), loglik = ll, loglik_var = v,
                                                     AIC = -2 * ll + 2 * np, n_pars = np, model = model,
                                                     method = "mcem", cond = FALSE), class = "emphasis_fit")
if (!is.null(res_B) && nrow(res_B) >= 2) {
  # take the cr near-truth E-step and the heaviest dd E-step as the two 'fits'
  i1 <- 1L; i2 <- nrow(res_B)
  f_pkg <- emphasis:::compare_models(CR = mk_fit(res_B$ell_LME[i1], res_B$se_pkg[i1]^2, 2L, c(0L,0L,0L)),
                          DD = mk_fit(res_B$ell_LME[i2], res_B$se_pkg[i2]^2, 4L, c(1L,0L,0L)))
  f_lme <- emphasis:::compare_models(CR = mk_fit(res_B$ell_LME[i1], res_B$se_boot_LME[i1]^2, 2L, c(0L,0L,0L)),
                          DD = mk_fit(res_B$ell_LME[i2], res_B$se_boot_LME[i2]^2, 4L, c(1L,0L,0L)))
  f_tru <- emphasis:::compare_models(CR = mk_fit(res_B$ell_LME[i1], res_B$se_8reps_LME[i1]^2, 2L, c(0L,0L,0L)),
                          DD = mk_fit(res_B$ell_LME[i2], res_B$se_8reps_LME[i2]^2, 4L, c(1L,0L,0L)))
  cat(sprintf("   dAIC = %.3f\n", diff(range(f_pkg$AIC))))
  cat(sprintf("   p (package V_pkg)      = %.4g\n", attr(f_pkg, "pairwise_p")[1, 2]))
  cat(sprintf("   p (bootstrap of LME)   = %.4g\n", attr(f_lme, "pairwise_p")[1, 2]))
  cat(sprintf("   p (8 fresh E-steps)   = %.4g\n", attr(f_tru, "pairwise_p")[1, 2]))
}

## ---------------------------------------------------------------------------
cat("\n=== D. end-to-end: estimate_rates cr on the tree, print() output ===\n")
fit <- tryCatch(estimate_rates(phy, model = "cr", method = "mcem", control = list(lower_bound = c(0, 0), upper_bound = c(2, 1), max_iter = 6L, num_trees = 100L, num_threads = 1L, max_time = 120)),
                error = function(e) e)
if (inherits(fit, "error")) cat("estimate_rates error:", conditionMessage(fit), "\n") else {
  print(fit)
  d <- fit$details
  fi <- if (!is.null(d$final_IS)) d$final_IS else d$details$final_IS
  if (!is.null(fi)) {
    lw <- fi$lw
    cat(sprintf("   final_IS: n=%d ESS=%.1f n_zero=%d  fhat=%.3f  reported loglik=%.3f  ell_BC=%.3f\n",
                length(lw), fi$ESS, fi$rejected_zero_weights, fi$fhat, fit$loglik, bcK2(lw)))
    cat(sprintf("   sqrt(loglik_var)=%.4f  recomputed pkg=%.4f  boot LME=%.4f\n",
                sqrt(fit$loglik_var), sqrt(pkg_var(lw)), sqrt(boot_var(lw, function(x) lme(x, fi$rejected_zero_weights)))))
  } else cat("   (no final_IS in details; names:", paste(names(d), collapse = ","), ")\n")
}

## ---------------------------------------------------------------------------
cat("\n=== C2. compare_models p-value at a small dAIC, using the dd b0=1.5 row's three s.e. values ===\n")
if (!is.null(res_B)) {
  i2 <- nrow(res_B); i1 <- 1L
  for (dA in c(1, 2, 4)) {
    pp <- sapply(c("se_pkg", "se_boot_LME", "se_8reps_LME"), function(col) {
      f <- emphasis:::compare_models(A = mk_fit(-60, res_B[[col]][i1]^2, 2L, c(0L,0L,0L)),
                                     B = mk_fit(-60 - dA/2 + 1, res_B[[col]][i2]^2, 4L, c(1L,0L,0L)))
      attr(f, "pairwise_p")[1, 2]
    })
    cat(sprintf("   dAIC=%.0f  p(pkg)=%.3f  p(boot LME)=%.3f  p(fresh E-steps)=%.3f\n", dA, pp[1], pp[2], pp[3]))
  }
}

## ---------------------------------------------------------------------------
cat("\n=== E. end-to-end default path (BDI sampler), model dd, rho=1, cond=NULL ===\n")
fit_dd <- tryCatch(estimate_rates(phy, model = "dd", method = "mcem",
                                  control = list(lower_bound = c(0, -0.1, 0, -0.01), upper_bound = c(2, 0.01, 1, 0.01),
                                                 max_iter = 6L, num_trees = 100L, num_threads = 1L, max_time = 150)),
                   error = function(e) e)
if (inherits(fit_dd, "error")) cat("estimate_rates error:", conditionMessage(fit_dd), "\n") else {
  print(fit_dd)
  fi <- fit_dd$details$final_IS
  if (!is.null(fi)) {
    lw <- fi$lw
    cat(sprintf("   final_IS: n=%d ESS=%.1f n_zero=%d  fhat=%.3f  reported loglik=%.3f  ell_BC(bootstrapped stat)=%.3f\n",
                length(lw), fi$ESS, .subset2(fi, "rejected_zero_weights"), fi$fhat, fit_dd$loglik, bcK2(lw)))
    cat(sprintf("   sqrt(loglik_var)=%.4f  recomputed pkg=%.4f  boot LME=%.4f\n",
                sqrt(fit_dd$loglik_var), sqrt(pkg_var(lw)), sqrt(boot_var(lw, lme))))
  } else cat("   (no final_IS; names:", paste(names(fit_dd$details), collapse = ","), ")\n")
}

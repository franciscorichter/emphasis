## H13 — The GAM stage reports as `loglik` the additive smooth's predicted value
## at the L-BFGS-B optimum (gam.R find_MLE: loglik = -opt$value; inference.R
## .run_gam passes it up; estimate_rates builds AIC from it). Hypothesis: an
## additive Σ s(p_j) cannot represent the lambda–mu ridge of the CR likelihood,
## so (i) the optimum need not sit on the ridge and (ii) the reported value is
## not an IS estimate of anything at that point.
##
## Method
##   Part A (noise-free): put the EXACT CR log-likelihood (DDD::bd_loglik, cond=0,
##     crown, unconditioned) on the same 150-point LHS grid .run_gam uses, fit the
##     package's own additive GAM (train_likelihood_GAM, spline_type="univariate",
##     the default) and run the package's own find_MLE. Compare the additive
##     optimum and its predicted value with the exact MLE (DDD::bd_ML, cond=0) and
##     the exact value at the additive optimum. Repeat with spline_type="bivariate".
##     Any gap here is the additive structure, not IS noise.
##   Part B (as shipped): estimate_rates(method="gam", model="cr", link="linear")
##     with n_grid=100, sample_size=100 (pipeline default 150/200; reduced to fit
##     the compute budget), two replicates per
##     tree. Then evaluate fhat with 2000 thinning trees AND the exact BDI fhat at
##     the GAM optimum and compare to the reported smooth value; compare pars with
##     DDD::bd_ML. The constant labelled-history offset between emphasis f and
##     DDD's branching-times likelihood is measured with BDI (whose lw is exact
##     for CR) at three points and verified constant.
##   Trees: bird.orders (README example, 23 tips) and a seeded 30-tip CR tree
##     with an interior MLE.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({library(emphasis); library(ape); library(DDD)})
options(width = 120)
set.seed(13)

CR  <- c(0L, 0L, 0L)
p8  <- function(l, m) c(l, 0, 0, 0, m, 0, 0, 0)
exact_ll <- function(l, m, brts)
  DDD::bd_loglik(pars1 = c(l, m, 0, 0), pars2 = c(0, 0, 0, 0, 2),
                 brts = brts, missnumspec = 0)
bdi_fhat <- function(l, m, brts, n = 20L)   # exact for CR (lw constant)
  emphasis:::.augment_tree_bdi(brts, p8(l, m), model_bin = CR,
                               sample_size = n, link = 0L, rho = 1)$fhat
if (!exists(".ess", envir = asNamespace("emphasis"))) {
  ess_fun <- function(lw) { w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }
} else ess_fun <- emphasis:::.ess
thin_fhat <- function(l, m, brts, n = 2000L) {
  a <- emphasis:::augment_trees(brts = brts, pars = p8(l, m), sample_size = n,
                                maxN = 200L * n, max_missing = 10000L,
                                max_lambda = 500, num_threads = 1L,
                                model = CR, link = 0L, rho = 1)
  lf <- emphasis:::eval_logf(p8(l, m), a$trees, model = CR, link = 0L, rho = 1)
  lw <- lf$logf - lf$logg
  c(fhat = emphasis:::.is_fhat(lf$logf, lf$logg,
                               n_zero_weight = emphasis:::.n0(a$rejected_zero_weights)),
    ess = ess_fun(lw), n = length(lw))
}

## ---- trees ------------------------------------------------------------------
data(bird.orders)
brts_bird <- sort(as.numeric(branching.times(bird.orders)), decreasing = TRUE)
sim <- TreeSim::sim.bd.taxa(n = 30, numbsim = 1, lambda = 0.5, mu = 0.25,
                            complete = FALSE)[[1]]
brts_sim <- sort(as.numeric(branching.times(sim)), decreasing = TRUE)
trees <- list(bird = brts_bird, sim30 = brts_sim)
## boxes sized so that thinning augmentation stays < ~1 s per grid point
## (cost grows with lambda * crown age: bird lambda=0.3,mu=0.06 already 3.4 s/pt)
boxes <- list(bird  = list(lb = c(0.01, 0.00), ub = c(0.15, 0.06)),
              sim30 = list(lb = c(0.10, 0.00), ub = c(0.80, 0.50)))

for (nm in names(trees)) {
  brts <- trees[[nm]]; lb <- boxes[[nm]]$lb; ub <- boxes[[nm]]$ub
  cat("\n==========", nm, ": n_tips =", length(brts) + 1,
      " crown =", round(brts[1], 3), " box lambda", lb[1], "-", ub[1],
      " mu", lb[2], "-", ub[2], "==========\n")

  ## exact MLE (unconditioned, crown, branching times)
  invisible(capture.output(
    ml <- DDD::bd_ML(brts = brts, initparsopt = c(0.2, 0.05), idparsopt = 1:2,
                     cond = 0, btorph = 0, soc = 2, verbose = FALSE)))
  l_ml <- ml$lambda0; m_ml <- ml$mu0; ll_ml <- ml$loglik
  cat(sprintf("DDD::bd_ML (cond=0): lambda=%.5f mu=%.5f loglik=%.4f\n", l_ml, m_ml, ll_ml))

  ## offset emphasis-f vs DDD, from BDI (exact for CR), at 3 points
  pts <- rbind(c(l_ml, max(m_ml, 0.01)), c(0.8 * l_ml + 0.1, 0.5 * l_ml),
               c(1.5 * l_ml, 0.2 * l_ml))
  offs <- apply(pts, 1, function(p) bdi_fhat(p[1], p[2], brts) - exact_ll(p[1], p[2], brts))
  cat(sprintf("offset (BDI fhat - bd_loglik) at 3 points: %s  -> spread %.2e\n",
              paste(round(offs, 6), collapse = " "), diff(range(offs))))
  off <- mean(offs)

  ## ---------------- Part A: additive GAM on the EXACT surface ----------------
  grid <- emphasis:::.lhs_sample(150L, lb, ub)
  colnames(grid) <- c("beta_0", "gamma_0")
  surf <- data.frame(grid, fhat = mapply(exact_ll, grid[, 1], grid[, 2],
                                         MoreArgs = list(brts = brts)),
                     n_trees = 1L)
  surf <- surf[is.finite(surf$fhat), ]
  for (st in c("univariate", "bivariate")) {
    g <- suppressMessages(emphasis:::train_likelihood_GAM(surf, par_names = c("beta_0", "gamma_0"),
                                               spline_type = st))
    res <- g$model$fhat - fitted(g)
    mle <- emphasis:::find_MLE(g, lower_bound = lb, upper_bound = ub,
                    par_names = c("beta_0", "gamma_0"))
    ll_at_opt <- exact_ll(mle$pars[1], mle$pars[2], brts)
    cat(sprintf(paste0("A[%s exact-surface] resid RMSE=%.3f max|resid|=%.3f | ",
                       "opt lambda=%.5f mu=%.5f (exact MLE %.5f %.5f) | ",
                       "reported=%.3f exact@opt=%.3f exact max=%.3f -> ",
                       "reported-exact@opt=%+.3f  reported-exactmax=%+.3f  ",
                       "exact@opt-exactmax=%+.3f\n"),
                st, sqrt(mean(res^2)), max(abs(res)), mle$pars[1], mle$pars[2],
                l_ml, m_ml, mle$loglik, ll_at_opt, ll_ml,
                mle$loglik - ll_at_opt, mle$loglik - ll_ml, ll_at_opt - ll_ml))
  }

  ## ---------------- Part B: as shipped --------------------------------------
  for (rep in 1:2) {
    t0 <- proc.time()[3]
    fit <- estimate_rates(brts, method = "gam", model = "cr", link = "linear",
                          control = list(lower_bound = lb, upper_bound = ub,
                                         n_grid = 100L, sample_size = 100L,   # pipeline default 150/200; reduced for budget
                                         num_threads = 1L, verbose = FALSE))
    el <- proc.time()[3] - t0
    lo <- fit$pars[["beta_0"]]; mo <- fit$pars[["gamma_0"]]
    th <- thin_fhat(lo, mo, brts, 2000L)
    bd <- bdi_fhat(lo, mo, brts)
    ex <- exact_ll(lo, mo, brts)
    g  <- fit$details$gam_fit
    sres <- g$model$fhat - fitted(g)
    cat(sprintf(paste0("B[rep %d, %.0fs] GAM pars lambda=%.5f mu=%.5f | reported loglik=%.3f ",
                       "AIC=%.3f | fhat2000@opt=%.3f (ESS %.0f) | BDI@opt=%.3f | ",
                       "exact@opt+off=%.3f | exact MLE+off=%.3f ",
                       "| reported-fhat2000=%+.3f reported-(exactmax+off)=%+.3f ",
                       "| IS-surface resid RMSE=%.3f | grid fhat range [%.1f, %.1f]\n"),
                rep, el, lo, mo, fit$loglik, fit$AIC, th[["fhat"]], th[["ess"]], bd,
                ex + off, ll_ml + off,
                fit$loglik - th[["fhat"]], fit$loglik - (ll_ml + off),
                sqrt(mean(sres^2)), min(fit$details$surface$fhat, na.rm = TRUE),
                max(fit$details$surface$fhat, na.rm = TRUE)))
    ## sanity: the reported loglik is literally the smooth's prediction
    pr <- as.numeric(predict(g, newdata = data.frame(beta_0 = lo, gamma_0 = mo)))
    cat(sprintf("     reported - predict(gam)@opt = %.2e (0 => it is the smooth value)\n",
                fit$loglik - pr))
  }
}

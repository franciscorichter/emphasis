## H13 replication (independent). Vary what the verifier did not:
##  - different trees (seeded 40-tip CR lambda=0.4 mu=0.2; 20-tip lambda=0.6 mu=0.1)
##  - the ACTUAL auto_bounds box for cr/linear (not a hand-chosen box)
##  - five different LHS grids on the EXACT surface (is the additive misplacement
##    systematic or grid noise?)
##  - the exponential link (surface in log lambda, log mu coordinates)
##  - as-shipped run with the pipeline default sample_size = 200 per grid point
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({library(emphasis); library(ape); library(DDD)})
options(width = 140)
t_start <- proc.time()[3]

CR <- c(0L, 0L, 0L)
p8 <- function(l, m) c(l, 0, 0, 0, m, 0, 0, 0)
exact_ll <- function(l, m, brts)
  DDD::bd_loglik(pars1 = c(l, m, 0, 0), pars2 = c(0, 0, 0, 0, 2),
                 brts = brts, missnumspec = 0)
ess_fun <- function(lw) { w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }
thin_fhat <- function(l, m, brts, n = 1000L) {
  a <- emphasis:::augment_trees(brts = brts, pars = p8(l, m), sample_size = n,
                                maxN = 200L * n, max_missing = 10000L,
                                max_lambda = 500, num_threads = 1L,
                                model = CR, link = 0L, rho = 1)
  lf <- emphasis:::eval_logf(p8(l, m), a$trees, model = CR, link = 0L, rho = 1)
  lw <- lf$logf - lf$logg
  c(fhat = emphasis:::.is_fhat(lf$logf, lf$logg,
                               n_zero_weight = emphasis:::.n0(a$rejected_zero_weights)),
    ess = ess_fun(lw))
}
bdi_fhat <- function(l, m, brts, n = 20L)
  emphasis:::.augment_tree_bdi(brts, p8(l, m), model_bin = CR,
                               sample_size = n, link = 0L, rho = 1)$fhat

set.seed(7)
t40 <- TreeSim::sim.bd.taxa(n = 40, numbsim = 1, lambda = 0.4, mu = 0.2, complete = FALSE)[[1]]
set.seed(99)
t20 <- TreeSim::sim.bd.taxa(n = 20, numbsim = 1, lambda = 0.6, mu = 0.1, complete = FALSE)[[1]]
trees <- list(sim40 = sort(as.numeric(branching.times(t40)), decreasing = TRUE),
              sim20 = sort(as.numeric(branching.times(t20)), decreasing = TRUE))

fitA <- function(surf, lb, ub, st, pn = c("beta_0", "gamma_0")) {
  g <- suppressMessages(capture.output(gg <- emphasis:::train_likelihood_GAM(
    surf, par_names = pn, spline_type = st)))
  res <- gg$model$fhat - fitted(gg)
  mle <- emphasis:::find_MLE(gg, lower_bound = lb, upper_bound = ub, par_names = pn)
  list(g = gg, rmse = sqrt(mean(res^2)), maxres = max(abs(res)), mle = mle)
}

for (nm in names(trees)) {
  brts <- trees[[nm]]
  cat("\n==========", nm, ": n_tips =", length(brts) + 1, " crown =", round(brts[1], 3), "==========\n")
  invisible(capture.output(
    ml <- DDD::bd_ML(brts = brts, initparsopt = c(0.3, 0.1), idparsopt = 1:2,
                     cond = 0, btorph = 0, soc = 2, verbose = FALSE)))
  l_ml <- ml$lambda0; m_ml <- ml$mu0; ll_ml <- ml$loglik
  cat(sprintf("DDD::bd_ML (cond=0): lambda=%.5f mu=%.5f loglik=%.4f\n", l_ml, m_ml, ll_ml))

  ## actual auto_bounds box for cr / linear (no survival GAM, to save time)
  invisible(capture.output(ab <- auto_bounds(brts, model = "cr", link = "linear",
                                             train_surv_gam = FALSE, verbose = FALSE)))
  lb <- as.numeric(ab$lower_bound); ub <- as.numeric(ab$upper_bound)
  cat(sprintf("auto_bounds box: lambda [%.4f, %.4f]  mu [%.4f, %.4f]  (MLE inside: %s)\n",
              lb[1], ub[1], lb[2], ub[2],
              l_ml >= lb[1] && l_ml <= ub[1] && m_ml >= lb[2] && m_ml <= ub[2]))

  ## ---- Part A: exact surface, 5 different LHS grids, additive vs bivariate ----
  for (st in c("univariate", "bivariate")) {
    out <- t(sapply(1:5, function(s) {
      set.seed(100 + s)
      grid <- emphasis:::.lhs_sample(150L, lb, ub); colnames(grid) <- c("beta_0", "gamma_0")
      surf <- data.frame(grid, fhat = mapply(exact_ll, grid[, 1], grid[, 2],
                                             MoreArgs = list(brts = brts)), n_trees = 1L)
      surf <- surf[is.finite(surf$fhat), ]
      f <- fitA(surf, lb, ub, st)
      ex <- exact_ll(f$mle$pars[1], f$mle$pars[2], brts)
      c(rmse = f$rmse, lambda = f$mle$pars[[1]], mu = f$mle$pars[[2]],
        reported = f$mle$loglik, exact_at_opt = ex,
        rep_minus_exact_at_opt = f$mle$loglik - ex,
        rep_minus_exactmax = f$mle$loglik - ll_ml,
        exact_at_opt_minus_max = ex - ll_ml)
    }))
    cat(sprintf("A[%s, linear, auto_bounds box] over 5 LHS grids:\n", st))
    print(round(out, 4))
  }

  ## ---- Part A': exponential link -- surface in (log lambda, log mu) ----
  lbe <- c(log(lb[1] + 1e-3), log(max(lb[2], 1e-3))); ube <- c(log(ub[1]), log(ub[2]))
  for (st in c("univariate", "bivariate")) {
    set.seed(11)
    grid <- emphasis:::.lhs_sample(150L, lbe, ube); colnames(grid) <- c("beta_0", "gamma_0")
    surf <- data.frame(grid, fhat = mapply(function(a, b) exact_ll(exp(a), exp(b), brts),
                                           grid[, 1], grid[, 2]), n_trees = 1L)
    surf <- surf[is.finite(surf$fhat), ]
    f <- fitA(surf, lbe, ube, st)
    ex <- exact_ll(exp(f$mle$pars[1]), exp(f$mle$pars[2]), brts)
    cat(sprintf(paste0("A'[%s, EXPONENTIAL link, log box lambda[%.3f,%.3f] mu[%.4f,%.3f]] rmse=%.3f | ",
                       "opt lambda=%.4f mu=%.4f | reported=%.3f exact@opt=%.3f exactmax=%.3f | ",
                       "rep-exact@opt=%+.3f rep-max=%+.3f exact@opt-max=%+.3f\n"),
                st, exp(lbe[1]), exp(ube[1]), exp(lbe[2]), exp(ube[2]), f$rmse,
                exp(f$mle$pars[1]), exp(f$mle$pars[2]), f$mle$loglik, ex, ll_ml,
                f$mle$loglik - ex, f$mle$loglik - ll_ml, ex - ll_ml))
  }

  ## ---- Part B: as shipped, pipeline-default 200 trees per grid point ----
  ## offset emphasis f vs DDD loglik (BDI exact for CR), checked at two points
  offs <- c(bdi_fhat(l_ml, max(m_ml, 0.01), brts) - exact_ll(l_ml, max(m_ml, 0.01), brts),
            bdi_fhat(1.3 * l_ml, 0.5 * l_ml, brts) - exact_ll(1.3 * l_ml, 0.5 * l_ml, brts))
  cat(sprintf("offset BDI-DDD at 2 pts: %s (spread %.1e)\n", paste(round(offs, 5), collapse = " "), diff(offs)))
  off <- mean(offs)
  if (nm == "sim40") {
    for (rep in 1:2) {
      t0 <- proc.time()[3]
      fit <- estimate_rates(brts, method = "gam", model = "cr", link = "linear",
                            control = list(lower_bound = lb, upper_bound = ub,
                                           n_grid = 60L, sample_size = 200L,
                                           num_threads = 1L, verbose = FALSE))
      el <- proc.time()[3] - t0
      lo <- fit$pars[["beta_0"]]; mo <- fit$pars[["gamma_0"]]
      th <- thin_fhat(lo, mo, brts, 1000L)
      ex <- exact_ll(lo, mo, brts)
      g <- fit$details$gam_fit; sres <- g$model$fhat - fitted(g)
      pr <- as.numeric(predict(g, newdata = data.frame(beta_0 = lo, gamma_0 = mo)))
      cat(sprintf(paste0("B[rep %d, %.0fs, 60x200] pars lambda=%.4f mu=%.4f | reported=%.3f AIC=%.3f | ",
                         "fhat1000@opt=%.3f (ESS %.0f) | exact@opt+off=%.3f | exactmax+off=%.3f | ",
                         "rep-fhat1000=%+.3f rep-(exactmax+off)=%+.3f | IS resid RMSE=%.3f | rep-predict=%.1e\n"),
                  rep, el, lo, mo, fit$loglik, fit$AIC, th[["fhat"]], th[["ess"]],
                  ex + off, ll_ml + off, fit$loglik - th[["fhat"]], fit$loglik - (ll_ml + off),
                  sqrt(mean(sres^2)), fit$loglik - pr))
    }
    ## same box, bivariate spline_type (the shipped alternative)
    t0 <- proc.time()[3]
    fitb <- estimate_rates(brts, method = "gam", model = "cr", link = "linear",
                           control = list(lower_bound = lb, upper_bound = ub,
                                          n_grid = 60L, sample_size = 200L, spline_type = "bivariate",
                                          num_threads = 1L, verbose = FALSE))
    lo <- fitb$pars[["beta_0"]]; mo <- fitb$pars[["gamma_0"]]
    th <- thin_fhat(lo, mo, brts, 1000L); ex <- exact_ll(lo, mo, brts)
    cat(sprintf(paste0("B[bivariate, %.0fs] pars lambda=%.4f mu=%.4f | reported=%.3f | fhat1000@opt=%.3f | ",
                       "exact@opt+off=%.3f | exactmax+off=%.3f | rep-fhat1000=%+.3f rep-(exactmax+off)=%+.3f\n"),
                proc.time()[3] - t0, lo, mo, fitb$loglik, th[["fhat"]], ex + off, ll_ml + off,
                fitb$loglik - th[["fhat"]], fitb$loglik - (ll_ml + off)))
  }
}
cat(sprintf("\nTotal elapsed %.0f s\n", proc.time()[3] - t_start))

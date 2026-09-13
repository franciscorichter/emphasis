## H73 — auto_bounds / .find_feasible_center / .wide_bounds put the gaussian-link
## intercept on a log scale although the C++ gaussian link uses beta_0 on the
## natural scale (rate = beta_0 * exp(-(eta_cov-1)^2/2)).
## Self-contained; ~1-2 min.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
data(bird.orders)
set.seed(1)

max_t  <- max(branching.times(bird.orders)); n_tips <- Ntip(bird.orders)
r_hat  <- log(n_tips / 2) / max_t; mu_hat <- 0.2 * r_hat; lam_hat <- r_hat + mu_hat
cat(sprintf("bird.orders: T=%.2f N=%d r_hat=%.4f lam_hat=%.4f log(lam_hat)=%.4f\n",
            max_t, n_tips, r_hat, lam_hat, log(lam_hat)))

## 1. Direct check of the link semantics: gaussian link, beta_0 on natural scale
##    (cr, so eta_cov = 0 and rate = beta_0 * exp(-1/2)).
sim_at <- function(b0, g0, link, model = "cr", n = 10) {
  pm <- if (model == "cr") cbind(b0, g0) else cbind(b0, 0, g0, 0)
  pm <- pm[rep(1, n), , drop = FALSE]
  s <- simulate_tree(pars = pm, max_t = max_t, model = model, link = link,
                     max_tries = 0, max_lin = 500L, num_threads = 1L)
  st <- sapply(s$simulations, `[[`, "status")
  nt <- sapply(s$simulations, function(x) if (x$status == "done" && !is.null(x$tes)) Ntip(x$tes) else 0L)
  list(status = st, ntips = nt)
}
cat("\n[1] gaussian link, cr, beta_0 = log(lam_hat) (the auto_bounds centre):\n")
a <- sim_at(log(lam_hat), log(mu_hat), "gaussian"); print(table(a$status)); print(a$ntips)
cat("[1b] gaussian link, cr, beta_0 = lam_hat (natural scale):\n")
b <- sim_at(lam_hat, mu_hat, "gaussian"); print(table(b$status)); print(b$ntips)
cat("[1c] exponential link, cr, beta_0 = log(lam_hat) (what the branch was written for):\n")
d <- sim_at(log(lam_hat), log(mu_hat), "exponential"); print(table(d$status)); print(d$ntips)

## 2. .test_feasibility at the negative centre (tip_lo = max(2, floor(0.1*23)) = 2)
tip_lo <- max(2, floor(n_tips * 0.1)); tip_hi <- ceiling(n_tips * 10)
centre_dd <- c(log(lam_hat), 0, log(mu_hat), 0)
feas <- emphasis:::.test_feasibility(centre_dd, "dd", "gaussian", max_t, 500L,
                                     n_test = 10L, tip_lo = tip_lo, tip_hi = tip_hi)
cat(sprintf("\n[2] tip_lo=%d; .test_feasibility(centre with beta_0=%.3f, gaussian) = %s\n",
            tip_lo, centre_dd[1], feas))

## 3. auto_bounds itself (small settings)
cat("\n[3] auto_bounds(bird.orders, 'dd', 'gaussian')\n")
ab <- auto_bounds(bird.orders, model = "dd", link = "gaussian",
                  n_test = 3L, bisect_steps = 4L, train_surv_gam = TRUE,
                  verbose = TRUE, num_threads = 1L)
cat("centre:\n"); print(ab$center)
cat("lower:\n");  print(ab$lower_bound)
cat("upper:\n");  print(ab$upper_bound)
cat(sprintf("ASSERT lower_bound['beta_0'] > 0 : %s\n", ab$lower_bound["beta_0"] > 0))
cat(sprintf("ASSERT upper_bound['beta_0'] > 0 : %s\n", ab$upper_bound["beta_0"] > 0))

## 3b. Where did the survival GAM get trained? Reproduce the LHS on the returned box.
gm <- emphasis:::.lhs_sample(500L, ab$lower_bound, ab$upper_bound)
colnames(gm) <- names(ab$lower_bound)
cat(sprintf("fraction of survival-GAM training points with beta_0 <= 0: %.3f\n",
            mean(gm[, "beta_0"] <= 0)))
if (!is.null(ab$survival_gam)) {
  nd_neg <- data.frame(beta_0 = -1, beta_N = 0, gamma_0 = -1, gamma_N = 0)
  nd_pos <- data.frame(beta_0 = lam_hat, beta_N = 0, gamma_0 = mu_hat, gamma_N = 0)
  cat(sprintf("predict_survival at beta_0=-1 (negative rate): %.3f\n",
              predict_survival(ab$survival_gam, nd_neg)))
  cat(sprintf("predict_survival at beta_0=lam_hat=%.3f (true-ish rate): %.3f\n",
              lam_hat, predict_survival(ab$survival_gam, nd_pos)))
}

## 4. .wide_bounds for gaussian vs linear; the pipeline's init_pars midpoint
wb_g <- emphasis:::.wide_bounds(c(1L,0L,0L), 2L, max_t, n_tips)
wb_l <- emphasis:::.wide_bounds(c(1L,0L,0L), 0L, max_t, n_tips)
cat("\n[4] .wide_bounds gaussian lb/ub:\n"); print(rbind(lb = wb_g$lb, ub = wb_g$ub))
cat(".wide_bounds linear lb/ub:\n");         print(rbind(lb = wb_l$lb, ub = wb_l$ub))
cat(sprintf("gaussian midpoint init beta_0 = %.4f (rate at eta_cov=0 would be %.4f)\n",
            mean(c(wb_g$lb[1], wb_g$ub[1])), mean(c(wb_g$lb[1], wb_g$ub[1])) * exp(-0.5)))

## 5. Control: linear link is on the link_int == 0 branch and unaffected
cat("\n[5] auto_bounds(bird.orders, 'cr', 'linear') — control\n")
abl <- auto_bounds(bird.orders, model = "cr", link = "linear",
                   n_test = 3L, bisect_steps = 4L, train_surv_gam = FALSE,
                   verbose = FALSE, num_threads = 1L)
print(rbind(lower = abl$lower_bound, upper = abl$upper_bound))
cat(sprintf("linear: lower_bound['beta_0'] > 0 : %s\n", abl$lower_bound["beta_0"] > 0))

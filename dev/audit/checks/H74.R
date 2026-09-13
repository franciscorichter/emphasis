## H74 — does auto_bounds() exclude high-turnover (mu ~ lambda) parameters
## that survival conditioning makes admissible?
##
## Design:
##   1. Simulate CR trees (linear link) with lambda = 1, mu = 0.9 (mu/lambda = 0.9),
##      crown age T = 20, conditioned on crown survival via retries, keeping trees
##      with 15..60 extant tips.
##   2. For each tree run auto_bounds(model = "cr", link = "linear",
##      train_surv_gam = FALSE) and record the gamma_0 (mu) box and beta_0 box.
##   3. Check whether (a) the generating point (1, 0.9) and (b) DDD::bd_ML(cond = 1)
##      lie inside the returned box.
##   4. Independently measure what .test_feasibility sees at the true point:
##      the fraction of unconditional forward sims that are "done" with
##      tips in [max(2, 0.1N), 10N], vs the analytic crown-survival probability.
##   5. Verify the box is a hard constraint: pass the box to estimate_rates(cr,
##      mcem) with init at the clamped edge and confirm the estimate cannot exceed ub.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(DDD) })

lam <- 1; mu <- 0.9; T0 <- 20
n_trees <- 6L

# ---- 1. simulate trees conditioned on survival, moderate size --------------
trees <- list(); tries <- 0L
while (length(trees) < n_trees && tries < 400L) {
  tries <- tries + 1L
  s <- simulate_tree(pars = c(lam, mu), max_t = T0, model = "cr", link = "linear",
                     max_tries = 200, max_lin = 5000)
  if (s$status == "done" && !is.null(s$tes)) {
    N <- Ntip(s$tes)
    if (N >= 15 && N <= 60) trees[[length(trees) + 1L]] <- s$tes
  }
}
cat(sprintf("simulated %d trees in %d attempts; tips: %s\n",
            length(trees), tries, paste(sapply(trees, Ntip), collapse = " ")))

# ---- analytic survival prob of the crown (both lineages) to T under CR ------
p_surv_lineage <- function(lam, mu, T) { r <- lam - mu; r / (lam - mu * exp(-r * T)) }
p_crown <- p_surv_lineage(lam, mu, T0)^2
cat(sprintf("analytic P(one lineage survives T=%g) = %.4f ; P(crown survives) = %.4f\n",
            T0, p_surv_lineage(lam, mu, T0), p_crown))

# ---- 2-4. per tree --------------------------------------------------------
res <- data.frame()
for (k in seq_along(trees)) {
  tr <- trees[[k]]; N <- Ntip(tr); Tk <- max(branching.times(tr))
  brts <- sort(branching.times(tr), decreasing = TRUE)

  ab <- auto_bounds(tr, model = "cr", link = "linear", train_surv_gam = FALSE,
                    verbose = FALSE)
  lb <- ab$lower_bound; ub <- ab$upper_bound

  # DDD MLE conditioned on survival (cond = 1), crown (soc = 2)
  ml <- tryCatch(bd_ML(brts = brts, initparsopt = c(0.5, 0.3), idparsopt = 1:2,
                       cond = 1, soc = 2, btorph = 1, verbose = FALSE),
                 error = function(e) NULL)
  ml_lam <- if (is.null(ml)) NA else ml$lambda0
  ml_mu  <- if (is.null(ml)) NA else ml$mu0

  # what the feasibility criterion sees at the true point (n = 200 sims)
  tip_lo <- max(2, floor(N * 0.1)); tip_hi <- ceiling(N * 10)
  max_lin <- as.integer(max(20 * N, 500))
  pm <- matrix(rep(c(lam, mu), 200), nrow = 200, byrow = TRUE)
  sims <- simulate_tree(pars = pm, max_t = Tk, model = "cr", link = "linear",
                        max_tries = 0, max_lin = max_lin)
  nt <- sapply(sims$simulations, function(s)
    if (s$status == "done" && !is.null(s$tes)) length(s$tes$tip.label) else 0L)
  feas_frac_true <- mean(nt >= tip_lo & nt <= tip_hi)
  done_frac_true <- mean(nt > 0)
  crit_true <- emphasis:::.test_feasibility(c(lam, mu), "cr", "linear", Tk, max_lin,
                                            n_test = 5L, tip_lo = tip_lo, tip_hi = tip_hi)

  # same at the DDD MLE
  feas_frac_ml <- NA
  if (!is.na(ml_mu)) {
    pm2 <- matrix(rep(c(ml_lam, ml_mu), 200), nrow = 200, byrow = TRUE)
    sims2 <- simulate_tree(pars = pm2, max_t = Tk, model = "cr", link = "linear",
                           max_tries = 0, max_lin = max_lin)
    nt2 <- sapply(sims2$simulations, function(s)
      if (s$status == "done" && !is.null(s$tes)) length(s$tes$tip.label) else 0L)
    feas_frac_ml <- mean(nt2 >= tip_lo & nt2 <= tip_hi)
  }

  inside_true <- all(c(lam, mu) >= lb) && all(c(lam, mu) <= ub)
  inside_ml   <- if (is.na(ml_mu)) NA else
    all(c(ml_lam, ml_mu) >= lb) && all(c(ml_lam, ml_mu) <= ub)

  res <- rbind(res, data.frame(
    tree = k, N = N, T = round(Tk, 2),
    lb_lam = lb[1], ub_lam = ub[1], lb_mu = lb[2], ub_mu = ub[2],
    max_turnover_in_box = ub[2] / lb[1],
    ml_lam = ml_lam, ml_mu = ml_mu,
    inside_true = inside_true, inside_ml = inside_ml,
    feas_frac_true = feas_frac_true, done_frac_true = done_frac_true,
    test_feas_true = crit_true, feas_frac_ml = feas_frac_ml,
    center_lam = ab$center[1], center_mu = ab$center[2]))
  cat(sprintf("tree %d (N=%d, T=%.1f): box lam=[%.3f,%.3f] mu=[%.3f,%.3f]; true (1,0.9) inside=%s; bd_ML=(%.3f,%.3f) inside=%s; feas.frac@true=%.3f (done=%.3f); .test_feasibility@true=%s\n",
              k, N, Tk, lb[1], ub[1], lb[2], ub[2], inside_true, ml_lam, ml_mu,
              inside_ml, feas_frac_true, done_frac_true, crit_true))
}
print(res, digits = 3)
cat(sprintf("\nSUMMARY: true point inside box in %d/%d trees; bd_ML(cond=1) inside box in %d/%d trees\n",
            sum(res$inside_true), nrow(res), sum(res$inside_ml, na.rm = TRUE), sum(!is.na(res$inside_ml))))
cat(sprintf("max mu/lambda ratio reachable inside box (ub_mu/lb_lam upper limit): %s\n",
            paste(round(res$max_turnover_in_box, 2), collapse = " ")))
cat(sprintf("max mu/lambda at the box corner (ub_mu / ub_lam): %s\n",
            paste(round(res$ub_mu / res$ub_lam, 2), collapse = " ")))

# ---- 5. the box is a hard constraint on the fit ---------------------------
tr <- trees[[1]]; ab1 <- auto_bounds(tr, model = "cr", link = "linear",
                                     train_surv_gam = FALSE, verbose = FALSE)
init <- c(ab1$upper_bound[1], ab1$upper_bound[2])  # start at the mu ceiling
fit <- tryCatch(estimate_rates(tr, model = "cr", link = "linear", method = "mcem",
                               init_pars = init, cond = NULL,
                               control = list(lower_bound = ab1$lower_bound,
                                              upper_bound = ab1$upper_bound,
                                              max_iter = 6, max_time = 60,
                                              sample_size = 50, num_threads = 1)),
                error = function(e) { cat("estimate_rates error:", conditionMessage(e), "\n"); NULL })
if (!is.null(fit)) {
  cat(sprintf("mcem fit from mu-ceiling init: pars=(%.3f, %.3f); ub=(%.3f, %.3f); mu <= ub_mu: %s\n",
              fit$pars[1], fit$pars[2], ab1$upper_bound[1], ab1$upper_bound[2],
              fit$pars[2] <= ab1$upper_bound[2] + 1e-9))
}

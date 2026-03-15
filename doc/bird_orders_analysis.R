#!/usr/bin/env Rscript
# ============================================================================
# Bird orders model comparison — emphasis package
#
# Fits CR, DD, PD, EP models to ape::bird.orders (23 tips) via GAM method
# (one-shot, non-iterative), compares via AIC.
#
# Usage:  Rscript bird_orders_analysis.R
# Output: console table + bird_analysis_results.RData
# ============================================================================

# ── Install emphasis from GitHub ─────────────────────────────────────────────
if (requireNamespace("emphasis", quietly = TRUE)) {
  remove.packages("emphasis")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("franciscorichter/emphasis", force = TRUE)

library(emphasis)
library(ape)

cat("emphasis version:", as.character(packageVersion("emphasis")), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

# ── Data ─────────────────────────────────────────────────────────────────────

data(bird.orders)
tree <- bird.orders
cat("Tree: bird.orders,", Ntip(tree), "tips\n\n")

# ── Model specifications ────────────────────────────────────────────────────

models <- list(
  CR = list(
    model = "cr", link = "linear",
    lb = c(0.01, 0), ub = c(2.0, 1.0)
  ),
  DD = list(
    model = "dd", link = "exponential",
    lb = c(-1, -0.3, -3, -0.1),
    ub = c( 2,  0.0,  0,  0.1)
  ),
  PD = list(
    model = "pd", link = "exponential",
    lb = c(-1, -0.5, -3, -0.5),
    ub = c( 2,  0.5,  0,  0.5)
  ),
  EP = list(
    model = "ep", link = "linear",
    lb = c(0.01, -0.2, 0, -0.2),
    ub = c(2.0,   0.5, 1,  0.5)
  )
)

# ── GAM control ──────────────────────────────────────────────────────────────
# CR: 2 params -> factorial grid (grid_points^2)
# DD/PD/EP: 4 params -> LHS grid (n_grid points)

gam_ctrl_2p <- list(grid_points = 10, sample_size = 50)   # 10^2 = 100 grid pts
gam_ctrl_4p <- list(n_grid = 100, sample_size = 50)        # 100 LHS points

# ── Fit all models ───────────────────────────────────────────────────────────

results <- list()
set.seed(2026)

for (mod_name in names(models)) {
  spec <- models[[mod_name]]
  n_pars <- length(spec$lb)
  ctrl <- if (n_pars == 2) gam_ctrl_2p else gam_ctrl_4p
  ctrl$lower_bound <- spec$lb
  ctrl$upper_bound <- spec$ub

  cat(sprintf("[%s] model=%s link=%s (%d pars) ... ", mod_name, spec$model, spec$link, n_pars))

  t0 <- proc.time()
  fit <- tryCatch(
    estimate_rates(tree, method = "gam", model = spec$model, link = spec$link,
                   control = ctrl),
    error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
  )
  elapsed <- (proc.time() - t0)[3]

  if (!is.null(fit)) {
    cat(sprintf("%.1fs\n", elapsed))
    cat(sprintf("  pars: %s\n", paste(names(fit$pars), "=", round(fit$pars, 4), collapse = ", ")))
    cat(sprintf("  loglik = %.2f,  AIC = %.2f\n\n", fit$loglik, fit$AIC))
    results[[mod_name]] <- list(
      model = mod_name, link = spec$link, n_pars = n_pars,
      pars = fit$pars, loglik = fit$loglik, AIC = fit$AIC, elapsed = elapsed
    )
  } else {
    cat(sprintf("FAILED (%.1fs)\n\n", elapsed))
    results[[mod_name]] <- list(
      model = mod_name, link = spec$link, n_pars = n_pars,
      pars = rep(NA, n_pars), loglik = NA, AIC = NA, elapsed = elapsed
    )
  }
}

# ── Summary table ────────────────────────────────────────────────────────────

cat("\n========================================\n")
cat("  bird.orders — Model comparison (GAM)\n")
cat("========================================\n\n")

aics <- sapply(results, function(r) r$AIC)
ord  <- order(aics)
best <- min(aics, na.rm = TRUE)

cat(sprintf("%-8s %5s %10s %10s %10s %8s\n", "Model", "nPar", "loglik", "AIC", "dAIC", "Time"))
cat(strrep("-", 60), "\n")
for (i in ord) {
  r <- results[[i]]
  cat(sprintf("%-8s %5d %10.2f %10.2f %10.2f %7.1fs\n",
              r$model, r$n_pars, r$loglik, r$AIC, r$AIC - best, r$elapsed))
}
cat("\nBest parameters:\n")
best_mod <- results[[ord[1]]]
for (j in seq_along(best_mod$pars)) {
  cat(sprintf("  %s = %.4f\n", names(best_mod$pars)[j], best_mod$pars[j]))
}

save(results, file = "bird_analysis_results.RData")
cat("\nSaved to bird_analysis_results.RData\n")

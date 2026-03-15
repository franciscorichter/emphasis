#!/usr/bin/env Rscript
# ============================================================================
# Bird orders model comparison — emphasis package
#
# Fits CR, DD, PD, EP models to ape::bird.orders (23 tips) and
# ape::bird.families (137 tips) via GAM method.
#
# DDD reference values (dd_ML, conditioned on survival):
#   bird.orders:   lambda_0=0.693, mu~0,    K=23
#   bird.families: lambda_0=0.795, mu=0.102, K=201
#
# Translation to emphasis linear-link DD:
#   beta_0=lambda_0, beta_N=-(lambda_0-mu)/K, gamma_0=mu, gamma_N=0
#   bird.orders:   (0.693, -0.030, ~0,    0)
#   bird.families: (0.795, -0.0035, 0.102, 0)
#
# Usage:  Rscript bird_orders_analysis.R
# Output: console table + bird_analysis_results.RData
# ============================================================================

# ── Install emphasis from GitHub ─────────────────────────────────
if (requireNamespace("emphasis", quietly = TRUE)) {
  remove.packages("emphasis")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github(
  "franciscorichter/emphasis",
  force = TRUE
)

library(emphasis)
library(ape)

cat("emphasis version:",
    as.character(packageVersion("emphasis")), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

# ── Data ─────────────────────────────────────────────────────────

data(bird.orders)
data(bird.families)

datasets <- list(
  orders   = bird.orders,
  families = bird.families
)

# ── Model specifications ─────────────────────────────────────────
#
# DD_lin: linear link, directly comparable to DDD parametrization
#   lambda(N) = max(0, beta_0 + beta_N * N)
#   mu(N)     = max(0, gamma_0 + gamma_N * N)
#
# DD_exp: exponential link (recommended for IS stability)
#   lambda(N) = exp(beta_0 + beta_N * N)
#
# Bounds are set to cover DDD reference values with margin.

models <- list(
  CR = list(
    model = "cr", link = "linear",
    lb = c(0.01, 0),
    ub = c(2.0, 1.0)
  ),
  DD_lin = list(
    model = "dd", link = "linear",
    lb = c(0.1, -0.15, 0, -0.01),
    ub = c(2.0, 0.0, 0.5, 0.01)
  ),
  DD_exp = list(
    model = "dd", link = "exponential",
    lb = c(-2, -0.5, -8, -0.1),
    ub = c(2, 0.0, 0, 0.1)
  ),
  PD = list(
    model = "pd", link = "exponential",
    lb = c(-2, -0.5, -8, -0.5),
    ub = c(2, 0.5, 0, 0.5)
  ),
  EP = list(
    model = "ep", link = "linear",
    lb = c(0.01, -0.2, 0, -0.2),
    ub = c(2.0, 0.5, 1, 0.5)
  )
)

# ── GAM control ──────────────────────────────────────────────────
# CR: 2 params -> factorial grid (grid_points^2)
# Others: 4 params -> LHS grid (n_grid points)

gam_ctrl_2p <- list(
  grid_points = 10,
  sample_size = 50
)
gam_ctrl_4p <- list(
  n_grid = 100,
  sample_size = 50
)

# ── Fit all models on all datasets ───────────────────────────────

all_results <- list()
set.seed(2026)

for (ds_name in names(datasets)) {
  tree <- datasets[[ds_name]]
  n_tips <- Ntip(tree)
  cat(sprintf(
    "=== %s (%d tips) ===\n\n",
    ds_name, n_tips
  ))

  for (mod_name in names(models)) {
    spec <- models[[mod_name]]
    n_pars <- length(spec$lb)
    ctrl <- if (n_pars == 2) gam_ctrl_2p else gam_ctrl_4p
    ctrl$lower_bound <- spec$lb
    ctrl$upper_bound <- spec$ub

    cat(sprintf(
      "  [%s] %s link=%s (%dp) ... ",
      mod_name, spec$model, spec$link, n_pars
    ))

    t0 <- proc.time()
    fit <- tryCatch(
      estimate_rates(
        tree,
        method  = "gam",
        model   = spec$model,
        link    = spec$link,
        control = ctrl
      ),
      error = function(e) {
        cat("ERROR:", e$message, "\n")
        NULL
      }
    )
    elapsed <- (proc.time() - t0)[3]

    if (!is.null(fit)) {
      cat(sprintf("%.1fs\n", elapsed))
      pars_str <- paste(
        names(fit$pars), "=",
        round(fit$pars, 4),
        collapse = ", "
      )
      cat(sprintf("    pars: %s\n", pars_str))
      cat(sprintf(
        "    loglik=%.2f  AIC=%.2f\n\n",
        fit$loglik, fit$AIC
      ))
      key <- paste(ds_name, mod_name, sep = "_")
      all_results[[key]] <- list(
        dataset = ds_name,
        model   = mod_name,
        link    = spec$link,
        n_tips  = n_tips,
        n_pars  = n_pars,
        pars    = fit$pars,
        loglik  = fit$loglik,
        AIC     = fit$AIC,
        elapsed = elapsed
      )
    } else {
      cat(sprintf("FAILED (%.1fs)\n\n", elapsed))
      key <- paste(ds_name, mod_name, sep = "_")
      all_results[[key]] <- list(
        dataset = ds_name,
        model   = mod_name,
        link    = spec$link,
        n_tips  = n_tips,
        n_pars  = n_pars,
        pars    = rep(NA, n_pars),
        loglik  = NA,
        AIC     = NA,
        elapsed = elapsed
      )
    }
  }
}

# ── Summary tables ───────────────────────────────────────────────

cat("\n")
cat("====================================\n")
cat("         RESULTS SUMMARY\n")
cat("====================================\n")

for (ds_name in names(datasets)) {
  n_tips <- Ntip(datasets[[ds_name]])
  cat(sprintf(
    "\n--- %s (%d tips) ---\n",
    ds_name, n_tips
  ))

  ds_res <- Filter(
    function(r) r$dataset == ds_name,
    all_results
  )
  aics <- sapply(ds_res, function(r) r$AIC)
  ord <- order(aics)
  best_aic <- min(aics, na.rm = TRUE)

  hdr <- sprintf(
    "%-10s %-6s %4s %10s %10s %10s %7s",
    "Model", "Link", "nP", "loglik",
    "AIC", "dAIC", "Time"
  )
  cat(hdr, "\n")
  cat(strrep("-", nchar(hdr)), "\n")

  for (i in ord) {
    r <- ds_res[[i]]
    daic <- r$AIC - best_aic
    cat(sprintf(
      "%-10s %-6s %4d %10.2f %10.2f %10.2f %6.1fs\n",
      r$model, r$link, r$n_pars,
      r$loglik, r$AIC, daic, r$elapsed
    ))
  }
}

# ── Save ─────────────────────────────────────────────────────────

save(all_results, file = "bird_analysis_results.RData")
cat("\nSaved to bird_analysis_results.RData\n")

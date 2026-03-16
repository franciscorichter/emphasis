#!/usr/bin/env Rscript
# ============================================================================
# Bird orders model comparison — two-stage pipeline (GAM → MCEM)
#
# Tree: bird.orders (23 tips, ape package)
# Models: DD_lin, DD_exp, PD_exp, EP_lin, EP_exp
# (PD_lin skipped — linear link + large P → IS collapse)
#
# Pipeline per model:
#   Stage 1: auto_bounds() → detect bounds + train survival GAM
#   Stage 2: estimate_rates(method="gam", cond=surv_gam) → initial params
#   Stage 3: estimate_rates(method="mcem", init_pars=gam$pars, cond=surv_gam)
#
# All results conditioned on survival (comparable to DDD cond=1).
# DDD reference: dd_ML(bird.orders) → lam0≈0.693, mu≈0, K≈23
#
# Usage:  Rscript doc/bird_orders_analysis.R
# Output: console + doc/bird_orders_results.RData
# ============================================================================

library(emphasis)
library(ape)

cat("emphasis version:", as.character(packageVersion("emphasis")), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

# ── Data ─────────────────────────────────────────────────────────

data(bird.orders)
tree <- bird.orders
n_tips <- Ntip(tree)
cat(sprintf("Tree: bird.orders (%d tips)\n\n", n_tips))

# ── Model configurations ────────────────────────────────────────
# Skip PD_lin (IS collapse with linear link + large PD values)

configs <- list(
  DD_lin = list(model = "dd", link = "linear"),
  DD_exp = list(model = "dd", link = "exponential"),
  PD_exp = list(model = "pd", link = "exponential"),
  EP_lin = list(model = "ep", link = "linear"),
  EP_exp = list(model = "ep", link = "exponential")
)

# ── Tuning ───────────────────────────────────────────────────────

n_cores <- max(1L, parallel::detectCores() - 1L)
gam_ctrl <- list(n_grid = 150, sample_size = 200)
mcem_ctrl <- list(
  sample_size = 200,
  maxN        = 5000,
  max_iter    = 100,
  tol         = 0.01,
  patience    = 3,
  verbose     = TRUE
)

# ── Two-stage pipeline ──────────────────────────────────────────

results <- list()
set.seed(2026)

for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  cat(sprintf("\n══════ %s (model=%s, link=%s) ══════\n",
              cfg_name, cfg$model, cfg$link))

  # --- Stage 1: auto_bounds + survival GAM ---
  cat("\n[Stage 1] auto_bounds...\n")
  t0 <- proc.time()
  ab <- tryCatch(
    auto_bounds(tree, model = cfg$model, link = cfg$link,
                num_threads = n_cores, verbose = TRUE),
    error = function(e) {
      cat("  FAILED:", e$message, "\n")
      NULL
    }
  )
  if (is.null(ab)) {
    results[[cfg_name]] <- list(
      config = cfg, n_tips = n_tips,
      gam_fit = NULL, mcem_fit = NULL,
      pars = NULL, loglik = NA, AIC = NA,
      elapsed = (proc.time() - t0)[3], status = "bounds_failed"
    )
    next
  }
  t_bounds <- (proc.time() - t0)[3]
  cat(sprintf("  bounds detected in %.0fs\n", t_bounds))
  cat(sprintf("  lb: %s\n", paste(round(ab$lower_bound, 4), collapse=", ")))
  cat(sprintf("  ub: %s\n", paste(round(ab$upper_bound, 4), collapse=", ")))

  # --- Stage 2: GAM fit (initial params) ---
  cat("\n[Stage 2] GAM fit...\n")
  ctrl_gam <- gam_ctrl
  ctrl_gam$lower_bound <- ab$lower_bound
  ctrl_gam$upper_bound <- ab$upper_bound
  ctrl_gam$num_threads <- n_cores

  t1 <- proc.time()
  gam_fit <- tryCatch(
    estimate_rates(tree, method = "gam",
                   model = cfg$model, link = cfg$link,
                   cond = ab$survival_gam,
                   control = ctrl_gam),
    error = function(e) {
      cat("  GAM FAILED:", e$message, "\n")
      NULL
    }
  )
  t_gam <- (proc.time() - t1)[3]

  if (!is.null(gam_fit)) {
    cat(sprintf("  GAM done in %.0fs\n", t_gam))
    cat(sprintf("  GAM pars: %s\n",
                paste(names(gam_fit$pars), "=", round(gam_fit$pars, 4),
                      collapse = ", ")))
    cat(sprintf("  GAM loglik=%.2f  AIC=%.2f\n", gam_fit$loglik, gam_fit$AIC))
  }

  # --- Stage 3: MCEM refinement ---
  cat("\n[Stage 3] MCEM refinement...\n")
  ctrl_mcem <- mcem_ctrl
  ctrl_mcem$lower_bound <- ab$lower_bound
  ctrl_mcem$upper_bound <- ab$upper_bound
  ctrl_mcem$num_threads <- n_cores

  init <- if (!is.null(gam_fit)) gam_fit$pars else NULL

  t2 <- proc.time()
  mcem_fit <- tryCatch(
    estimate_rates(tree, method = "mcem",
                   model = cfg$model, link = cfg$link,
                   init_pars = init,
                   cond = ab$survival_gam,
                   control = ctrl_mcem),
    error = function(e) {
      cat("  MCEM FAILED:", e$message, "\n")
      NULL
    }
  )
  t_mcem <- (proc.time() - t2)[3]
  elapsed <- t_bounds + t_gam + t_mcem

  # Use MCEM if available and has finite loglik, else fall back to GAM
  mcem_ok <- !is.null(mcem_fit) && is.finite(mcem_fit$loglik)
  gam_ok  <- !is.null(gam_fit) && is.finite(gam_fit$loglik)

  if (mcem_ok) {
    final_fit <- mcem_fit
    stage     <- "mcem"
  } else if (gam_ok) {
    final_fit <- gam_fit
    stage     <- "gam_fallback"
    if (!is.null(mcem_fit))
      cat("  MCEM returned NA loglik — falling back to GAM result\n")
  } else {
    final_fit <- NULL
    stage     <- "failed"
  }

  if (!is.null(mcem_fit)) {
    cat(sprintf("  MCEM done in %.0fs\n", t_mcem))
    if (mcem_ok) {
      cat(sprintf("  MCEM pars: %s\n",
                  paste(names(mcem_fit$pars), "=", round(mcem_fit$pars, 4),
                        collapse = ", ")))
      cat(sprintf("  MCEM loglik=%.2f  AIC=%.2f\n",
                  mcem_fit$loglik, mcem_fit$AIC))
    } else {
      cat("  MCEM loglik=NA (E-step collapse)\n")
    }
  }

  cat(sprintf("  Final result from: %s\n", stage))

  results[[cfg_name]] <- list(
    config   = cfg,
    n_tips   = n_tips,
    gam_fit  = gam_fit,
    mcem_fit = mcem_fit,
    pars     = if (!is.null(final_fit)) final_fit$pars else NULL,
    loglik   = if (!is.null(final_fit)) final_fit$loglik else NA,
    AIC      = if (!is.null(final_fit)) final_fit$AIC else NA,
    n_pars   = if (!is.null(final_fit)) final_fit$n_pars else NA,
    elapsed  = elapsed,
    status   = stage
  )
}

# ── DDD reference validation ────────────────────────────────────

cat("\n\n══════ DDD REFERENCE VALIDATION ══════\n")
ddd_ref <- NULL
if (requireNamespace("DDD", quietly = TRUE)) {
  cat("Running DDD::dd_ML on bird.orders...\n")
  ddd_fit <- tryCatch(
    DDD::dd_ML(brts = ape::branching.times(tree), cond = 1,
               btorph = 1, initparsopt = c(0.7, 0.01, 30)),
    error = function(e) { cat("DDD failed:", e$message, "\n"); NULL }
  )
  if (!is.null(ddd_fit)) {
    ddd_ref <- list(
      lambda = ddd_fit$lambda,
      mu     = ddd_fit$mu,
      K      = ddd_fit$K,
      loglik = ddd_fit$loglik,
      AIC    = -2 * ddd_fit$loglik + 2 * 3
    )
    cat(sprintf("  DDD: lambda=%.4f  mu=%.4f  K=%.1f  loglik=%.2f  AIC=%.2f\n",
                ddd_ref$lambda, ddd_ref$mu, ddd_ref$K,
                ddd_ref$loglik, ddd_ref$AIC))
  }
} else {
  cat("DDD package not installed — skipping reference validation.\n")
  cat("Known DDD reference: lambda=0.693, mu~0, K=23, loglik~-49\n")
  ddd_ref <- list(lambda = 0.693, mu = 0, K = 23, loglik = -49, AIC = -92)
}

# ── Summary table ────────────────────────────────────────────────

cat("\n\n============================================\n")
cat("      BIRD ORDERS — MODEL COMPARISON\n")
cat("============================================\n\n")

aics <- sapply(results, function(r) r$AIC)
ord  <- order(aics, na.last = TRUE)
best_aic <- min(aics, na.rm = TRUE)

hdr <- sprintf("%-10s %-6s %3s %8s %10s %10s %8s %8s",
               "Model", "Link", "nP", "Status",
               "loglik", "AIC", "dAIC", "Time")
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")

for (i in ord) {
  r <- results[[i]]
  nm <- names(results)[i]
  if (is.na(r$AIC)) {
    cat(sprintf("%-10s %-6s %3s %8s %10s %10s %8s %7.0fs\n",
                nm, r$config$link, "---", r$status,
                "---", "---", "---", r$elapsed))
  } else {
    cat(sprintf("%-10s %-6s %3d %8s %10.2f %10.2f %8.2f %7.0fs\n",
                nm, r$config$link, r$n_pars, r$status,
                r$loglik, r$AIC, r$AIC - best_aic, r$elapsed))
  }
}

# Best model details
valid <- is.finite(aics)
if (any(valid)) {
  best_key <- names(results)[which.min(aics)]
  best <- results[[best_key]]
  cat(sprintf("\nBest model: %s (AIC=%.2f)\n", best_key, best$AIC))
  cat("Parameters:\n")
  for (j in seq_along(best$pars)) {
    cat(sprintf("  %s = %.6f\n", names(best$pars)[j], best$pars[j]))
  }

  # DDD comparison for DD models
  if (grepl("^DD", best_key) && !is.null(ddd_ref)) {
    cat(sprintf("\nDDD reference: lambda=%.4f, mu=%.4f, K=%.1f, loglik=%.2f\n",
                ddd_ref$lambda, ddd_ref$mu, ddd_ref$K, ddd_ref$loglik))
  }
}

# Print all parameter estimates
cat("\n\nAll parameter estimates:\n")
cat(strrep("-", 60), "\n")
for (nm in names(results)) {
  r <- results[[nm]]
  if (!is.null(r$pars)) {
    cat(sprintf("%-10s: %s\n", nm,
                paste(names(r$pars), "=", round(r$pars, 4), collapse=", ")))
  } else {
    cat(sprintf("%-10s: (failed)\n", nm))
  }
}

# ── Save ─────────────────────────────────────────────────────────

save(results, ddd_ref, file = "doc/bird_orders_results.RData")
cat("\nSaved to doc/bird_orders_results.RData\n")

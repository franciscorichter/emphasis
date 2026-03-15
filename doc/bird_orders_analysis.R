#!/usr/bin/env Rscript
# ============================================================================
# Bird phylogenies model comparison — emphasis package
#
# 4 models (CR, DD, PD, EP) x 2 links (linear, exponential) x
# 2 phylogenies (bird.orders 23 tips, bird.families 137 tips)
# = 16 fits via GAM method, ranked by AIC per dataset.
#
# EP + exponential link is not yet supported in inference and
# will be skipped gracefully.
#
# DDD reference (dd_ML, cond=1):
#   bird.orders:   lam0=0.693, mu~0,    K=23
#   bird.families: lam0=0.795, mu=0.102, K=201
#
# Usage:  Rscript bird_orders_analysis.R
# Output: console tables + bird_analysis_results.RData
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

# ── 4 models x 2 links = 8 configurations ───────────────────────

models <- list(
  CR_lin = list(
    model = "cr", link = "linear",
    lb = c(0.01, 0),
    ub = c(2.0, 1.0)
  ),
  CR_exp = list(
    model = "cr", link = "exponential",
    lb = c(-5, -8),
    ub = c(1, 0)
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
  PD_lin = list(
    model = "pd", link = "linear",
    lb = c(0.01, -0.5, 0, -0.5),
    ub = c(2.0, 0.5, 1.0, 0.5)
  ),
  PD_exp = list(
    model = "pd", link = "exponential",
    lb = c(-2, -0.5, -8, -0.5),
    ub = c(2, 0.5, 0, 0.5)
  ),
  EP_lin = list(
    model = "ep", link = "linear",
    lb = c(0.01, -0.2, 0, -0.2),
    ub = c(2.0, 0.5, 1.0, 0.5)
  ),
  EP_exp = list(
    model = "ep", link = "exponential",
    lb = c(-2, -0.5, -8, -0.5),
    ub = c(2, 0.5, 0, 0.5)
  )
)

# ── GAM control ──────────────────────────────────────────────────

gam_ctrl_2p <- list(
  grid_points = 10,
  sample_size = 50
)
gam_ctrl_4p <- list(
  n_grid = 100,
  sample_size = 50
)

# ── Fit all 8 models on both datasets (16 fits) ─────────────────

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
        cat("SKIP:", e$message, "\n")
        NULL
      }
    )
    elapsed <- (proc.time() - t0)[3]

    key <- paste(ds_name, mod_name, sep = "_")
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
      cat(sprintf("  (%.1fs)\n\n", elapsed))
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

# ── Summary tables (per dataset, ranked by AIC) ─────────────────

cat("\n")
cat("============================================\n")
cat("           RESULTS SUMMARY\n")
cat("============================================\n")

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
  valid <- is.finite(aics)
  ord <- order(aics)
  best_aic <- min(aics, na.rm = TRUE)

  hdr <- sprintf(
    "%-10s %-5s %3s %10s %10s %10s %7s",
    "Model", "Link", "nP", "loglik",
    "AIC", "dAIC", "Time"
  )
  cat(hdr, "\n")
  cat(strrep("-", nchar(hdr)), "\n")

  for (i in ord) {
    r <- ds_res[[i]]
    if (is.na(r$AIC)) {
      cat(sprintf(
        "%-10s %-5s %3d %10s %10s %10s %6.1fs\n",
        r$model, r$link, r$n_pars,
        "---", "---", "---", r$elapsed
      ))
    } else {
      cat(sprintf(
        "%-10s %-5s %3d %10.2f %10.2f %10.2f %6.1fs\n",
        r$model, r$link, r$n_pars,
        r$loglik, r$AIC, r$AIC - best_aic,
        r$elapsed
      ))
    }
  }

  # Best model summary
  if (any(valid)) {
    best_key <- names(ds_res)[which.min(aics)]
    best <- ds_res[[best_key]]
    cat(sprintf(
      "\nBest model: %s (link=%s)\n",
      best$model, best$link
    ))
    for (j in seq_along(best$pars)) {
      cat(sprintf(
        "  %s = %.4f\n",
        names(best$pars)[j], best$pars[j]
      ))
    }
  }
}

# ── Save ─────────────────────────────────────────────────────────

save(all_results, file = "bird_analysis_results.RData")
cat("\nSaved to bird_analysis_results.RData\n")

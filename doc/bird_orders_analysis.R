#!/usr/bin/env Rscript
# ============================================================================
# Bird phylogenies model comparison — emphasis package
#
# 3 models (DD, PD, EP) x 2 links (linear, exponential) x
# 2 phylogenies (bird.orders 23 tips, bird.families 137 tips)
# = 12 fits via GAM method, ranked by AIC per dataset.
#
# Bounds are auto-detected via auto_bounds().
# EP + exponential link is not yet supported and will be skipped.
#
# DDD reference (dd_ML, cond=1):
#   bird.orders:   lam0=0.693, mu~0,    K=23
#   bird.families: lam0=0.795, mu=0.102, K=201
#
# Usage:  Rscript bird_orders_analysis.R
# Output: console tables + bird_analysis_results.RData
# ============================================================================

# ── Install dependencies ─────────────────────────────────────────
for (pkg in c("devtools", "ape", "mgcv")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install emphasis from GitHub (always fresh)
if (requireNamespace("emphasis", quietly = TRUE)) {
  try(remove.packages("emphasis"), silent = TRUE)
}
devtools::install_github(
  "franciscorichter/emphasis",
  force = TRUE, upgrade = "never"
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

# ── 3 models x 2 links = 6 configurations ───────────────────────

configs <- list(
  DD_lin = list(model = "dd", link = "linear"),
  DD_exp = list(model = "dd", link = "exponential"),
  PD_lin = list(model = "pd", link = "linear"),
  PD_exp = list(model = "pd", link = "exponential"),
  EP_lin = list(model = "ep", link = "linear"),
  EP_exp = list(model = "ep", link = "exponential")
)

# ── GAM control ──────────────────────────────────────────────────

gam_ctrl <- list(n_grid = 100, sample_size = 50)

# ── Fit all models on both datasets ─────────────────────────────

all_results <- list()
set.seed(2026)

for (ds_name in names(datasets)) {
  tree <- datasets[[ds_name]]
  n_tips <- Ntip(tree)
  cat(sprintf(
    "\n=== %s (%d tips) ===\n\n",
    ds_name, n_tips
  ))

  for (cfg_name in names(configs)) {
    cfg <- configs[[cfg_name]]

    cat(sprintf(
      "  [%s] model=%s link=%s\n",
      cfg_name, cfg$model, cfg$link
    ))

    # Auto-detect bounds
    cat("    detecting bounds... ")
    t0 <- proc.time()
    ab <- tryCatch(
      auto_bounds(
        tree, model = cfg$model, link = cfg$link,
        n_probe = 300, verbose = FALSE
      ),
      error = function(e) {
        cat("SKIP:", e$message, "\n")
        NULL
      }
    )

    if (is.null(ab)) {
      key <- paste(ds_name, cfg_name, sep = "_")
      all_results[[key]] <- list(
        dataset = ds_name, model = cfg_name,
        link = cfg$link, n_tips = n_tips,
        n_pars = NA, pars = NA,
        loglik = NA, AIC = NA, elapsed = 0
      )
      next
    }

    t_bounds <- (proc.time() - t0)[3]
    cat(sprintf(
      "done (%.0fs, %d/%d feasible)\n",
      t_bounds, ab$n_feasible, ab$n_probed
    ))
    cat(sprintf(
      "    lb: %s\n    ub: %s\n",
      paste(round(ab$lower_bound, 3), collapse = ", "),
      paste(round(ab$upper_bound, 3), collapse = ", ")
    ))

    # Fit via GAM
    cat("    fitting GAM... ")
    ctrl <- gam_ctrl
    ctrl$lower_bound <- ab$lower_bound
    ctrl$upper_bound <- ab$upper_bound

    t1 <- proc.time()
    fit <- tryCatch(
      estimate_rates(
        tree,
        method  = "gam",
        model   = cfg$model,
        link    = cfg$link,
        cond    = ab$survival_gam,
        control = ctrl
      ),
      error = function(e) {
        cat("ERROR:", e$message, "\n")
        NULL
      }
    )
    t_fit <- (proc.time() - t1)[3]
    elapsed <- t_bounds + t_fit

    key <- paste(ds_name, cfg_name, sep = "_")
    if (!is.null(fit)) {
      cat(sprintf("done (%.0fs)\n", t_fit))
      pars_str <- paste(
        names(fit$pars), "=",
        round(fit$pars, 4),
        collapse = ", "
      )
      cat(sprintf(
        "    pars: %s\n",
        pars_str
      ))
      cat(sprintf(
        "    loglik=%.2f  AIC=%.2f\n\n",
        fit$loglik, fit$AIC
      ))
      all_results[[key]] <- list(
        dataset = ds_name,
        model   = cfg_name,
        link    = cfg$link,
        n_tips  = n_tips,
        n_pars  = length(fit$pars),
        pars    = fit$pars,
        loglik  = fit$loglik,
        AIC     = fit$AIC,
        elapsed = elapsed
      )
    } else {
      all_results[[key]] <- list(
        dataset = ds_name,
        model   = cfg_name,
        link    = cfg$link,
        n_tips  = n_tips,
        n_pars  = length(ab$lower_bound),
        pars    = rep(NA, length(ab$lower_bound)),
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
        "%-10s %-5s %3s %10s %10s %10s %6.0fs\n",
        r$model, r$link, "---",
        "---", "---", "---", r$elapsed
      ))
    } else {
      cat(sprintf(
        "%-10s %-5s %3d %10.2f %10.2f %10.2f %6.0fs\n",
        r$model, r$link, r$n_pars,
        r$loglik, r$AIC, r$AIC - best_aic,
        r$elapsed
      ))
    }
  }

  # Best model
  valid <- is.finite(aics)
  if (any(valid)) {
    best_key <- names(ds_res)[which.min(aics)]
    best <- ds_res[[best_key]]
    cat(sprintf(
      "\nBest: %s (link=%s)\n",
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

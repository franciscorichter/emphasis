#!/usr/bin/env Rscript
# ============================================================================
# Bird families: 6 models (DD/PD/EP x linear/exponential)
# Pipeline: auto_bounds -> GAM(sample_size=1) -> MCEM(sample_size=1)
# All conditioned on survival. Detailed logging.
# ============================================================================

library(emphasis)
library(ape)

LOG <- function(...) {
  msg <- sprintf(...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush.console()
}

LOG("bird.families analysis — emphasis %s", as.character(packageVersion("emphasis")))

data(bird.families)
tree <- bird.families
LOG("Tree: %d tips, crown age=%.1f", Ntip(tree), max(branching.times(tree)))

configs <- list(
  DD_lin = list(model = "dd", link = "linear"),
  DD_exp = list(model = "dd", link = "exponential"),
  PD_lin = list(model = "pd", link = "linear"),
  PD_exp = list(model = "pd", link = "exponential"),
  EP_lin = list(model = "ep", link = "linear"),
  EP_exp = list(model = "ep", link = "exponential")
)

n_cores <- max(1L, parallel::detectCores() - 1L)
LOG("Threads: %d, Models: %d", n_cores, length(configs))

gam_ctrl  <- list(n_grid = 200, sample_size = 1, verbose = FALSE)
mcem_ctrl <- list(sample_size = 1, maxN = 200, max_iter = 200,
                  tol = 0.02, patience = 5, verbose = FALSE)

results <- list()
set.seed(2026)
t0_total <- proc.time()[3]

for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  LOG("======== %s (model=%s, link=%s) ========", cfg_name, cfg$model, cfg$link)
  t0_model <- proc.time()[3]

  # Stage 1: bounds
  LOG("  Stage 1: auto_bounds")
  t0 <- proc.time()[3]
  ab <- tryCatch(
    auto_bounds(tree, model = cfg$model, link = cfg$link,
                num_threads = n_cores, verbose = FALSE),
    error = function(e) { LOG("    FAILED: %s", e$message); NULL }
  )
  t_bounds <- proc.time()[3] - t0
  if (is.null(ab)) {
    LOG("  Bounds failed (%.0fs) — skipping", t_bounds)
    results[[cfg_name]] <- list(config = cfg, pars = NULL, loglik = NA,
                                AIC = NA, elapsed = t_bounds, status = "bounds_failed")
    next
  }
  LOG("    lb: %s", paste(round(ab$lower_bound, 4), collapse = ", "))
  LOG("    ub: %s", paste(round(ab$upper_bound, 4), collapse = ", "))
  LOG("    Done (%.0fs)", t_bounds)

  # Stage 2: GAM
  LOG("  Stage 2: GAM (n_grid=200, sample_size=1)")
  ctrl_g <- gam_ctrl
  ctrl_g$lower_bound <- ab$lower_bound
  ctrl_g$upper_bound <- ab$upper_bound
  ctrl_g$num_threads <- n_cores

  t0 <- proc.time()[3]
  gam_fit <- tryCatch(
    estimate_rates(tree, method = "gam", model = cfg$model, link = cfg$link,
                   cond = ab$survival_gam, control = ctrl_g),
    error = function(e) { LOG("    GAM FAILED: %s", e$message); NULL }
  )
  t_gam <- proc.time()[3] - t0
  gam_ok <- !is.null(gam_fit) && is.finite(gam_fit$loglik)
  if (gam_ok) {
    LOG("    loglik=%.2f AIC=%.2f (%.0fs)", gam_fit$loglik, gam_fit$AIC, t_gam)
    LOG("    pars: %s", paste(names(gam_fit$pars), "=", round(gam_fit$pars, 5), collapse = ", "))
  } else {
    LOG("    No valid result (%.0fs)", t_gam)
  }

  # Stage 3: MCEM
  LOG("  Stage 3: MCEM (sample_size=1, max_iter=200)")
  ctrl_m <- mcem_ctrl
  ctrl_m$lower_bound <- ab$lower_bound
  ctrl_m$upper_bound <- ab$upper_bound
  ctrl_m$num_threads <- n_cores
  init <- if (gam_ok) gam_fit$pars else NULL

  t0 <- proc.time()[3]
  mcem_fit <- tryCatch(
    estimate_rates(tree, method = "mcem", model = cfg$model, link = cfg$link,
                   init_pars = init, cond = ab$survival_gam, control = ctrl_m),
    error = function(e) { LOG("    MCEM FAILED: %s", e$message); NULL }
  )
  t_mcem <- proc.time()[3] - t0
  mcem_ok <- !is.null(mcem_fit) && is.finite(mcem_fit$loglik)
  if (mcem_ok) {
    LOG("    loglik=%.2f AIC=%.2f (%.0fs)", mcem_fit$loglik, mcem_fit$AIC, t_mcem)
    LOG("    pars: %s", paste(names(mcem_fit$pars), "=", round(mcem_fit$pars, 5), collapse = ", "))
    if (!is.null(mcem_fit$details) && !is.null(mcem_fit$details$mcem)) {
      tr <- mcem_fit$details$mcem
      LOG("    MCEM iterations: %d, stop: %s", nrow(tr), mcem_fit$details$stop_reason)
    }
  } else {
    LOG("    No valid result (%.0fs)", t_mcem)
  }

  # Pick best
  if (mcem_ok) {
    final <- mcem_fit; stage <- "mcem"
  } else if (gam_ok) {
    final <- gam_fit; stage <- "gam"
  } else {
    final <- NULL; stage <- "failed"
  }

  elapsed <- proc.time()[3] - t0_model
  LOG("  RESULT: %s | loglik=%s | AIC=%s | %.0fs",
      stage,
      if (!is.null(final)) sprintf("%.2f", final$loglik) else "NA",
      if (!is.null(final)) sprintf("%.2f", final$AIC) else "NA",
      elapsed)

  results[[cfg_name]] <- list(
    config = cfg, gam_fit = gam_fit, mcem_fit = mcem_fit,
    pars = if (!is.null(final)) final$pars else NULL,
    loglik = if (!is.null(final)) final$loglik else NA,
    AIC = if (!is.null(final)) final$AIC else NA,
    n_pars = if (!is.null(final)) final$n_pars else NA,
    elapsed = elapsed, status = stage,
    t_bounds = t_bounds, t_gam = t_gam, t_mcem = t_mcem
  )
  save(results, file = "doc/bird_families_results.RData")
}

# DDD reference
LOG("======== DDD REFERENCE ========")
ddd_ref <- NULL
if (requireNamespace("DDD", quietly = TRUE)) {
  t0 <- proc.time()[3]
  ddd_fit <- tryCatch(
    DDD::dd_ML(brts = branching.times(tree), cond = 1, btorph = 1,
               initparsopt = c(0.8, 0.1, 200)),
    error = function(e) { LOG("DDD failed: %s", e$message); NULL }
  )
  if (!is.null(ddd_fit)) {
    ddd_ref <- list(lambda = ddd_fit$lambda, mu = ddd_fit$mu, K = ddd_fit$K,
                    loglik = ddd_fit$loglik, AIC = -2 * ddd_fit$loglik + 2 * 3)
    LOG("  lambda=%.4f mu=%.4f K=%.1f loglik=%.2f AIC=%.2f (%.0fs)",
        ddd_ref$lambda, ddd_ref$mu, ddd_ref$K, ddd_ref$loglik, ddd_ref$AIC,
        proc.time()[3] - t0)
  }
} else {
  LOG("DDD not installed")
  ddd_ref <- list(lambda = 0.795, mu = 0.102, K = 201, loglik = -458.8, AIC = 923.6)
}

# Summary
total_time <- proc.time()[3] - t0_total
LOG("")
LOG("============ SUMMARY (%.1f min) ============", total_time / 60)

aics <- sapply(results, function(r) r$AIC)
ord  <- order(aics, na.last = TRUE)
best_aic <- min(aics, na.rm = TRUE)

cat(sprintf("\n%-10s %-6s %3s %8s %10s %10s %8s %5s %5s %5s\n",
            "Model", "Link", "nP", "Status", "loglik", "AIC", "dAIC", "Bnd", "GAM", "MCEM"))
cat(strrep("-", 80), "\n")
for (i in ord) {
  r <- results[[i]]; nm <- names(results)[i]
  if (is.na(r$AIC)) {
    cat(sprintf("%-10s %-6s %3s %8s %10s %10s %8s\n",
                nm, r$config$link, "---", r$status, "---", "---", "---"))
  } else {
    cat(sprintf("%-10s %-6s %3d %8s %10.2f %10.2f %8.2f %4.0fs %4.0fs %4.0fs\n",
                nm, r$config$link, r$n_pars, r$status,
                r$loglik, r$AIC, r$AIC - best_aic,
                r$t_bounds, r$t_gam, r$t_mcem))
  }
}
if (!is.null(ddd_ref))
  cat(sprintf("\n%-10s %-6s %3d %8s %10.2f %10.2f\n",
              "DDD_ref", "---", 3, "analytic", ddd_ref$loglik, ddd_ref$AIC))

cat("\nParameter estimates:\n")
for (nm in names(results)) {
  r <- results[[nm]]
  if (!is.null(r$pars))
    cat(sprintf("  %-10s [%s]: %s\n", nm, r$status,
                paste(names(r$pars), "=", round(r$pars, 5), collapse = ", ")))
}

save(results, ddd_ref, file = "doc/bird_families_results.RData")
LOG("Saved doc/bird_families_results.RData")

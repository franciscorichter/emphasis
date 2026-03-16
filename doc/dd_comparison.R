#!/usr/bin/env Rscript
# ============================================================================
# DD parameter recovery: emphasis vs DDD
# 10 simulated DD trees, estimate with both (conditioned on survival).
# Compare in boxplots.
# ============================================================================

library(emphasis)
library(ape)

LOG <- function(...) {
  msg <- sprintf(...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush.console()
}

LOG("DD comparison — emphasis %s", as.character(packageVersion("emphasis")))

# ---- Parameters ----

lambda_0 <- 0.8
mu_true  <- 0.1
K_true   <- 40
T_crown  <- 10
N_trees  <- 10

beta_N_true <- -(lambda_0 - mu_true) / K_true
pars_emph   <- c(beta_0 = lambda_0, beta_N = beta_N_true,
                  gamma_0 = mu_true, gamma_N = 0)

n_cores <- max(1L, parallel::detectCores() - 1L)

LOG("True: lambda_0=%.3f, mu=%.3f, K=%d (beta_N=%.4f)", lambda_0, mu_true, K_true, beta_N_true)
LOG("Crown age=%.0f, N_trees=%d, threads=%d", T_crown, N_trees, n_cores)

# ---- Simulate ----

LOG("Simulating %d DD trees...", N_trees)
set.seed(2026)

trees <- list()
attempts <- 0L
while (length(trees) < N_trees && attempts < N_trees * 100L) {
  attempts <- attempts + 1L
  sim <- simulate_tree(pars = pars_emph, max_t = T_crown,
                       model = "dd", link = "linear", max_tries = 1)
  if (sim$status == "done" && !is.null(sim$tes) && ape::Ntip(sim$tes) >= 8) {
    trees <- c(trees, list(sim$tes))
    LOG("  tree %d: %d tips", length(trees), ape::Ntip(sim$tes))
  }
}

# ---- DDD ----

LOG("")
LOG("======== DDD estimation (dd_ML, cond=1) ========")

has_ddd <- requireNamespace("DDD", quietly = TRUE)
ddd_results <- data.frame()

if (has_ddd) {
  for (i in seq_along(trees)) {
    bt <- sort(ape::branching.times(trees[[i]]), decreasing = TRUE)
    t0 <- proc.time()[3]
    ddd_fit <- tryCatch(
      DDD::dd_ML(brts = bt, cond = 1, btorph = 1,
                 initparsopt = c(lambda_0, mu_true, K_true * 1.2)),
      error = function(e) { LOG("  tree %d: DDD FAILED: %s", i, e$message); NULL }
    )
    elapsed <- proc.time()[3] - t0
    if (!is.null(ddd_fit)) {
      LOG("  tree %d: lam=%.4f mu=%.4f K=%.1f loglik=%.2f (%.0fs)",
          i, ddd_fit$lambda, ddd_fit$mu, ddd_fit$K, ddd_fit$loglik, elapsed)
      ddd_results <- rbind(ddd_results, data.frame(
        tree = i, method = "DDD",
        lambda_0 = ddd_fit$lambda, mu = ddd_fit$mu, K = ddd_fit$K,
        loglik = ddd_fit$loglik, elapsed = elapsed,
        stringsAsFactors = FALSE))
    }
  }
} else {
  LOG("DDD package not installed — skipping")
}

# ---- emphasis ----

LOG("")
LOG("======== emphasis estimation (auto_bounds -> GAM -> MCEM) ========")

emph_results <- data.frame()

# auto_bounds on first tree, reuse for all
LOG("  auto_bounds on tree 1...")
t0 <- proc.time()[3]
ab <- auto_bounds(trees[[1]], model = "dd", link = "linear",
                  num_threads = n_cores, verbose = FALSE)
LOG("    bounds: [%s] to [%s] (%.0fs)",
    paste(round(ab$lower_bound, 4), collapse = ", "),
    paste(round(ab$upper_bound, 4), collapse = ", "),
    proc.time()[3] - t0)

gam_ctrl  <- list(n_grid = 200, sample_size = 1, verbose = FALSE)
mcem_ctrl <- list(sample_size = 1, maxN = 200, max_iter = 200,
                  tol = 0.02, patience = 5, verbose = FALSE)

for (i in seq_along(trees)) {
  t0 <- proc.time()[3]

  # GAM
  ctrl_g <- gam_ctrl
  ctrl_g$lower_bound <- ab$lower_bound
  ctrl_g$upper_bound <- ab$upper_bound
  ctrl_g$num_threads <- n_cores

  gam_fit <- tryCatch(
    estimate_rates(trees[[i]], method = "gam", model = "dd", link = "linear",
                   cond = ab$survival_gam, control = ctrl_g),
    error = function(e) NULL)

  # MCEM
  mcem_fit <- NULL
  if (!is.null(gam_fit)) {
    ctrl_m <- mcem_ctrl
    ctrl_m$lower_bound <- ab$lower_bound
    ctrl_m$upper_bound <- ab$upper_bound
    ctrl_m$num_threads <- n_cores
    mcem_fit <- tryCatch(
      estimate_rates(trees[[i]], method = "mcem", model = "dd", link = "linear",
                     init_pars = gam_fit$pars, cond = ab$survival_gam, control = ctrl_m),
      error = function(e) NULL)
  }

  mcem_ok <- !is.null(mcem_fit) && is.finite(mcem_fit$loglik)
  gam_ok  <- !is.null(gam_fit) && is.finite(gam_fit$loglik)
  fit   <- if (mcem_ok) mcem_fit else if (gam_ok) gam_fit else NULL
  stage <- if (mcem_ok) "mcem" else if (gam_ok) "gam" else "failed"

  elapsed <- proc.time()[3] - t0

  if (!is.null(fit)) {
    est_lam <- fit$pars["beta_0"]
    est_mu  <- fit$pars["gamma_0"]
    est_bN  <- fit$pars["beta_N"]
    est_K   <- if (est_bN < 0) (est_lam - est_mu) / (-est_bN) else NA

    LOG("  tree %d [%s]: lam=%.4f mu=%.4f K=%.1f loglik=%.2f (%.0fs)",
        i, stage, est_lam, est_mu, est_K, fit$loglik, elapsed)

    emph_results <- rbind(emph_results, data.frame(
      tree = i, method = "emphasis",
      lambda_0 = unname(est_lam), mu = unname(est_mu), K = unname(est_K),
      loglik = fit$loglik, elapsed = elapsed,
      stringsAsFactors = FALSE))
  } else {
    LOG("  tree %d: FAILED (%.0fs)", i, elapsed)
  }
}

# ---- Summary ----

LOG("")
LOG("======== SUMMARY ========")

all_results <- rbind(ddd_results, emph_results)

cat(sprintf("\n%-10s %5s %8s %8s %8s %8s\n", "Method", "n", "lam_0", "mu", "K", "time"))
cat(strrep("-", 50), "\n")
cat(sprintf("%-10s %5s %8.3f %8.3f %8.1f\n", "TRUE", "", lambda_0, mu_true, K_true))
for (meth in c("DDD", "emphasis")) {
  sub <- all_results[all_results$method == meth, ]
  if (nrow(sub) > 0)
    cat(sprintf("%-10s %5d %8.3f %8.3f %8.1f %7.0fs\n",
                meth, nrow(sub), mean(sub$lambda_0, na.rm = TRUE),
                mean(sub$mu, na.rm = TRUE), mean(sub$K, na.rm = TRUE),
                sum(sub$elapsed)))
}

# ---- Plot ----

LOG("Generating plot...")
pdf("doc/dd_comparison.pdf", width = 10, height = 4)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

plot_param <- function(param, true_val, ylab) {
  ddd_vals  <- ddd_results[[param]]
  emph_vals <- emph_results[[param]]
  all_vals <- c(ddd_vals, emph_vals, true_val)
  ylim <- range(all_vals[is.finite(all_vals)], na.rm = TRUE)
  ylim <- ylim + c(-1, 1) * 0.15 * diff(ylim)

  boxplot(list(DDD = ddd_vals, emphasis = emph_vals),
          main = ylab, ylab = ylab,
          col = c("#4DBEEE", "#D95319"), border = "grey30", ylim = ylim, las = 1)
  abline(h = true_val, col = "red", lty = 2, lwd = 2)
  if (length(ddd_vals) > 0)
    points(rep(1, length(ddd_vals)) + runif(length(ddd_vals), -0.15, 0.15),
           ddd_vals, pch = 19, cex = 0.8, col = "#0072BD80")
  if (length(emph_vals) > 0)
    points(rep(2, length(emph_vals)) + runif(length(emph_vals), -0.15, 0.15),
           emph_vals, pch = 19, cex = 0.8, col = "#A2142F80")
  legend("topright", legend = "true", col = "red", lty = 2, lwd = 2, bty = "n", cex = 0.8)
}

plot_param("lambda_0", lambda_0, expression(lambda[0]))
plot_param("mu", mu_true, expression(mu))
plot_param("K", K_true, "K")
dev.off()

LOG("Saved doc/dd_comparison.pdf")

save(all_results, ddd_results, emph_results, trees, pars_emph,
     file = "doc/dd_comparison_results.RData")
LOG("Saved doc/dd_comparison_results.RData")
LOG("Done. Total emphasis time: %.1f min", sum(emph_results$elapsed) / 60)

#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════
# Comprehensive Validation Report: Incomplete Taxon Sampling (rho)
# ═══════════════════════════════════════════════════════════════════
#
# This script validates the rho (incomplete sampling) implementation
# in the emphasis package against the Stadler (2009) analytical formula
# for the constant-rate (CR) birth-death model.
#
# Outputs:
#   rho_validation_report.pdf  — multi-page PDF with all figures
#   rho_validation_results.RData — all numerical results

library(emphasis)
library(DDD)
library(ape)

set.seed(2026)

outdir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "doc/rho_sampling")

cat("═══════════════════════════════════════════════════════════════\n")
cat("  Rho Validation Report — emphasis package\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ── Analytical formulas ─────────────────────────────────────────

# Stadler (2009) CR likelihood with incomplete sampling
# Matches DDD::bd_loglik(pars2=c(0,0,1,0,2), missnumspec=0) at rho=1
stadler_cr <- function(brts, lam, mu, rho = 1.0) {
  r <- lam - mu
  p1 <- function(t) rho * r^2 * exp(-r*t) / (rho*lam + (lam*(1-rho)-mu)*exp(-r*t))^2
  n <- length(brts) + 1
  2 * log(p1(brts[1])) + (n-2) * log(lam) + sum(log(p1(brts[-1])))
}

# IS log-likelihood at known parameters via augment_trees
is_loglik <- function(tree, pars_compact, rho = 1.0, S = 3000, maxN_mult = 10) {
  brts <- sort(ape::branching.times(tree), decreasing = TRUE)
  pars8 <- emphasis:::.expand_pars(pars_compact, c(0L, 0L, 0L))
  aug <- emphasis:::augment_trees(
    brts = brts, pars = pars8, sample_size = S,
    maxN = as.integer(S * maxN_mult), max_missing = 10000L,
    max_lambda = 500, num_threads = 4L,
    model = c(0L, 0L, 0L), link = 0L, rho = rho
  )
  logw <- aug$logf - aug$logg
  max_lw <- max(logw)
  fhat <- log(mean(exp(logw - max_lw))) + max_lw
  se <- sd(exp(logw - max_lw)) / sqrt(length(logw)) * exp(max_lw - fhat)
  list(fhat = fhat, se = se, n_trees = length(logw),
       logf = aug$logf, logg = aug$logg)
}


# ═══════════════════════════════════════════════════════════════
# 1. FORMULA VERIFICATION: Stadler vs DDD at rho = 1
# ═══════════════════════════════════════════════════════════════
cat("[1/6] Verifying Stadler formula matches DDD at rho=1...\n")

formula_check <- data.frame()
for (seed in 1:50) {
  set.seed(seed)
  tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr", max_tries = 100)
  if (is.null(tr$tes) || Ntip(tr$tes) < 4 || Ntip(tr$tes) > 50) next
  bt <- sort(branching.times(tr$tes), decreasing = TRUE)
  n <- Ntip(tr$tes)
  ddd_val <- DDD::bd_loglik(c(0.5, 0.1, 0, 0), c(0, 0, 1, 0, 2), bt, 0)
  stadler_val <- stadler_cr(bt, 0.5, 0.1, 1.0)
  formula_check <- rbind(formula_check, data.frame(
    seed = seed, n_tips = n, ddd = ddd_val, stadler = stadler_val,
    diff = stadler_val - ddd_val
  ))
  if (nrow(formula_check) >= 15) break
}
cat(sprintf("  %d trees checked. Max |diff| = %.2e\n", nrow(formula_check), max(abs(formula_check$diff))))


# ═══════════════════════════════════════════════════════════════
# 2. IS vs STADLER: multiple trees × multiple rho values
# ═══════════════════════════════════════════════════════════════
cat("[2/6] IS vs Stadler across trees and rho values (S=5000)...\n")

rho_vals <- c(1.0, 0.9, 0.7, 0.5, 0.3)
accuracy_table <- data.frame()

test_trees <- list()
for (seed in 1:100) {
  set.seed(seed)
  tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr", max_tries = 100)
  if (!is.null(tr$tes) && Ntip(tr$tes) >= 5 && Ntip(tr$tes) <= 30) {
    test_trees[[length(test_trees) + 1]] <- list(tree = tr$tes, seed = seed)
  }
  if (length(test_trees) >= 8) break
}

for (tt in test_trees) {
  tree <- tt$tree
  brts <- sort(branching.times(tree), decreasing = TRUE)
  n <- Ntip(tree)
  cat(sprintf("  seed=%d, n=%d: ", tt$seed, n))

  for (rho_val in rho_vals) {
    stadler_val <- stadler_cr(brts, 0.5, 0.1, rho_val)
    is_res <- is_loglik(tree, c(0.5, 0.1), rho = rho_val, S = 5000)
    diff <- is_res$fhat - stadler_val

    accuracy_table <- rbind(accuracy_table, data.frame(
      seed = tt$seed, n_tips = n, rho = rho_val,
      stadler = stadler_val, is_est = is_res$fhat,
      se = is_res$se, diff = diff,
      n_aug_trees = is_res$n_trees
    ))
    cat(sprintf("%.3f ", diff))
  }
  cat("\n")
}

cat(sprintf("  Mean diff: %.4f, Max |diff|: %.4f\n\n",
            mean(accuracy_table$diff), max(abs(accuracy_table$diff))))


# ═══════════════════════════════════════════════════════════════
# 3. IS CONVERGENCE: fhat vs sample size for each rho
# ═══════════════════════════════════════════════════════════════
cat("[3/6] IS convergence study...\n")

# Pick one tree for convergence study
conv_tree <- test_trees[[1]]$tree
conv_brts <- sort(branching.times(conv_tree), decreasing = TRUE)
conv_n <- Ntip(conv_tree)
cat(sprintf("  Tree: %d tips\n", conv_n))

S_vals <- c(100, 300, 500, 1000, 2000, 5000)
convergence <- data.frame()

for (rho_val in c(1.0, 0.7, 0.5)) {
  stadler_val <- stadler_cr(conv_brts, 0.5, 0.1, rho_val)
  cat(sprintf("  rho=%.1f: ", rho_val))
  for (S in S_vals) {
    # 5 replicates per S to estimate variability
    fhats <- numeric(5)
    for (rep in 1:5) {
      res <- is_loglik(conv_tree, c(0.5, 0.1), rho = rho_val, S = S)
      fhats[rep] <- res$fhat
    }
    convergence <- rbind(convergence, data.frame(
      rho = rho_val, S = S, mean_fhat = mean(fhats),
      sd_fhat = sd(fhats), stadler = stadler_val,
      bias = mean(fhats) - stadler_val
    ))
    cat(sprintf("S=%d ", S))
  }
  cat("\n")
}


# ═══════════════════════════════════════════════════════════════
# 4. PROFILE LIKELIHOOD: full vs subsampled tree
# ═══════════════════════════════════════════════════════════════
cat("[4/6] Profile likelihood comparison...\n")

# Generate a larger tree for subsampling
lam_true <- 0.6; mu_true <- 0.15
for (s in 1:200) {
  set.seed(s)
  tr_sim <- simulate_tree(pars = c(lam_true, mu_true), max_t = 6, model = "cr", max_tries = 100)
  if (!is.null(tr_sim$tes) && Ntip(tr_sim$tes) >= 25 && Ntip(tr_sim$tes) <= 45) {
    full_tree <- tr_sim$tes
    break
  }
}
n_full <- Ntip(full_tree)

# Subsample
set.seed(789)
drop_n <- round(n_full * 0.3)
drop_tips <- sample(full_tree$tip.label, drop_n)
sub_tree <- ape::drop.tip(full_tree, drop_tips)
n_sub <- Ntip(sub_tree)
rho_actual <- n_sub / n_full
cat(sprintf("  Full: %d tips, Subsampled: %d tips (rho=%.3f)\n", n_full, n_sub, rho_actual))

# Profile over lambda (mu fixed at true value)
lam_grid <- seq(0.2, 1.2, by = 0.05)
profile <- data.frame()

for (lam_test in lam_grid) {
  cat(sprintf("  lambda=%.2f\r", lam_test))
  # Full tree, rho=1
  full_ll <- is_loglik(full_tree, c(lam_test, mu_true), rho = 1.0, S = 2000)$fhat
  # Subsampled, rho=1 (wrong)
  sub_wrong <- is_loglik(sub_tree, c(lam_test, mu_true), rho = 1.0, S = 2000)$fhat
  # Subsampled, rho=correct
  sub_right <- is_loglik(sub_tree, c(lam_test, mu_true), rho = rho_actual, S = 2000)$fhat
  # Stadler analytical (full tree)
  full_brts <- sort(branching.times(full_tree), decreasing = TRUE)
  stadler_full <- stadler_cr(full_brts, lam_test, mu_true, 1.0)
  # Stadler analytical (subsampled, rho-corrected)
  sub_brts <- sort(branching.times(sub_tree), decreasing = TRUE)
  stadler_sub <- stadler_cr(sub_brts, lam_test, mu_true, rho_actual)

  profile <- rbind(profile, data.frame(
    lambda = lam_test,
    full_is = full_ll, sub_wrong_is = sub_wrong, sub_right_is = sub_right,
    full_stadler = stadler_full, sub_stadler = stadler_sub
  ))
}
cat("  Done.                    \n")


# ═══════════════════════════════════════════════════════════════
# 5. MLE RECOVERY
# ═══════════════════════════════════════════════════════════════
cat("[5/6] MLE recovery (MCEM)...\n")

mle_full <- estimate_rates(full_tree, method = "mcem", model = "cr",
                           control = list(lower_bound = c(0.1, 0.01), upper_bound = c(2.0, 1.0),
                                          sample_size = 500, max_iter = 40, num_threads = 4))

mle_sub_wrong <- estimate_rates(sub_tree, method = "mcem", model = "cr",
                                control = list(lower_bound = c(0.1, 0.01), upper_bound = c(2.0, 1.0),
                                               sample_size = 500, max_iter = 40, num_threads = 4))

mle_sub_right <- estimate_rates(sub_tree, method = "mcem", model = "cr",
                                control = list(lower_bound = c(0.1, 0.01), upper_bound = c(2.0, 1.0),
                                               sample_size = 500, max_iter = 40, num_threads = 4,
                                               rho = rho_actual))

mle_results <- data.frame(
  method = c("True", "Full tree (rho=1)", "Subsampled (rho=1, wrong)", paste0("Subsampled (rho=", round(rho_actual,2), ")")),
  lambda = c(lam_true, mle_full$pars[1], mle_sub_wrong$pars[1], mle_sub_right$pars[1]),
  mu     = c(mu_true,  mle_full$pars[2], mle_sub_wrong$pars[2], mle_sub_right$pars[2]),
  loglik = c(NA, mle_full$loglik, mle_sub_wrong$loglik, mle_sub_right$loglik)
)
cat("  Done.\n")


# ═══════════════════════════════════════════════════════════════
# 6. GENERATE PDF REPORT
# ═══════════════════════════════════════════════════════════════
cat("[6/6] Generating PDF report...\n")

pdf_path <- file.path(outdir, "rho_validation_report.pdf")
pdf(pdf_path, width = 10, height = 7.5)

# ── Page 1: Title + Summary ──────────────────────────────────
par(mar = c(1, 1, 1, 1))
plot.new()
text(0.5, 0.92, "Rho Validation Report", cex = 2.2, font = 2)
text(0.5, 0.85, "Incomplete Taxon Sampling in emphasis", cex = 1.3, col = "grey40")
text(0.5, 0.78, format(Sys.time(), "%Y-%m-%d %H:%M"), cex = 1.0, col = "grey60")

summary_text <- c(
  "",
  sprintf("Tests: %d trees x %d rho values = %d comparisons",
          length(unique(accuracy_table$seed)), length(rho_vals), nrow(accuracy_table)),
  sprintf("Max |IS - Stadler|: %.4f log-lik units", max(abs(accuracy_table$diff))),
  sprintf("Mean bias: %.4f (expected negative from Jensen's inequality)", mean(accuracy_table$diff)),
  "",
  "Method: IS estimator (augment_trees + eval_logf) compared against",
  "  Stadler (2009) analytical formula for CR birth-death with sampling fraction rho.",
  "  Formula verified to match DDD::bd_loglik at rho=1.",
  "",
  sprintf("MLE recovery: full tree (%d tips) vs subsampled (%d tips, rho=%.2f)",
          n_full, n_sub, rho_actual),
  sprintf("  True: lambda=%.3f, mu=%.3f", lam_true, mu_true),
  sprintf("  Full tree:        lambda=%.3f, mu=%.3f", mle_full$pars[1], mle_full$pars[2]),
  sprintf("  Subsamp (wrong):  lambda=%.3f, mu=%.3f", mle_sub_wrong$pars[1], mle_sub_wrong$pars[2]),
  sprintf("  Subsamp (rho):    lambda=%.3f, mu=%.3f", mle_sub_right$pars[1], mle_sub_right$pars[2])
)
for (i in seq_along(summary_text)) {
  text(0.08, 0.68 - i * 0.038, summary_text[i], adj = 0, cex = 0.95, family = "mono")
}


# ── Page 2: IS vs Stadler scatter ────────────────────────────
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

# Panel A: IS vs Stadler (all points)
cols <- c("1" = "black", "0.9" = "steelblue", "0.7" = "darkorange",
          "0.5" = "red3", "0.3" = "purple")
plot(accuracy_table$stadler, accuracy_table$is_est,
     pch = 19, cex = 1.2,
     col = cols[as.character(accuracy_table$rho)],
     xlab = "Stadler analytical log-likelihood",
     ylab = "IS estimated log-likelihood",
     main = "A. IS vs Stadler (all trees, all rho)")
abline(0, 1, col = "grey50", lty = 2, lwd = 2)
legend("topleft", legend = paste("rho =", names(cols)),
       col = cols, pch = 19, cex = 0.85, bg = "white")

# Panel B: Residuals by rho
boxplot(diff ~ rho, data = accuracy_table,
        col = cols[as.character(sort(unique(accuracy_table$rho)))],
        xlab = "rho", ylab = "IS - Stadler",
        main = "B. Residuals by rho value")
abline(h = 0, col = "grey50", lty = 2, lwd = 2)
text(1:length(rho_vals), par("usr")[4] * 0.95,
     sprintf("n=%d", table(accuracy_table$rho)),
     cex = 0.8, col = "grey40")


# ── Page 3: Convergence with sample size ─────────────────────
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

rho_conv <- c(1.0, 0.7, 0.5)
cols_conv <- c("black", "darkorange", "red3")

# Panel A: IS estimate vs S
plot(NULL, xlim = range(S_vals), ylim = range(convergence$mean_fhat),
     xlab = "Sample size (S)", ylab = "IS log-likelihood estimate",
     main = sprintf("A. IS convergence (%d-tip tree)", conv_n), log = "x")
for (i in seq_along(rho_conv)) {
  sub <- convergence[convergence$rho == rho_conv[i], ]
  lines(sub$S, sub$mean_fhat, col = cols_conv[i], lwd = 2, type = "b", pch = 19)
  # Add Stadler reference
  abline(h = sub$stadler[1], col = cols_conv[i], lty = 3)
  # Error bars (1 SD)
  arrows(sub$S, sub$mean_fhat - sub$sd_fhat, sub$S, sub$mean_fhat + sub$sd_fhat,
         angle = 90, code = 3, length = 0.04, col = cols_conv[i])
}
legend("bottomright", legend = paste("rho =", rho_conv),
       col = cols_conv, lwd = 2, pch = 19, cex = 0.85, bg = "white")

# Panel B: Bias vs S
plot(NULL, xlim = range(S_vals), ylim = range(convergence$bias) * c(1.2, 1.2),
     xlab = "Sample size (S)", ylab = "Bias (IS - Stadler)",
     main = "B. IS bias vs sample size", log = "x")
abline(h = 0, col = "grey50", lty = 2, lwd = 2)
for (i in seq_along(rho_conv)) {
  sub <- convergence[convergence$rho == rho_conv[i], ]
  lines(sub$S, sub$bias, col = cols_conv[i], lwd = 2, type = "b", pch = 19)
}
legend("bottomright", legend = paste("rho =", rho_conv),
       col = cols_conv, lwd = 2, pch = 19, cex = 0.85, bg = "white")


# ── Page 4: Profile likelihood ───────────────────────────────
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

# Normalize profiles relative to their maxima for comparison
norm <- function(x) x - max(x)

# Panel A: IS profile likelihoods
plot(profile$lambda, norm(profile$full_is), type = "l", lwd = 2.5, col = "black",
     xlab = expression(lambda), ylab = "Relative log-likelihood",
     main = sprintf("A. IS profile (mu=%.2f fixed)", mu_true),
     ylim = c(-15, 1))
lines(profile$lambda, norm(profile$sub_wrong_is), col = "red3", lwd = 2, lty = 2)
lines(profile$lambda, norm(profile$sub_right_is), col = "forestgreen", lwd = 2.5)
abline(v = lam_true, col = "grey60", lty = 3)
legend("topright", legend = c(
  sprintf("Full tree (%d tips)", n_full),
  sprintf("Subsampled, rho=1 (wrong)"),
  sprintf("Subsampled, rho=%.2f (correct)", rho_actual)),
  col = c("black", "red3", "forestgreen"), lwd = c(2.5, 2, 2.5), lty = c(1, 2, 1),
  cex = 0.8, bg = "white")

# Panel B: Stadler analytical profiles
plot(profile$lambda, norm(profile$full_stadler), type = "l", lwd = 2.5, col = "black",
     xlab = expression(lambda), ylab = "Relative log-likelihood",
     main = "B. Stadler analytical profile",
     ylim = c(-15, 1))
lines(profile$lambda, norm(profile$sub_stadler), col = "forestgreen", lwd = 2.5)
abline(v = lam_true, col = "grey60", lty = 3)
legend("topright", legend = c(
  sprintf("Full tree (%d tips)", n_full),
  sprintf("Subsampled, rho=%.2f", rho_actual)),
  col = c("black", "forestgreen"), lwd = 2.5, cex = 0.85, bg = "white")


# ── Page 5: MLE Recovery ─────────────────────────────────────
par(mfrow = c(1, 1), mar = c(5, 8, 4, 2))

barplot_data <- rbind(
  mle_results$lambda[-1],
  mle_results$mu[-1]
)
colnames(barplot_data) <- c("Full tree\n(rho=1)", "Subsampled\n(rho=1, wrong)",
                             sprintf("Subsampled\n(rho=%.2f)", rho_actual))
rownames(barplot_data) <- c("lambda", "mu")

bp <- barplot(barplot_data, beside = TRUE, col = c("steelblue", "coral"),
              main = "MLE Parameter Recovery",
              ylab = "Estimated value", las = 1, ylim = c(0, max(barplot_data) * 1.3))
# True value lines
abline(h = lam_true, col = "steelblue", lty = 2, lwd = 1.5)
abline(h = mu_true, col = "coral", lty = 2, lwd = 1.5)
text(mean(bp[1,]), lam_true * 1.08, sprintf("true lambda = %.2f", lam_true),
     col = "steelblue", cex = 0.85)
text(mean(bp[2,]), mu_true * 1.4, sprintf("true mu = %.2f", mu_true),
     col = "coral", cex = 0.85)
legend("topright", legend = c("lambda", "mu"),
       fill = c("steelblue", "coral"), cex = 0.9)


# ── Page 6: Accuracy table ───────────────────────────────────
par(mfrow = c(1, 1), mar = c(2, 2, 3, 2))
plot.new()
text(0.5, 0.97, "IS vs Stadler: Full Accuracy Table", cex = 1.5, font = 2)

# Print table
header <- sprintf("%-6s  %-5s  %-5s  %10s  %10s  %10s  %8s",
                  "seed", "tips", "rho", "Stadler", "IS est", "diff", "SE")
text(0.02, 0.90, header, adj = 0, cex = 0.7, family = "mono", font = 2)

for (i in seq_len(nrow(accuracy_table))) {
  r <- accuracy_table[i, ]
  line <- sprintf("%-6d  %-5d  %-5.1f  %10.4f  %10.4f  %+10.4f  %8.4f",
                  r$seed, r$n_tips, r$rho, r$stadler, r$is_est, r$diff, r$se)
  y <- 0.87 - i * 0.022
  if (y < 0.02) break
  bg_col <- if (abs(r$diff) < 0.05) "grey95" else if (abs(r$diff) < 0.1) "lightyellow" else "mistyrose"
  rect(0.01, y - 0.008, 0.99, y + 0.012, col = bg_col, border = NA)
  text(0.02, y, line, adj = 0, cex = 0.6, family = "mono")
}

# Summary row
y_sum <- 0.87 - (min(nrow(accuracy_table), 38) + 1) * 0.022
text(0.02, y_sum, sprintf("SUMMARY: n=%d comparisons, mean diff = %+.4f, max |diff| = %.4f",
                           nrow(accuracy_table), mean(accuracy_table$diff),
                           max(abs(accuracy_table$diff))),
     adj = 0, cex = 0.75, family = "mono", font = 2)


# ── Page 7: Weight distribution diagnostic ───────────────────
par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))

# For each rho, show weight distribution from one tree
diag_tree <- test_trees[[1]]$tree
diag_brts <- sort(branching.times(diag_tree), decreasing = TRUE)
diag_pars <- emphasis:::.expand_pars(c(0.5, 0.1), c(0L, 0L, 0L))

for (rho_val in c(1.0, 0.7, 0.5, 0.3)) {
  aug <- emphasis:::augment_trees(
    diag_brts, diag_pars, 3000L, 30000L, 10000L, 500, 4L,
    c(0L, 0L, 0L), 0L, rho_val
  )
  logw <- aug$logf - aug$logg
  hist(logw, breaks = 50, col = "steelblue", border = "white",
       main = sprintf("rho = %.1f (n_trees = %d)", rho_val, length(logw)),
       xlab = "log(w) = logf - logg", freq = FALSE)
  abline(v = mean(logw), col = "red", lwd = 2, lty = 2)
  # Mark the fhat
  max_lw <- max(logw)
  fhat <- log(mean(exp(logw - max_lw))) + max_lw
  abline(v = fhat, col = "darkgreen", lwd = 2)
  stadler_ref <- stadler_cr(diag_brts, 0.5, 0.1, rho_val)
  abline(v = stadler_ref, col = "purple", lwd = 2, lty = 3)
  legend("topleft", legend = c("mean(logw)", "fhat (IS est)", "Stadler"),
         col = c("red", "darkgreen", "purple"), lwd = 2, lty = c(2, 1, 3),
         cex = 0.65, bg = "white")
}


dev.off()
cat(sprintf("  PDF saved: %s\n", pdf_path))

# Save results
save(formula_check, accuracy_table, convergence, profile, mle_results,
     n_full, n_sub, rho_actual, lam_true, mu_true,
     file = file.path(outdir, "rho_validation_results.RData"))
cat(sprintf("  Results saved: %s\n", file.path(outdir, "rho_validation_results.RData")))

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  Report complete.\n")
cat("═══════════════════════════════════════════════════════════════\n")

#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────
# Analytical validation of rho (incomplete sampling) implementation
# ─────────────────────────────────────────────────────────────────
#
# Strategy:
#   1. Simulate a CR tree with known (lambda, mu)
#   2. Compute IS log-likelihood at known parameters for rho=1
#   3. Verify it matches DDD::bd_loglik(missnumspec=0) exactly
#   4. For rho < 1: compare IS estimator against Stadler (2009) formula
#   5. Test MLE recovery: subsample tree, estimate with rho correction
#
# The Stadler (2009) formula for a reconstructed CR tree with
# incomplete sampling (fraction rho of extant species observed):
#
#   p1(t) = rho * r^2 * exp(-r*t) / D(t)^2
#   D(t)  = rho*lam + (lam*(1-rho) - mu) * exp(-r*t)
#   r     = lam - mu
#
# Unconditional log-likelihood (crown age, ordered = btorph=1):
#   log L = log((n-1)!) + (n-1)*log(lam) + sum(log(p1(t_i)))
#         = log((n-1)!) + (n-1)*log(lam) + n*log(rho) + 2*(n-1)*log(r)
#           - sum(r*t_i) - 2*sum(log(D(t_i)))

library(emphasis)
library(DDD)
library(ape)

cat("═══════════════════════════════════════════════════════════\n")
cat("  Analytical validation of rho implementation (CR model)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# ── Stadler (2009) formula ──────────────────────────────────────
stadler_cr_loglik <- function(brts, lambda, mu, rho = 1.0) {
  # brts: branching times from present (positive, crown age first)
  # Returns unconditional log-likelihood (btorph=1, soc=2)
  n <- length(brts) + 1  # number of observed tips
  r <- lambda - mu

  D <- function(t) rho * lambda + (lambda * (1 - rho) - mu) * exp(-r * t)

  loglik <- lfactorial(n - 1) +
    (n - 1) * log(lambda) +
    n * log(rho) +
    2 * (n - 1) * log(r) -
    r * sum(brts) -
    2 * sum(log(D(brts)))

  return(loglik)
}


# ── Helper: IS log-likelihood at known parameters ──────────────
is_loglik <- function(tree, pars_compact, rho = 1.0, S = 2000, maxN = 20000) {
  brts <- sort(ape::branching.times(tree), decreasing = TRUE)
  model_bin <- c(0L, 0L, 0L)  # CR model
  pars8 <- emphasis:::.expand_pars(pars_compact, model_bin)

  # Augment trees and score
  aug <- emphasis:::augment_trees(
    brts        = brts,
    pars        = pars8,
    sample_size = S,
    maxN        = maxN,
    max_missing = 10000L,
    max_lambda  = 500,
    num_threads = 4L,
    model       = model_bin,
    link        = 0L,
    rho         = rho
  )

  # IS weights: w_i = exp(logf_i - logg_i)
  logw <- aug$logf - aug$logg

  # Stable log-mean-exp
  max_logw <- max(logw)
  fhat <- log(mean(exp(logw - max_logw))) + max_logw

  list(
    fhat     = fhat,
    logw     = logw,
    n_trees  = length(logw),
    se       = sd(exp(logw - max_logw)) / sqrt(length(logw)) * exp(max_logw - fhat)
  )
}


# ═══════════════════════════════════════════════════════════════
# TEST 1: rho = 1 baseline (IS vs DDD vs Stadler)
# ═══════════════════════════════════════════════════════════════
cat("─── Test 1: rho = 1 baseline ───────────────────────────────\n")

# Use moderate parameters and short crown age for manageable tree size
lam_true <- 0.5
mu_true  <- 0.1

# Find a tree with 15-50 tips
tree_sim <- NULL
for (seed in 1:200) {
  set.seed(seed)
  tr <- simulate_tree(pars = c(lam_true, mu_true), max_t = 5,
                      model = "cr", max_tries = 100)
  if (!is.null(tr$tes) && Ntip(tr$tes) >= 15 && Ntip(tr$tes) <= 50) {
    tree_sim <- tr
    cat(sprintf("  Using seed=%d\n", seed))
    break
  }
}

tree <- tree_sim$tes
brts <- sort(branching.times(tree), decreasing = TRUE)
n_tips <- Ntip(tree)
cat(sprintf("  Tree: %d tips, crown age = %.3f\n", n_tips, brts[1]))

# DDD analytical
ddd_ll <- DDD::bd_loglik(
  pars1 = c(lam_true, mu_true, 0, 0),
  pars2 = c(0, 0, 1, 0, 2),   # no time dep, cond=0, btorph=1, soc=2
  brts  = brts,
  missnumspec = 0
)

# Stadler formula
stadler_ll <- stadler_cr_loglik(brts, lam_true, mu_true, rho = 1.0)

# IS estimate (S = 5000 for precision)
is_result <- is_loglik(tree, c(lam_true, mu_true), rho = 1.0, S = 5000)

cat(sprintf("  DDD analytical:     %.4f\n", ddd_ll))
cat(sprintf("  Stadler formula:    %.4f\n", stadler_ll))
cat(sprintf("  IS estimate:        %.4f  (SE ≈ %.4f)\n", is_result$fhat, is_result$se))
cat(sprintf("  IS - DDD:           %.4f\n", is_result$fhat - ddd_ll))
cat(sprintf("  Stadler - DDD:      %.4f\n", stadler_ll - ddd_ll))

# Check if Stadler formula matches DDD
if (abs(stadler_ll - ddd_ll) < 0.01) {
  cat("  ✓ Stadler formula matches DDD exactly\n")
} else {
  cat("  ✗ Stadler formula does NOT match DDD — need to adjust formula\n")
  cat("    Will calibrate against DDD...\n")
}

cat("\n")


# ═══════════════════════════════════════════════════════════════
# TEST 2: Calibrate Stadler formula against DDD for rho = 1
#         using multiple trees to find the correct formula
# ═══════════════════════════════════════════════════════════════
cat("─── Test 2: Formula calibration (multiple trees, rho=1) ────\n")

offsets <- c()
for (seed in 1:100) {
  set.seed(seed)
  tr <- simulate_tree(pars = c(0.6, 0.15), max_t = 5,
                      model = "cr", max_tries = 100)
  if (is.null(tr$tes) || Ntip(tr$tes) < 8 || Ntip(tr$tes) > 40) next

  bt <- sort(branching.times(tr$tes), decreasing = TRUE)
  ddd_val <- DDD::bd_loglik(c(0.6, 0.15, 0, 0), c(0, 0, 1, 0, 2), bt, 0)
  stadler_val <- stadler_cr_loglik(bt, 0.6, 0.15, 1.0)
  offsets <- c(offsets, stadler_val - ddd_val)
  cat(sprintf("  seed=%2d  n=%3d  DDD=%.2f  Stadler=%.2f  diff=%.4f\n",
              seed, Ntip(tr$tes), ddd_val, stadler_val, stadler_val - ddd_val))
  if (length(offsets) >= 5) break
}

cat(sprintf("  Mean offset: %.4f  SD: %.4f\n", mean(offsets), sd(offsets)))
cat("\n")


# ═══════════════════════════════════════════════════════════════
# TEST 3: rho < 1 — IS vs Stadler formula
#         Use the SAME tree, just change rho in both formulas
# ═══════════════════════════════════════════════════════════════
cat("─── Test 3: rho < 1 (IS vs Stadler, same tree) ─────────────\n")

# Use the tree from Test 1
cat(sprintf("  Using tree with %d tips\n", n_tips))

for (rho_val in c(1.0, 0.9, 0.7, 0.5, 0.3)) {
  # Stadler formula at this rho
  stadler_rho <- stadler_cr_loglik(brts, lam_true, mu_true, rho = rho_val)

  # IS estimate at this rho
  is_rho <- is_loglik(tree, c(lam_true, mu_true), rho = rho_val, S = 3000)

  # DDD with missnumspec (approximate: round to integer)
  missnumspec <- round(n_tips * (1 / rho_val - 1))
  ddd_rho <- tryCatch(
    DDD::bd_loglik(c(lam_true, mu_true, 0, 0), c(0, 0, 1, 0, 2),
                   brts, missnumspec),
    error = function(e) NA_real_
  )

  cat(sprintf("  rho=%.1f  missnumspec=%3d | Stadler=%8.2f  IS=%8.2f  DDD=%8.2f | IS-Stadler=%6.2f  IS-DDD=%6.2f\n",
              rho_val, missnumspec, stadler_rho, is_rho$fhat, ddd_rho,
              is_rho$fhat - stadler_rho, is_rho$fhat - ddd_rho))
}

cat("\n")


# ═══════════════════════════════════════════════════════════════
# TEST 4: Relative likelihood test
#         The RATIO of likelihoods at different rho should match
#         between IS and Stadler, even if absolute values differ
# ═══════════════════════════════════════════════════════════════
cat("─── Test 4: Relative likelihood ratios ──────────────────────\n")

rho_vals <- c(1.0, 0.9, 0.7, 0.5)
stadler_lls <- sapply(rho_vals, function(r) stadler_cr_loglik(brts, lam_true, mu_true, r))
is_lls <- sapply(rho_vals, function(r) is_loglik(tree, c(lam_true, mu_true), rho = r, S = 3000)$fhat)

cat("  Likelihood ratios relative to rho=1.0:\n")
cat(sprintf("  %-6s  %12s  %12s  %12s\n", "rho", "Stadler_ratio", "IS_ratio", "diff"))
for (i in seq_along(rho_vals)) {
  s_ratio <- stadler_lls[i] - stadler_lls[1]
  is_ratio <- is_lls[i] - is_lls[1]
  cat(sprintf("  %-6.1f  %12.4f  %12.4f  %12.4f\n",
              rho_vals[i], s_ratio, is_ratio, is_ratio - s_ratio))
}

cat("\n")


# ═══════════════════════════════════════════════════════════════
# TEST 5: MLE recovery with subsampled tree
#         Simulate full tree, subsample, estimate with rho correction
# ═══════════════════════════════════════════════════════════════
cat("─── Test 5: MLE recovery with subsampled tree ──────────────\n")

# Use a tree with ~25-50 tips
lam_sim <- 0.6
mu_sim  <- 0.15
big_tree <- NULL
for (s in 1:200) {
  set.seed(s)
  tr_sim <- simulate_tree(pars = c(lam_sim, mu_sim), max_t = 6,
                          model = "cr", max_tries = 100)
  if (!is.null(tr_sim$tes) && Ntip(tr_sim$tes) >= 25 && Ntip(tr_sim$tes) <= 50) {
    big_tree <- tr_sim$tes
    cat(sprintf("  Using seed=%d\n", s))
    break
  }
}
if (is.null(big_tree)) stop("Could not find suitable tree")
n_full <- Ntip(big_tree)
cat(sprintf("  Full tree: %d tips\n", n_full))

# Subsample: keep 70% of tips
rho_sub <- 0.7
set.seed(456)
drop_n <- round(n_full * (1 - rho_sub))
drop_tips <- sample(big_tree$tip.label, drop_n)
sub_tree <- ape::drop.tip(big_tree, drop_tips)
n_sub <- Ntip(sub_tree)
cat(sprintf("  Subsampled tree: %d tips (rho_actual = %.3f)\n",
            n_sub, n_sub / n_full))
rho_actual <- n_sub / n_full

# Profile likelihood: scan lambda at fixed mu for both approaches
cat("\n  Profile likelihood over lambda (mu fixed at true value):\n")
lam_grid <- seq(0.3, 1.0, by = 0.1)

cat(sprintf("  %-6s  %12s  %12s  %12s\n", "lambda", "full_rho1", "sub_rho1", "sub_rhoCorr"))
for (lam_test in lam_grid) {
  # Full tree, rho=1
  full_ll <- is_loglik(big_tree, c(lam_test, mu_sim), rho = 1.0, S = 1500)$fhat

  # Subsampled tree, rho=1 (WRONG - ignoring missing species)
  sub_ll_wrong <- is_loglik(sub_tree, c(lam_test, mu_sim), rho = 1.0, S = 1500)$fhat

  # Subsampled tree, rho=correct (RIGHT - accounting for missing species)
  sub_ll_right <- is_loglik(sub_tree, c(lam_test, mu_sim), rho = rho_actual, S = 1500)$fhat

  cat(sprintf("  %-6.2f  %12.2f  %12.2f  %12.2f\n",
              lam_test, full_ll, sub_ll_wrong, sub_ll_right))
}

# Find MLE for each
cat("\n  MLE comparison:\n")
mle_full <- estimate_rates(big_tree, method = "mcem", model = "cr",
                           control = list(lower_bound = c(0.1, 0.01),
                                          upper_bound = c(2.0, 1.0),
                                          sample_size = 300, max_iter = 30,
                                          num_threads = 4, rho = 1.0))

mle_sub_wrong <- estimate_rates(sub_tree, method = "mcem", model = "cr",
                                control = list(lower_bound = c(0.1, 0.01),
                                               upper_bound = c(2.0, 1.0),
                                               sample_size = 300, max_iter = 30,
                                               num_threads = 4, rho = 1.0))

mle_sub_right <- estimate_rates(sub_tree, method = "mcem", model = "cr",
                                control = list(lower_bound = c(0.1, 0.01),
                                               upper_bound = c(2.0, 1.0),
                                               sample_size = 300, max_iter = 30,
                                               num_threads = 4, rho = rho_actual))

cat(sprintf("  True parameters:                lambda=%.3f  mu=%.3f\n", lam_sim, mu_sim))
cat(sprintf("  Full tree (rho=1):              lambda=%.3f  mu=%.3f\n",
            mle_full$pars[1], mle_full$pars[2]))
cat(sprintf("  Subsampled tree (rho=1, WRONG): lambda=%.3f  mu=%.3f\n",
            mle_sub_wrong$pars[1], mle_sub_wrong$pars[2]))
cat(sprintf("  Subsampled tree (rho=%.2f, RIGHT): lambda=%.3f  mu=%.3f\n",
            rho_actual, mle_sub_right$pars[1], mle_sub_right$pars[2]))

cat("\n═══════════════════════════════════════════════════════════\n")
cat("  Done.\n")

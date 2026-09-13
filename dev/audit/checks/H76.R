## H76: batch survival_prob = mean(1/n_attempts), not the `done` fraction.
## Self-contained; ~30 s.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))

lam <- 0.5; mu <- 0.3; Tc <- 5; n <- 1000L
pars_mat <- cbind(beta_0 = rep(lam, n), gamma_0 = rep(mu, n))

## Closed-form single-attempt survival of a crown process (both crown
## lineages leave descendants at T; Kendall 1948).
p_ext <- mu * (1 - exp(-(lam - mu) * Tc)) / (lam - mu * exp(-(lam - mu) * Tc))
p     <- (1 - p_ext)^2
E_inv_n <- function(p, k) sum(p * (1 - p)^(0:k) / (1:(k + 1)))   # E[1/n_attempts * 1{done}]
E_done  <- function(p, k) 1 - (1 - p)^(k + 1)

cat(sprintf("Closed-form per-attempt crown survival p = %.4f\n\n", p))
cat(sprintf("%-9s %-14s %-14s %-14s %-14s %-16s\n",
            "max_tries", "reported", "E[mean(1/n)]", "done_frac", "E[done_frac]", "values seen"))
res <- list()
for (k in c(0L, 1L, 5L)) {
  b   <- simulate_tree(pars = pars_mat, max_t = Tc, model = "cr",
                       max_tries = k, useDDD = FALSE, num_threads = 1L)
  sp  <- vapply(b$simulations, `[[`, 0.0, "survival_prob")
  done <- mean(vapply(b$simulations, function(s) s$status == "done", TRUE))
  res[[as.character(k)]] <- c(reported = b$survival_prob, done = done)
  cat(sprintf("%-9d %-14.4f %-14.4f %-14.4f %-14.4f %-16s\n",
              k, b$survival_prob, E_inv_n(p, k), done, E_done(p, k),
              paste(sort(unique(round(sp, 3))), collapse = ",")))
  ## what train_GAM() would see for this batch (it reads status, not survival_prob)
  surv_gam_view <- mean(sapply(b$simulations, function(s) as.integer(s$status == "done")))
  stopifnot(isTRUE(all.equal(surv_gam_view, done)))
  if (k == 0L) stopifnot(isTRUE(all.equal(b$survival_prob, done)))  # identical at max_tries = 0
}

## Decision: at max_tries >= 1 the reported batch value must differ from the
## done fraction by more than Monte Carlo noise (binomial SE ~ 0.016 at n = 1000).
gap1 <- res[["1"]]["done"] - res[["1"]]["reported"]
gap5 <- res[["5"]]["done"] - res[["5"]]["reported"]
cat(sprintf("\ngap done - reported: max_tries=1: %.4f (expected %.4f); max_tries=5: %.4f (expected %.4f)\n",
            gap1, E_done(p, 1) - E_inv_n(p, 1), gap5, E_done(p, 5) - E_inv_n(p, 5)))
cat("Consumers of survival_prob in R/ outside simulate.R:\n")
print(system("grep -rn survival_prob /Users/pancho/Code/emphasis/R | grep -v 'R/simulate.R'", intern = TRUE))
cat("(train_GAM/.test_feasibility use status == 'done' with max_tries = 0)\n")

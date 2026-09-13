.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
# Analytic P(a lineage alive at T) under CR (Kendall): P_ext = mu(1-e^{-rT})/(lam - mu e^{-rT})
p_ext <- function(lam, mu, T) { r <- lam - mu; mu * (1 - exp(-r * T)) / (lam - mu * exp(-r * T)) }
p_both <- function(lam, mu, T) (1 - p_ext(lam, mu, T))^2   # both crown lineages alive at T
# DDD's cond = 1 normaliser for CR (soc = 2): log-lik minus log of P(both crown lineages survive)
ddd_cond_term <- function(lam, mu, brts) {
  l0 <- DDD::bd_loglik(pars1 = c(lam, mu, 0, 0), pars2 = c(0, 0, 1, 0, 2),
                       brts = brts, missnumspec = 0)
  l1 <- DDD::bd_loglik(pars1 = c(lam, mu, 0, 0), pars2 = c(0, 1, 1, 0, 2),
                       brts = brts, missnumspec = 0)
  exp(l0 - l1)   # = normaliser P(event)
}
T <- 5; n_sim <- 400
grid <- rbind(c(0.5, 0.1), c(0.5, 0.3), c(0.5, 0.45), c(1.0, 0.8))
cat(sprintf("%-6s %-6s %-9s %-9s %-9s %-10s\n", "lam", "mu", "sim_done", "analytic", "DDD_norm", "sim_tooLarge"))
for (i in seq_len(nrow(grid))) {
  lam <- grid[i, 1]; mu <- grid[i, 2]
  pm <- matrix(rep(c(lam, mu), each = n_sim), ncol = 2)
  s <- simulate_tree(pars = pm, max_t = T, model = "cr", max_tries = 0,
                     useDDD = FALSE, max_lin = 500)
  st <- vapply(s$simulations, function(x) x$status, "")
  # a synthetic brts just to obtain DDD's normaliser (it depends only on T)
  brts <- c(T, 3, 1)
  cat(sprintf("%-6.2f %-6.2f %-9.3f %-9.3f %-9.3f %-10.3f\n", lam, mu,
              mean(st == "done"), p_both(lam, mu, T), ddd_cond_term(lam, mu, brts),
              mean(st == "too_large")))
}
# And the object actually used as the conditioning term: the survival GAM
set.seed(3)
pars_mat <- cbind(beta_0 = runif(300, 0.2, 1.2), gamma_0 = runif(300, 0.0, 0.9))
pars_mat <- pars_mat[pars_mat[, 2] < pars_mat[, 1], ]
sims <- simulate_tree(pars = pars_mat, max_t = T, model = "cr", max_tries = 0,
                      useDDD = FALSE, max_lin = 500)
g <- suppressMessages(train_GAM(sims$simulations, pars_mat, model = "cr"))
nd <- data.frame(beta_0 = grid[, 1], gamma_0 = grid[, 2])
cat("\nGAM P_surv vs analytic (1-P_ext)^2 at the grid:\n")
print(cbind(nd, gam = round(predict_survival(g, nd), 3),
            analytic = round(p_both(grid[, 1], grid[, 2], T), 3)))

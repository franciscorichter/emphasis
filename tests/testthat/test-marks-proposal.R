# The mark-sampling proposal (R/marks-proposal.R).  With one mark it samples the exact posterior of
# cr: every importance weight equals the reconstructed-tree likelihood, which the p1 form gives.

p1 <- function(x, l, m) { d <- l - m; d^2 * exp(-d * x) / (l - m * exp(-d * x))^2 }
one <- list(S = 1L, u = 0, log_clam = 0, log_cmu = 0, C = 1, P = 1, pi = 1)

test_that("the closed-form branch integrals match quadrature", {
  for (lam in c(0.35, 3, 120)) {
    mu <- lam / 5; T <- 14; a <- 2; b <- 9
    expect_equal(.pm_int_l1q(a, b, lam, mu, T),
                 stats::integrate(function(t) lam * (1 - .pm_q(t, lam, mu, T)), a, b)$value, tolerance = 1e-6)
    expect_equal(.pm_int_muq(a, b, lam, mu, T),
                 stats::integrate(function(t) mu / .pm_q(t, lam, mu, T), a, b)$value, tolerance = 1e-6)
  }
})

test_that("with one mark every weight is the exact cr likelihood", {
  set.seed(2)
  repeat {
    r <- .simulate_marks(c(log(0.4), log(0.1)), "cr", 6, one, max_N = 2000, max_tries = 1L)
    if (r$status == "done") { nt <- sum(r$Ltable[, 4] == -1); if (nt >= 8 && nt <= 30) break }
  }
  phy <- DDD::L2phylo(r$Ltable, dropextinct = TRUE)
  x <- as.numeric(sort(ape::branching.times(phy), decreasing = TRUE))
  mb <- .resolve_model("cr")
  for (th in list(c(0.4, 0.1), c(0.6, 0.3))) {
    ref <- 2 * log(p1(x[1], th[1], th[2])) + sum(log(th[1] * p1(x[-1], th[1], th[2])))
    set.seed(3)
    d <- .pm_draws(phy, th[1], th[2], one, 8L)
    lw <- .eval_logf_marks(.expand_pars(log(th), mb), d$trees, mb, one) - d$logg
    expect_equal(d$log_pobs, ref, tolerance = 1e-9)
    expect_equal(lw, rep(ref, 8L), tolerance = 1e-9)
  }
})

test_that("the proposal needs cr in the marks MCEM", {
  set.seed(1)
  phy <- ape::rcoal(6)
  expect_error(.mcem_marks(phy, "dd", c(log(0.3), 0, log(0.1)), log(0.3),
                           function(x) one, n_iter = 1L, N = 2L, aware = "sample"), "needs model cr")
})

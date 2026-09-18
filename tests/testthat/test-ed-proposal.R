# The ED-aware proposal of the thinning sampler (Model::ed_proposal_seg_t,
# inst/include/model.hpp): for an ED model the missing births are drawn at
# the model's total speciation rate over the alive lineages, each at its own
# ED, and a birth is attached to a lineage in proportion to that lineage's
# rate.  The mean-field proposal (nh_rate: N * lambda at ED = 0, uniform
# attachment) stays reachable through model slot 4 == 2, so the two can be
# run on the same trees.  Four checks:
#   1. at beta_ED = gamma_ED = 0 the two densities agree on trees drawn by
#      either proposal, at rho = 1 and rho < 1, on both links;
#   2. the envelope dominates the intensity: no thinning violations;
#   3. the two proposals estimate the same marginal likelihood of an
#      observed tree (importance sampling is unbiased under either), and where
#      the effect is strong the ED-aware one does it with the larger ESS;
#   4. estimate_rates runs under both, through control$proposal.

mb_ed <- emphasis:::.resolve_model("ned")
mb_mf <- { m <- mb_ed; m[4L] <- 2L; m }

aug_with <- function(brts, pars, mb, n, seed, link = 0L, rho = 1.0, max_missing = 5000L) {
  emphasis:::augment_trees(brts, emphasis:::.expand_pars(pars, mb_ed), sample_size = n,
                           maxN = 50L * n, max_missing = max_missing, max_lambda = 1e6,
                           num_threads = 1L, model = mb, link = link, rho = rho,
                           parent_tip_start = emphasis:::.pts(brts), seed = seed,
                           parent_id = emphasis:::.pid(brts))
}

logfg <- function(pars, trees, mb, link = 0L, rho = 1.0) {
  emphasis:::eval_logf(emphasis:::.expand_pars(pars, mb_ed), trees, model = mb, link = link, rho = rho)
}

# log of the mean of exp(lw), its standard error on the log scale, and the ESS
is_summary <- function(lw) {
  lw  <- lw[is.finite(lw) | lw == -Inf]
  m   <- max(lw)
  w   <- exp(lw - m)
  est <- m + log(mean(w))
  se  <- stats::sd(w) / (sqrt(length(w)) * mean(w))
  list(est = est, se = se, ess = sum(w)^2 / sum(w^2), n = length(w))
}

test_that("at beta_ED = gamma_ED = 0 the ED-aware density is the mean-field density", {
  set.seed(21)
  phy  <- ape::rcoal(12)
  brts <- emphasis:::.extract_brts(phy)
  ned0 <- c(0.7, -0.01, 0.0, 0.2, 0.0, 0.0)   # beta_0, beta_N, beta_ED, gamma_0, gamma_N, gamma_ED
  for (link in c(0L, 1L)) for (rho in c(1.0, 0.7)) {
    p <- if (link == 1L) c(log(0.7), -0.01, 0.0, log(0.2), 0.0, 0.0) else ned0
    for (mb_draw in list(mb_ed, mb_mf)) {
      aug <- aug_with(brts, p, mb_draw, n = 15L, seed = 3L + link, link = link, rho = rho)
      expect_gt(length(aug$trees), 0L)
      g_ed <- logfg(p, aug$trees, mb_ed, link = link, rho = rho)
      g_mf <- logfg(p, aug$trees, mb_mf, link = link, rho = rho)
      expect_true(all(is.finite(g_ed$logg)))
      expect_equal(g_ed$logg, g_mf$logg, tolerance = 1e-9)
      expect_equal(g_ed$logf, g_mf$logf, tolerance = 1e-12)
    }
  }
})

test_that("the ED-aware envelope dominates the intensity, for either sign of beta_ED and both links", {
  set.seed(22)
  phy  <- ape::rcoal(15)
  brts <- emphasis:::.extract_brts(phy)
  emphasis:::thinning_envelope_violations(TRUE)
  for (bED in c(-0.15, 0.1)) for (link in c(0L, 1L)) {
    p <- if (link == 1L) c(log(0.8), -0.01, bED, log(0.3), 0.0, 0.5 * bED)
         else c(0.8, -0.01, bED, 0.3, 0.0, 0.5 * bED)
    aug <- aug_with(brts, p, mb_ed, n = 30L, seed = 11L, link = link, rho = 0.8)
    expect_gt(length(aug$trees), 0L)
    g <- logfg(p, aug$trees, mb_ed, link = link, rho = 0.8)
    expect_true(all(is.finite(g$logg)))
  }
  expect_equal(emphasis:::thinning_envelope_violations(TRUE), 0)
})

test_that("both proposals estimate the same marginal likelihood; the ED-aware one with the larger ESS", {
  set.seed(23)
  truth <- c(0.6, -0.004, -0.06, 0.2, 0.0, 0.0)
  s <- simulate_tree(pars = truth, max_t = 8, model = "ned", max_lin = 3000L, num_threads = 1L)
  expect_equal(s$status, "done")
  brts <- emphasis:::.extract_brts(s$tes)
  n <- 400L
  a_ed <- aug_with(brts, truth, mb_ed, n = n, seed = 5L)
  a_mf <- aug_with(brts, truth, mb_mf, n = n, seed = 6L)
  g_ed <- logfg(truth, a_ed$trees, mb_ed)
  g_mf <- logfg(truth, a_mf$trees, mb_mf)
  s_ed <- is_summary(g_ed$logf - g_ed$logg)
  s_mf <- is_summary(g_mf$logf - g_mf$logg)
  expect_true(is.finite(s_ed$est) && is.finite(s_mf$est))
  # the same target: the two estimates agree within their Monte Carlo error
  expect_lt(abs(s_ed$est - s_mf$est), 4 * sqrt(s_ed$se^2 + s_mf$se^2))
  # and the ED-aware proposal is the sharper one
  expect_gt(s_ed$ess, s_mf$ess)
  # a stronger effect with turnover: the mean-field weights collapse, the
  # ED-aware ones do not
  truth2 <- c(0.9, -0.006, -0.2, 0.45, 0.0, 0.0)
  s2 <- simulate_tree(pars = truth2, max_t = 8, model = "ned", max_lin = 3000L, num_threads = 1L)
  if (identical(s2$status, "done") && length(s2$tes$tip.label) >= 8L) {
    b2 <- emphasis:::.extract_brts(s2$tes)
    e2 <- aug_with(b2, truth2, mb_ed, n = 200L, seed = 7L, max_missing = 20000L)
    m2 <- aug_with(b2, truth2, mb_mf, n = 200L, seed = 8L, max_missing = 20000L)
    if (length(e2$trees) == 200L && length(m2$trees) == 200L) {
      w_e2 <- is_summary(with(logfg(truth2, e2$trees, mb_ed), logf - logg))
      w_m2 <- is_summary(with(logfg(truth2, m2$trees, mb_mf), logf - logg))
      expect_gt(w_e2$ess, 2 * w_m2$ess)
    }
  }
})

test_that("estimate_rates runs an ED model under either proposal", {
  set.seed(24)
  phy <- ape::rcoal(12)
  for (proposal in c("ed", "meanfield")) {
    fit <- estimate_rates(phy, method = "mcem", model = "ned",
                          init_pars = c(0.8, -0.01, 0.0, 0.2, 0.0, 0.0),
                          control = list(sampling = "dynamic_fresh", proposal = proposal,
                                         sample_size = 20L, max_iter = 2L, num_threads = 1L,
                                         lower_bound = c(0.01, -0.2, -1, 0.0, -0.2, -1),
                                         upper_bound = c(3, 0.2, 1, 2, 0.2, 1),
                                         verbose = FALSE))
    expect_s3_class(fit, "emphasis_fit")
    expect_length(fit$pars, 6L)
    expect_true(is.finite(fit$loglik))
  }
  expect_error(estimate_rates(phy, method = "mcem", model = "ned",
                              init_pars = c(0.8, -0.01, 0.0, 0.2, 0.0, 0.0),
                              control = list(sampling = "dynamic_fresh", proposal = "uniform",
                                             sample_size = 5L, max_iter = 1L, num_threads = 1L,
                                             lower_bound = c(0.01, -0.2, -1, 0.0, -0.2, -1),
                                             upper_bound = c(3, 0.2, 1, 2, 0.2, 1))),
               "proposal")
})

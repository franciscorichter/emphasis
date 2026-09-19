# A log of a Monte Carlo mean is not the mean of a log.  The EM objective is
# safe -- it averages log f directly (src/M_step.cpp) -- but the
# log-likelihood a fit reports is log(mean(w)), which Jensen puts below the
# truth, by more the heavier the weights.  Unequal precision between two fits
# therefore biases the AIC comparison, always against the harder-to-sample
# model.

test_that("the gap is zero for equal weights and grows as ESS falls", {
  expect_equal(emphasis:::.jensen_gap(200, 200)$gap, 0)
  expect_gt(emphasis:::.jensen_gap(20, 200)$gap,
            emphasis:::.jensen_gap(100, 200)$gap)
  expect_equal(emphasis:::.jensen_gap(1, 200)$gap, 0.5 * (1 - 1 / 200))
  # a degenerate or missing draw yields no claim rather than a wrong one
  expect_true(is.na(emphasis:::.jensen_gap(NA, 200)$gap))
  expect_true(is.na(emphasis:::.jensen_gap(10, 0)$gap))
  expect_false(emphasis:::.jensen_gap(150, 200)$heavy)
  expect_true(emphasis:::.jensen_gap(9, 200)$heavy)
})

test_that("the gap tracks the real bias, and `heavy` marks where it stops", {
  skip_on_cran()
  set.seed(4)
  n <- 200L; R <- 6000L
  for (sdw in c(0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)) {
    lw <- matrix(stats::rnorm(R * n, -sdw^2 / 2, sdw), R, n)   # E[w] = 1
    w  <- exp(lw)
    truth <- -mean(log(rowMeans(w)))                           # > 0, downward
    ess   <- stats::median(rowSums(w)^2 / rowSums(w^2))
    claim <- emphasis:::.jensen_gap(ess, n)

    # where it overstates, the bias is negligible in absolute terms anyway
    if (claim$gap > truth) expect_lt(claim$gap, 2e-3)
    # where it understates by more than half, the draw is flagged heavy
    if (claim$gap < truth / 2) expect_true(claim$heavy)
    # and a flagged draw is one where the bias is no longer ignorable
    if (claim$heavy) expect_gt(truth, 0.01)
  }
})

test_that(".is_summary carries the gap alongside the estimate", {
  s <- emphasis:::.is_summary(rep(0, 50))            # uniform weights
  expect_equal(s$ess, 50)
  expect_equal(s$gap, 0)
  expect_false(s$gap_heavy)
  s2 <- emphasis:::.is_summary(c(0, rep(-40, 49)))   # one weight dominates
  expect_lt(s2$ess, 1.01)
  expect_gt(s2$gap, 0.48)
  expect_true(s2$gap_heavy)
})

mk_fit <- function(loglik, n_pars, ess, n = 200) {
  structure(list(
    pars = c(beta_0 = 0.5), loglik = loglik, loglik_var = NA_real_,
    n_pars = as.integer(n_pars), AIC = 2 * n_pars - 2 * loglik,
    method = "mcem", model = c(0L, 0L, 0L), cond = FALSE,
    details = list(final_IS = list(
      lw = rep(0, n), ESS = ess,
      gap = emphasis:::.jensen_gap(ess, n)$gap,
      gap_heavy = emphasis:::.jensen_gap(ess, n)$heavy))
  ), class = "emphasis_fit")
}

test_that("compare_models is quiet when both fits are sampled equally well", {
  expect_silent(tb <- compare_models(a = mk_fit(-100, 2, 190),
                                     b = mk_fit(-102, 2, 190)))
  expect_true(all(c("ESS", "loglik_gap") %in% names(tb)))
  expect_length(attr(tb, "ranking_at_risk"), 0L)
})

test_that("compare_models warns when unequal precision could flip the ranking", {
  # the loser is the badly sampled one, and its handicap exceeds the AIC gap
  expect_warning(tb <- compare_models(a = mk_fit(-100, 2, 190),
                                      b = mk_fit(-100.2, 2, 2)),
                 "not safe")
  expect_identical(attr(tb, "ranking_at_risk"), "b")
  expect_gt(tb$AIC_swing[tb$model == "b"], tb$delta_AIC[tb$model == "b"])
})

test_that("a heavy draw is unquantified, not safe, even behind a wide AIC gap", {
  # swing alone would look small against a 40-point gap, but the loser's draw
  # is heavy, so its gap is an underestimate of unknown size
  expect_warning(tb <- compare_models(a = mk_fit(-100, 2, 195),
                                      b = mk_fit(-120, 2, 5)),
                 "not safe")
  expect_identical(attr(tb, "ranking_at_risk"), "b")
})

test_that("a fit that predates the fields still compares safely", {
  old <- structure(list(
    pars = c(beta_0 = 0.5), loglik = -100, loglik_var = NA_real_,
    n_pars = 2L, AIC = 204, method = "mcem", model = c(0L, 0L, 0L),
    cond = FALSE,
    details = list(final_IS = list(lw = rep(0, 100), ESS = 4))  # no gap fields
  ), class = "emphasis_fit")
  expect_equal(emphasis:::.fit_ess(old), 4)
  expect_equal(emphasis:::.fit_gap(old), 0.5 * (1 / 4 - 1 / 100))
  expect_true(emphasis:::.fit_gap_heavy(old))
})

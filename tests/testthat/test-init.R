brts_fixed <- c(10, 8.7, 7.9, 7.1, 6.6, 5.8, 5.2, 4.9, 4.1, 3.8, 3.3, 2.9, 2.6,
                2.2, 1.9, 1.7, 1.4, 1.2, 0.9, 0.7, 0.5, 0.35, 0.2, 0.1)

test_that(".cr_mle agrees with DDD::bd_ML, unconditioned and conditioned", {
  for (cnd in 0:1) {
    mine <- .cr_mle(brts_fixed, cond = cnd == 1)
    ref <- suppressMessages(DDD::bd_ML(brts = brts_fixed, cond = cnd, btorph = 0, soc = 2,
                                       verbose = FALSE, tdmodel = 0, idparsopt = 1:2,
                                       initparsopt = c(0.4, 0.1)))
    expect_equal(mine[1], ref$lambda0, tolerance = 1e-3)
    expect_equal(mine[2], ref$mu0, tolerance = 1e-3)
  }
})

test_that(".cr_mle's likelihood is the cr log-likelihood the BDI E-step reports", {
  for (p in list(c(0.5, 0.2), c(0.4, 0.05))) {
    e <- .augment_tree_bdi(tree = brts_fixed, pars = .expand_pars(log(p), c(0L, 0L, 0L)),
                           model_bin = c(0L, 0L, 0L), sample_size = 5L, max_missing = 10000L,
                           link = 1L, rho = 1)
    l <- p[1]; m <- p[2]; r <- l - m; x <- brts_fixed
    lp1 <- 2 * log(r) - r * x - 2 * log(l - m * exp(-r * x))
    expect_equal(e$fhat, 2 * lp1[1] + sum(log(l) + lp1[-1]), tolerance = 1e-8)
  }
})

test_that(".default_init places the cr estimate in the intercepts, slopes at zero", {
  est <- .cr_mle(brts_fixed)
  lb <- c(-5, -0.2, -8, -0.2); ub <- c(2, 0.2, 2, 0.2)
  ini <- .default_init(brts_fixed, c(1L, 0L, 0L), 1L, lb, ub)
  expect_equal(ini[c(2, 4)], c(0, 0))
  expect_equal(ini[1], log(est[1]))
  expect_equal(ini[3], max(log(est[2]), lb[3] + 1e-3 * (ub[3] - lb[3])))
})

test_that(".em_remaining projects a slow EM and passes Monte Carlo noise", {
  th <- 2 + 1.5 * 0.97^(0:30)              # geometric approach to 2, rate 0.97
  w <- utils::tail(th, 8)
  expect_equal(.em_remaining(cbind(w), 0.01), abs(w[8] - 2) / w[8], tolerance = 1e-8)
  expect_gt(.em_remaining(cbind(w), 0.01), 1e-2)  # the step alone is below 1e-2
  expect_lt(abs(diff(utils::tail(w, 2))) / w[7], 1e-2)
  set.seed(1)
  noise <- 2 + stats::rnorm(8, sd = 0.002)
  expect_lt(.em_remaining(cbind(noise), 0.01), 1e-2)
  expect_identical(.em_remaining(cbind(c(1, 1, 1)), 0.01), 0)
  expect_identical(.em_remaining(cbind(c(1, 1.1, 1.3, 1.6)), 0.01), Inf)
})

# Tests for auto_bounds() and its helpers: the scale the intercepts live on
# under each link, the feasibility criterion, and the observed covariate values.

# The branching times of a tree simulated from cr (lambda = 0.40, mu = 0.06,
# crown age 10), kept verbatim so the exact MLE below is the same number on
# every run.
ab_test_brts <- c(
  10.0000000000, 9.9509298541, 6.4895174503, 5.9586496353, 5.7669186592,
  5.2027902603, 5.1180691719, 5.1167187691, 4.5694217682, 4.5278587341,
  4.3009586334, 4.1878695488, 3.8652491570, 3.7823085785, 3.2337293625,
  3.2106838226, 2.8642201424, 2.1649961472, 2.0914907455, 2.0042562485,
  1.9750614166, 1.7716503143, 1.7539339066, 1.7431106567, 1.4805088043,
  1.4731340408, 1.2839851379, 1.2531337738, 1.2497835159, 1.0027036667,
  0.8967571259, 0.8942689896, 0.8594207764, 0.7203598022, 0.6892919540,
  0.6423854828, 0.5741643906, 0.3768663406, 0.2965583801, 0.2340002060,
  0.1837100983, 0.1440057755, 0.0903158188, 0.0758924484)

# Intercept that puts the per-lineage rate at `rate` when the covariate term of
# the linear predictor is zero, written out independently of .rate_intercept().
ab_intercept <- function(rate, link) {
  switch(link,
         linear      = rate,
         exponential = log(rate),
         gaussian    = rate * exp(0.5))
}


# --- the scale each link's intercept lives on ------------------------------

test_that(".rate_intercept inverts each link at a zero covariate term", {
  rate <- 0.37
  # linear: rate = max(0, beta_0)
  expect_equal(emphasis:::.rate_intercept(rate, 0L), rate)
  # exponential: rate = exp(beta_0)
  expect_equal(exp(emphasis:::.rate_intercept(rate, 1L)), rate)
  # gaussian: rate = beta_0 * exp(-(eta_cov - 1)^2 / 2), eta_cov = 0
  expect_equal(emphasis:::.rate_intercept(rate, 2L) * exp(-0.5), rate)
})

test_that(".wide_bounds puts the gaussian intercepts on the natural scale", {
  model_bin <- c(0L, 0L, 0L)
  lin <- emphasis:::.wide_bounds(model_bin, 0L, max_t = 20, n_tips = 23)
  gau <- emphasis:::.wide_bounds(model_bin, 2L, max_t = 20, n_tips = 23)

  # beta_0 above zero, gamma_0 at zero
  expect_gt(gau$lb[1], 0)
  expect_equal(gau$lb[2], 0)

  # Under cr the covariate term is zero on both links and both entries are
  # intercepts, so the gaussian box is the linear box scaled by exp(1/2).
  expect_equal(gau$lb, lin$lb * exp(0.5))
  expect_equal(gau$ub, lin$ub * exp(0.5))
})

test_that(".wide_bounds leaves the linear and exponential boxes alone", {
  # The linear box is the rate box itself; the exponential box is its log,
  # with the extra 3 log-units of headroom below on gamma_0.
  max_t <- 20; n_tips <- 23
  r_hat  <- log(n_tips / 2) / max_t
  lam_hi <- max(log(50 * n_tips / 2) / max_t + 0.5 * r_hat, 0.1)
  lam_lo <- max(0.1 * r_hat, 1e-4)

  lin <- emphasis:::.wide_bounds(c(0L, 0L, 0L), 0L, max_t, n_tips)
  expect_equal(unname(lin$lb), c(lam_lo, 0))
  expect_equal(unname(lin$ub), c(lam_hi, lam_hi * 0.9))

  exps <- emphasis:::.wide_bounds(c(0L, 0L, 0L), 1L, max_t, n_tips)
  expect_equal(unname(exps$lb), c(log(lam_lo), log(lam_lo) - 3))
  expect_equal(unname(exps$ub), c(log(lam_hi), log(lam_hi)))
})

test_that("the gaussian covariate bounds are untouched by the rescaling", {
  # Only the intercepts are rates; the slopes are slopes on both links.
  lin <- emphasis:::.wide_bounds(c(1L, 0L, 0L), 0L, max_t = 20, n_tips = 23)
  gau <- emphasis:::.wide_bounds(c(1L, 0L, 0L), 2L, max_t = 20, n_tips = 23)
  expect_equal(gau$lb[c(2, 4)], lin$lb[c(2, 4)])   # beta_N, gamma_N
  expect_equal(gau$ub[c(2, 4)], lin$ub[c(2, 4)])
})


# --- the feasibility hole --------------------------------------------------

test_that(".crown_rate_positive rejects a non-positive total rate", {
  cr <- c(0L, 0L, 0L)
  # gaussian: beta_0 is the peak rate, so beta_0 <= 0 is a non-positive rate
  expect_false(emphasis:::.crown_rate_positive(c(-2.257, -4.049), cr, 2L))
  expect_false(emphasis:::.crown_rate_positive(c(0, 0), cr, 2L))
  expect_true(emphasis:::.crown_rate_positive(c(0.1726, 0.0288), cr, 2L))

  # linear: max(0, eta) on both rates
  expect_false(emphasis:::.crown_rate_positive(c(-0.1, -0.02), cr, 0L))
  expect_true(emphasis:::.crown_rate_positive(c(0.1, 0), cr, 0L))

  # exponential: exp() is positive everywhere
  expect_true(emphasis:::.crown_rate_positive(c(-20, -30), cr, 1L))
})

test_that(".crown_rate_positive reads the N slope at the crown", {
  dd <- c(1L, 0L, 0L)
  # At the crown N = 2, so the linear predictor is beta_0 + 2 * beta_N.
  expect_false(emphasis:::.crown_rate_positive(c(0.4, -0.3, 0.1, -0.1), dd, 0L))
  expect_true(emphasis:::.crown_rate_positive(c(0.4, -0.1, 0.1, -0.1), dd, 0L))
})

test_that(".test_feasibility rejects a gaussian negative-rate point", {
  skip_on_cran()
  # This is the point auto_bounds used to pick as its centre on bird.orders:
  # the log-scale intercept read as a natural-scale peak rate. Every draw ends
  # at once with the 2-row crown L-table and status "done".
  expect_false(emphasis:::.test_feasibility(
    c(-2.257, -4.049), "cr", "gaussian", max_t = 20, max_lin = 500L,
    n_test = 5L, tip_lo = 2, tip_hi = 230, num_threads = 1L))
})

test_that(".test_feasibility rejects a rate that yields only the crown", {
  skip_on_cran()
  # Positive but negligible rates: the simulator does reach the present, and
  # does it with 2 tips, so there is no branching event past the crown.
  expect_false(emphasis:::.test_feasibility(
    c(1e-5, 1e-6), "cr", "linear", max_t = 20, max_lin = 500L,
    n_test = 5L, tip_lo = 2, tip_hi = 230, num_threads = 1L))
})


# --- auto_bounds end to end ------------------------------------------------

test_that("auto_bounds returns a positive gaussian box", {
  skip_on_cran()
  skip_if_not_installed("ape")
  data(bird.orders, package = "ape", envir = environment())

  for (md in c("cr", "dd")) {
    ab <- auto_bounds(bird.orders, model = md, link = "gaussian",
                      n_test = 3L, bisect_steps = 4L, train_surv_gam = FALSE,
                      verbose = FALSE, num_threads = 1L)
    expect_gt(ab$lower_bound["beta_0"], 0)
    expect_gt(ab$upper_bound["beta_0"], ab$lower_bound["beta_0"])
    expect_gte(ab$lower_bound["gamma_0"], 0)
    # The centre is the observed net rate carried onto the peak scale.
    max_t   <- max(ape::branching.times(bird.orders))
    r_hat   <- log(ape::Ntip(bird.orders) / 2) / max_t
    lam_hat <- r_hat + 0.2 * r_hat
    expect_equal(unname(ab$center[1]) * exp(-0.5), lam_hat)
  }
})

test_that("the exact MLE lies inside the box on all three links", {
  skip_on_cran()
  skip_if_not_installed("DDD")

  # DDD::bd_ML is the exact constant-rate MLE conditioned on crown survival
  # and crown age; it is deterministic given the branching times.
  invisible(utils::capture.output(
    mle <- DDD::bd_ML(brts = ab_test_brts, cond = 1, btorph = 1, soc = 2,
                      initparsopt = c(0.4, 0.06), idparsopt = 1:2,
                      tdmodel = 0, verbose = FALSE)))
  expect_gt(mle$mu0, 0)   # a log-link box cannot contain mu = 0

  for (lk in c("linear", "exponential", "gaussian")) {
    ab <- auto_bounds(ab_test_brts, model = "cr", link = lk,
                      n_test = 5L, bisect_steps = 6L, train_surv_gam = FALSE,
                      verbose = FALSE, num_threads = 1L)
    b0 <- ab_intercept(mle$lambda0, lk)
    g0 <- ab_intercept(mle$mu0, lk)
    expect_true(b0 >= ab$lower_bound[1] && b0 <= ab$upper_bound[1],
                info = sprintf("%s: beta_0 %.4f outside [%.4f, %.4f]",
                               lk, b0, ab$lower_bound[1], ab$upper_bound[1]))
    expect_true(g0 >= ab$lower_bound[2] && g0 <= ab$upper_bound[2],
                info = sprintf("%s: gamma_0 %.4f outside [%.4f, %.4f]",
                               lk, g0, ab$lower_bound[2], ab$upper_bound[2]))
  }
})


# --- the M covariate at the observed tree ----------------------------------

test_that(".mean_pendant_age is the mean terminal edge length", {
  skip_if_not_installed("ape")
  data(bird.orders, package = "ape", envir = environment())
  brts <- emphasis:::.extract_brts(bird.orders)
  tip_edge <- bird.orders$edge[, 2] <= ape::Ntip(bird.orders)
  expect_equal(emphasis:::.mean_pendant_age(brts),
               mean(bird.orders$edge.length[tip_edge]))
})

test_that(".mean_pendant_age is the crown age for a two-tip tree", {
  phy <- ape::read.tree(text = "(a:3,b:3);")
  brts <- emphasis:::.extract_brts(phy)
  expect_equal(emphasis:::.mean_pendant_age(brts), 3)
})

test_that(".mean_pendant_age falls back without topology", {
  # A bare branching-time vector carries no parent_tip_start attribute, so the
  # tips cannot be told from the internal nodes.
  brts <- c(10, 7, 4, 1)
  expect_equal(emphasis:::.pts(brts), numeric(0))
  expect_equal(emphasis:::.mean_pendant_age(brts), 5)
})

test_that(".observed_covariates reports the mean pendant age for M", {
  skip_if_not_installed("ape")
  data(bird.orders, package = "ape", envir = environment())
  brts <- emphasis:::.extract_brts(bird.orders)
  tip_edge <- bird.orders$edge[, 2] <= ape::Ntip(bird.orders)
  mean_pendant <- mean(bird.orders$edge.length[tip_edge])

  oc <- emphasis:::.observed_covariates(brts, c(0L, 1L, 0L))
  expect_equal(oc[[1]]$name, "M")
  expect_equal(oc[[1]]$X_obs, mean_pendant)
  # The branching-time mean is a different quantity, and was what this used.
  expect_gt(mean(brts), mean_pendant)
})

test_that(".observed_covariates reports N as the tip count", {
  skip_if_not_installed("ape")
  data(bird.orders, package = "ape", envir = environment())
  brts <- emphasis:::.extract_brts(bird.orders)
  oc <- emphasis:::.observed_covariates(brts, c(1L, 0L, 0L))
  expect_equal(oc[[1]]$name, "N")
  expect_equal(oc[[1]]$X_obs, ape::Ntip(bird.orders))
})

test_that("the feasibility test judges survivors, not survival (H74)", {
  skip_on_cran()
  # The observed tree survived to the present, so feasibility is a question
  # about the surviving clades a parameter vector produces, not about how often
  # it produces one.  Judging unconditionally made the test a survival test and
  # excluded high-turnover regions from the box before any likelihood was
  # evaluated.  Measured over turnover 0 to 0.9: containment of the exact MLE
  # rose from 0.75 to 0.96 overall and from 0.33 to 1.00 at turnover 0.9, with
  # the mu axis widening (0.55 to 1.65 mean) and the lambda axis narrowing.
  #
  # Containment is improved, not guaranteed: at high turnover the MLE is itself
  # unstable and can sit well outside any box built from the observed tree's
  # own rate scale.  The mechanism is what is pinned here.
  set.seed(4)
  phy  <- TreeSim::sim.bd.taxa(n = 30, numbsim = 1, lambda = 1, mu = 0.9,
                               complete = FALSE)[[1]]
  brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
  set.seed(99)
  ab <- auto_bounds(brts, model = "cr", link = "linear", num_threads = 1L,
                    verbose = FALSE)
  mu_w <- unname(ab$upper_bound[2] - ab$lower_bound[2])
  lam_w <- unname(ab$upper_bound[1] - ab$lower_bound[1])
  # the mu axis must be able to hold a high-turnover optimum at all: before the
  # change it was a fraction of the lambda axis on trees like this one
  expect_gt(mu_w, 0.5)
  expect_gt(mu_w / lam_w, 0.5)
  # and the test itself must accept a vector whose survivors match the data
  # even when most of its draws die out
  expect_true(emphasis:::.test_feasibility(
    c(1.0, 0.9), model = "cr", link = "linear",
    max_t = brts[1], max_lin = 600L, n_test = 8L,
    tip_lo = 3, tip_hi = 300, num_threads = 1L))
})

# Tests for the inference module (estimate_rates, compare_models, etc.)

# ---------- Internal helpers ----------

test_that(".par_names returns correct names for each model", {
  expect_equal(emphasis:::.par_names(c(0L, 0L, 0L)), c("beta_0", "gamma_0"))
  expect_equal(emphasis:::.par_names(c(1L, 0L, 0L)),
               c("beta_0", "beta_N", "gamma_0", "gamma_N"))
  expect_equal(emphasis:::.par_names(c(0L, 0L, 1L)),
               c("beta_0", "beta_D", "gamma_0", "gamma_D"))
  expect_equal(emphasis:::.par_names(c(1L, 0L, 1L)),
               c("beta_0", "beta_N", "beta_D", "gamma_0", "gamma_N", "gamma_D"))
})

test_that(".contract_pars and .expand_pars are inverses", {
  pars8 <- c(0.5, -0.01, 0, 0, 0.1, 0.02, 0, 0)
  model <- c(1L, 0L, 0L)
  compact <- emphasis:::.contract_pars(pars8, model)
  expect_equal(compact, c(0.5, -0.01, 0.1, 0.02))
  back <- emphasis:::.expand_pars(compact, model)
  expect_equal(back, pars8)
})

test_that(".contract_pars and .expand_pars round-trip on the partial layouts", {
  # The round trip in test-covariates.R uses mb = c(1, 1, 1), where every slot
  # is active and .expand_pars is the identity (audit id H102).  The layouts a
  # user actually fits are partial: the expansion has to put each free
  # parameter in its own slot and leave the rest at zero.
  layouts <- list(
    list(mb = c(1L, 0L, 0L), compact = c(0.5, -0.01, 0.1, 0.02),
         full = c(0.5, -0.01, 0, 0, 0.1, 0.02, 0, 0)),
    list(mb = c(0L, 0L, 1L), compact = c(0.5, 0.05, 0.1, 0.01),
         full = c(0.5, 0, 0, 0.05, 0.1, 0, 0, 0.01)),
    list(mb = c(1L, 0L, 1L), compact = c(0.5, -0.01, 0.05, 0.1, 0.02, 0.01),
         full = c(0.5, -0.01, 0, 0.05, 0.1, 0.02, 0, 0.01)),
    list(mb = c(0L, 1L, 0L), compact = c(0.5, 0.03, 0.1, -0.02),
         full = c(0.5, 0, 0.03, 0, 0.1, 0, -0.02, 0))
  )
  for (lay in layouts) {
    expect_equal(emphasis:::.expand_pars(lay$compact, lay$mb), lay$full)
    expect_equal(emphasis:::.contract_pars(lay$full, lay$mb), lay$compact)
    expect_length(emphasis:::.par_names(lay$mb), length(lay$compact))
    # the slots the model does not select stay zero whatever the input
    expect_true(all(emphasis:::.expand_pars(lay$compact, lay$mb)[
      c(1 + which(lay$mb == 0L), 5 + which(lay$mb == 0L))] == 0))
  }
})

test_that(".model_label returns correct labels", {
  expect_equal(emphasis:::.model_label(c(0L, 0L, 0L)), "CR")
  expect_equal(emphasis:::.model_label(c(1L, 0L, 0L)), "N")
  expect_equal(emphasis:::.model_label(c(0L, 0L, 1L)), "D")
  expect_equal(emphasis:::.model_label(c(1L, 0L, 1L)), "N + D")
})

test_that(".extract_brts handles different input types", {
  # Numeric vector
  brts <- c(5.0, 3.0, 1.0)
  expect_equal(emphasis:::.extract_brts(brts), c(5.0, 3.0, 1.0))

  # Unsorted numeric vector
  expect_equal(emphasis:::.extract_brts(c(1.0, 5.0, 3.0)), c(5.0, 3.0, 1.0))
})

test_that("estimate_rates_control returns expected keys", {
  ctrl_mcem <- emphasis:::estimate_rates_control("mcem")
  expect_true(all(c("lower_bound", "upper_bound", "sample_size",
                     "max_iter", "maxN") %in% names(ctrl_mcem)))

  ctrl_cem <- emphasis:::estimate_rates_control("cem")
  expect_true(all(c("lower_bound", "upper_bound", "num_particles",
                     "max_iter") %in% names(ctrl_cem)))

  ctrl_gam <- emphasis:::estimate_rates_control("gam")
  expect_true(all(c("lower_bound", "upper_bound", "grid_points",
                     "sample_size", "spline_type") %in% names(ctrl_gam)))
})

test_that("estimate_rates errors on missing bounds", {
  expect_error(
    estimate_rates(c(5, 3, 1), method = "mcem", model = "cr"),
    "lower_bound.*upper_bound.*must be supplied"
  )
})

test_that("estimate_rates errors on wrong bound length", {
  expect_error(
    estimate_rates(c(5, 3, 1), method = "mcem", model = "cr",
                   control = list(lower_bound = c(0, 0, 0),
                                  upper_bound = c(1, 1, 1))),
    "requires 2 parameters"
  )
})

# ---------- the three methods end to end (audit ids H94, H98, H104) ---------
#
# These blocks carried an unconditional skip() behind a skip_on_cran(), so the
# reason reported was "On CRAN" even locally, and their body built the tree
# with a positional simulate_tree(c(0.5, 0.1), ...) call, which binds the
# parameters to `tree`.  estimate_rates takes branching times directly, so a
# fixed vector removes both the stale call and the chance of a run drawing an
# extinct tree.

# Fixed 12-tip branching times, crown age 5, decreasing.
brts_cr <- c(5, 4.5, 4.0, 3.5, 3.0, 2.6, 2.2, 1.8, 1.4, 1.0, 0.5)

test_that("estimate_rates CR with mcem runs end-to-end", {
  fit <- estimate_rates(brts_cr, method = "mcem", model = "cr",
    control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                   max_iter = 3L, sample_size = 20L, num_threads = 1L))
  expect_s3_class(fit, "emphasis_fit")
  expect_equal(length(fit$pars), 2L)
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$method, "mcem")
  expect_equal(fit$model, c(0L, 0L, 0L))
  expect_equal(fit$AIC, -2 * fit$loglik + 2 * fit$n_pars)
  expect_true(all(fit$pars >= c(0, 0) - 1e-8))
  expect_true(all(fit$pars <= c(2, 1) + 1e-8))
})

test_that("estimate_rates CR with cem runs end-to-end", {
  fit <- estimate_rates(brts_cr, method = "cem", model = "cr",
    control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                   max_iter = 3L, num_particles = 20L, num_threads = 1L))
  expect_s3_class(fit, "emphasis_fit")
  expect_equal(length(fit$pars), 2L)
  expect_equal(fit$method, "cem")
  expect_true(all(fit$pars >= c(0, 0) - 1e-8))
  expect_true(all(fit$pars <= c(2, 1) + 1e-8))
})

test_that("estimate_rates CR with gam runs end-to-end", {
  skip_on_cran()                       # trains two GAMs; ~3 s
  skip_if_not_installed("mgcv")
  fit <- suppressMessages(estimate_rates(brts_cr, method = "gam", model = "cr",
    control = list(lower_bound = c(0.1, 0.01), upper_bound = c(1.5, 0.5),
                   grid_points = 8L, sample_size = 50L, num_threads = 1L)))
  expect_s3_class(fit, "emphasis_fit")
  expect_equal(length(fit$pars), 2L)
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$method, "gam")
  expect_true(all(fit$pars >= c(0.1, 0.01) - 1e-8))
  expect_true(all(fit$pars <= c(1.5, 0.5) + 1e-8))
})

test_that("the three methods agree on the same tree to within Monte Carlo error", {
  skip_on_cran()
  # A cr fit has a closed form: DDD::bd_loglik at the MLE.  Two Monte Carlo
  # methods on the same tree must land in the same region of the box, which is
  # what a user reading compare_models() relies on.  The tolerance is the width
  # of the box, not a pinned number.
  mc  <- estimate_rates(brts_cr, method = "mcem", model = "cr",
    control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                   max_iter = 5L, sample_size = 50L, num_threads = 1L))
  ce  <- estimate_rates(brts_cr, method = "cem", model = "cr",
    control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                   max_iter = 5L, num_particles = 30L, num_threads = 1L))
  expect_lt(abs(unname(mc$pars["beta_0"]) - unname(ce$pars["beta_0"])), 0.5)
  expect_true(is.finite(mc$loglik) && is.finite(ce$loglik))
  expect_lt(abs(mc$loglik - ce$loglik), 5)
})

test_that("compare_models errors with fewer than 2 fits", {
  fake_fit <- structure(list(pars = c(0.5, 0.1), loglik = -10,
    n_pars = 2L, AIC = 24, model = c(0L, 0L, 0L)), class = "emphasis_fit")
  expect_error(emphasis:::compare_models(fake_fit), "at least two")
})

test_that("compare_models produces correct table structure", {
  fake_cr <- structure(list(pars = c(beta_0 = 0.5, gamma_0 = 0.1),
    loglik = -10, loglik_var = NA_real_, n_pars = 2L, AIC = 24,
    model = c(0L, 0L, 0L)), class = "emphasis_fit")
  fake_dd <- structure(list(pars = c(beta_0 = 0.5, beta_N = -0.01,
    gamma_0 = 0.1, gamma_N = 0.005),
    loglik = -8, loglik_var = NA_real_, n_pars = 4L, AIC = 24,
    model = c(1L, 0L, 0L)), class = "emphasis_fit")
  tab <- emphasis:::compare_models(CR = fake_cr, DD = fake_dd)
  expect_true(all(c("model", "n_pars", "loglik", "AIC",
                     "delta_AIC", "AICw") %in% names(tab)))
  expect_equal(nrow(tab), 2L)
  expect_equal(min(tab$delta_AIC), 0)
})

test_that("print.emphasis_fit works", {
  fake_fit <- structure(list(pars = c(beta_0 = 0.5, gamma_0 = 0.1),
    loglik = -10, loglik_var = NA_real_, n_pars = 2L, AIC = 24,
    method = "mcem", model = c(0L, 0L, 0L)), class = "emphasis_fit")
  out <- capture.output(print(fake_fit))
  expect_true(any(grepl("emphasis fit", out)))
  expect_true(any(grepl("beta_0", out)))
})

# ---------- IS fhat denominator correction ----------

test_that(".is_fhat includes zero-weight trees in denominator", {
  # With 10 valid trees and 5 zero-weight trees:
  # fhat should use S_completed = 15, not 10
  set.seed(1)
  logf <- rnorm(10, -5, 1)
  logg <- rnorm(10, -5, 1)

  fhat_without <- emphasis:::.is_fhat(logf, logg, n_zero_weight = 0L)
  fhat_with    <- emphasis:::.is_fhat(logf, logg, n_zero_weight = 5L)

  # With more zero-weight trees in denominator, fhat should be lower
  expect_true(fhat_with < fhat_without)

  # The difference should be log(10/15)
  expected_diff <- log(10 / 15)
  expect_equal(fhat_with - fhat_without, expected_diff, tolerance = 1e-10)
})

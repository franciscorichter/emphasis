# Tests for the inference module (estimate_rates, compare_models, etc.)

# ---------- Internal helpers ----------

test_that(".par_names returns correct names for each model", {
  expect_equal(emphasis:::.par_names(c(0L, 0L, 0L)), c("beta_0", "gamma_0"))
  expect_equal(emphasis:::.par_names(c(1L, 0L, 0L)),
               c("beta_0", "beta_N", "gamma_0", "gamma_N"))
  expect_equal(emphasis:::.par_names(c(0L, 1L, 0L)),
               c("beta_0", "beta_P", "gamma_0", "gamma_P"))
  expect_equal(emphasis:::.par_names(c(0L, 0L, 1L)),
               c("beta_0", "beta_E", "gamma_0", "gamma_E"))
  expect_equal(emphasis:::.par_names(c(1L, 1L, 0L)),
               c("beta_0", "beta_N", "beta_P", "gamma_0", "gamma_N", "gamma_P"))
})

test_that(".contract_pars and .expand_pars are inverses", {
  pars8 <- c(0.5, -0.01, 0, 0, 0.1, 0.02, 0, 0)
  model <- c(1L, 0L, 0L)
  compact <- emphasis:::.contract_pars(pars8, model)
  expect_equal(compact, c(0.5, -0.01, 0.1, 0.02))
  back <- emphasis:::.expand_pars(compact, model)
  expect_equal(back, pars8)
})

test_that(".model_label returns correct labels", {
  expect_equal(emphasis:::.model_label(c(0L, 0L, 0L)), "CR")
  expect_equal(emphasis:::.model_label(c(1L, 0L, 0L)), "N")
  expect_equal(emphasis:::.model_label(c(0L, 1L, 0L)), "PD")
  expect_equal(emphasis:::.model_label(c(0L, 0L, 1L)), "EP")
  expect_equal(emphasis:::.model_label(c(1L, 1L, 0L)), "N + PD")
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

test_that("estimate_rates CR with mcem runs end-to-end", {
  skip_on_cran()
  skip("C++ integration test: run locally with devtools::test(filter='inference')")
  set.seed(42)
  sim <- simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr")
  fit <- estimate_rates(sim, method = "mcem", model = "cr",
    control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                   max_iter = 5, sample_size = 50))
  expect_s3_class(fit, "emphasis_fit")
  expect_equal(length(fit$pars), 2L)
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$method, "mcem")
})

test_that("estimate_rates CR with cem runs end-to-end", {
  skip_on_cran()
  skip("C++ integration test: run locally with devtools::test(filter='inference')")
  set.seed(42)
  sim <- simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr")
  fit <- estimate_rates(sim, method = "cem", model = "cr",
    control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                   max_iter = 5, num_particles = 20))
  expect_s3_class(fit, "emphasis_fit")
  expect_equal(length(fit$pars), 2L)
  expect_equal(fit$method, "cem")
})

test_that("estimate_rates CR with gam runs end-to-end", {
  skip_on_cran()
  skip("C++ integration test: run locally with devtools::test(filter='inference')")
  set.seed(42)
  sim <- simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr")
  fit <- estimate_rates(sim, method = "gam", model = "cr",
    control = list(lower_bound = c(0.1, 0.01), upper_bound = c(1.5, 0.5),
                   grid_points = 8, sample_size = 50))
  expect_s3_class(fit, "emphasis_fit")
  expect_equal(length(fit$pars), 2L)
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$method, "gam")
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

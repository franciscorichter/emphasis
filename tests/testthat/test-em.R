# Tests for EM / E-step functions
# C++ integration tests are slow and run locally only

test_that("mc_loglik returns a list with expected fields", {
  skip("C++ integration test: run locally with devtools::test(filter='em')")
  set.seed(42)
  brts <- c(0.8, 0.6, 0.4, 0.2)
  # 8-param layout: beta_0, beta_N, beta_P, beta_E, gamma_0, gamma_N, gamma_P, gamma_E
  pars8 <- c(0.5, -0.01, 0.01, 0, 0.1, 0, 0, 0)
  lb8 <- c(0, -0.1, -0.1, 0, 0, 0, 0, 0)
  ub8 <- c(2, 0.1, 0.1, 0, 0.5, 0, 0, 0)
  result <- mc_loglik(
    brts        = brts,
    pars        = pars8,
    sample_size = 10,
    maxN        = 100,
    max_missing = 1000,
    max_lambda  = 500,
    lower_bound = lb8,
    upper_bound = ub8,
    xtol_rel    = 1e-3,
    num_threads = 1,
    model       = c(1L, 1L, 0L)
  )

  expect_type(result, "list")
  expect_true("fhat" %in% names(result))
})

test_that("estimate_rates CR smoke test", {
  skip("C++ integration test: run locally with devtools::test(filter='em')")
  set.seed(42)
  tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
  fit <- estimate_rates(tr, method = "mcem", model = "cr",
    lower_bound = c(0, 0), upper_bound = c(2, 1),
    control = list(sample_size = 20, tol = 0.5, burnin = 2))
  expect_named(fit$pars, c("beta_0", "gamma_0"))
  expect_length(fit$pars, 2)
  expect_true(is.numeric(fit$loglik))
})

test_that("estimate_rates DD smoke test", {
  skip("C++ integration test: run locally with devtools::test(filter='em')")
  set.seed(42)
  tr <- simulate_tree(pars = c(0.5, -0.005, 0.1, 0), max_t = 8, model = "dd")
  fit <- estimate_rates(tr, method = "mcem", model = "dd",
    lower_bound = c(0.1, -0.1, 0, -0.01),
    upper_bound = c(2, 0.01, 0.5, 0.01),
    control = list(sample_size = 50, tol = 0.5, burnin = 5))
  expect_named(fit$pars, c("beta_0", "beta_N", "gamma_0", "gamma_N"))
  expect_length(fit$pars, 4)
})

test_that("estimate_rates EP smoke test", {
  skip("C++ integration test: run locally with devtools::test(filter='em')")
  set.seed(1)
  tr <- simulate_tree(pars = c(0.5, 0.05, 0.1, 0.01), max_t = 5, model = "ep")
  expect_equal(tr$status, "done")
  fit <- estimate_rates(tr, method = "mcem", model = "ep",
    lower_bound = c(0.1, -0.5, 0,  -0.5),
    upper_bound = c(2,    0.5, 0.5, 0.5),
    control = list(sample_size = 50, tol = 0.5, burnin = 3))
  expect_named(fit$pars, c("beta_0", "beta_E", "gamma_0", "gamma_E"))
  expect_true(is.numeric(fit$loglik))
})

test_that("estimate_rates EP + exponential link errors", {
  skip("C++ integration test: run locally with devtools::test(filter='em')")
  set.seed(1)
  tr <- simulate_tree(pars = c(0.5, 0.05, 0.1, 0.01), max_t = 5, model = "ep")
  expect_error(
    estimate_rates(tr, method = "mcem", model = "ep", link = "exponential",
      lower_bound = c(-5, -5, -5, -5),
      upper_bound = c(2, 2, 2, 2)),
    "EP model with exponential link not yet supported"
  )
})

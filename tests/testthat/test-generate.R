# Tests for rate / generation functions

# ---------- nhExpRand ----------

test_that("nhExpRand generates the correct number of variates", {
  rate_func <- function(t) 1 + 0.5 * sin(t)
  result <- nhExpRand(n = 1, rate_func = rate_func, now = 0, tMax = 5)
  expect_length(result, 1)
  expect_true(result >= 0)
  expect_true(result <= 5)
})

test_that("nhExpRand errors on invalid inputs", {
  expect_error(
    nhExpRand(1, function(t) t, tMax = Inf),
    "valid finite tMax"
  )
  expect_error(
    nhExpRand(1, 42),
    "rate_func must be a function"
  )
  expect_error(
    nhExpRand(1, function(t) t, now = -1, tMax = 5),
    "now must be greater"
  )
  expect_error(
    nhExpRand(0, function(t) t, tMax = 5),
    "'n' must be >= 1"
  )
})

test_that("nhExpRand returns NA with warning on impossible root-finding", {
  # A rate function that is essentially zero — uniroot can't find a crossing
  expect_warning(
    result <- nhExpRand(n = 1, rate_func = function(t) 1e-300, now = 0, tMax = 0.001),
    "could not be computed"
  )
  expect_length(result, 1)
})

# ---------- generateNonHomogeneousExp ----------

test_that("generateNonHomogeneousExp errors when max_time <= start_time", {
  expect_error(
    generateNonHomogeneousExp(1, function(t) t + 1, start_time = 5, max_time = 3),
    "max_time.*must be greater"
  )
})

test_that("generateNonHomogeneousExp errors on non-function rate_func", {
  expect_error(
    generateNonHomogeneousExp(1, 42, start_time = 0, max_time = 5),
    "'rate_func' must be a function"
  )
})

test_that("generateNonHomogeneousExp errors on non-positive rate", {
  expect_error(
    generateNonHomogeneousExp(1, function(t) -1, start_time = 0, max_time = 5),
    "non-positive"
  )
})

test_that("generateNonHomogeneousExp returns correct count", {
  rate_func <- function(t) 2 * t + 1
  result <- generateNonHomogeneousExp(3, rate_func, 0, 10)
  expect_length(result, 3)
  expect_true(all(result > 0 & result < 10))
})

# ---------- rate_t ----------

test_that("rate_t returns correct scalar", {
  cov_func1 <- function(t) sin(t)
  cov_func2 <- function(t) cos(t)
  params <- c(1, 0.5, -0.5)

  r <- rate_t(t = 0, params = params, cov_funcs = list(cov_func1, cov_func2))
  # At t=0: sin(0)=0, cos(0)=1 => 1*1 + 0.5*0 + (-0.5)*1 = 0.5
  expect_equal(r, 0.5)

  r_exp <- rate_t(
    t = 0, params = params,
    cov_funcs = list(cov_func1, cov_func2),
    use_exponential = TRUE
  )
  expect_equal(r_exp, exp(0.5))
})

test_that("rate_t errors on params length mismatch", {
  expect_error(
    rate_t(t = 0, params = c(1, 2, 3), cov_funcs = list(function(t) t)),
    "must have length 2"
  )
})

test_that("rate_t errors when cov_funcs is not a list", {
  expect_error(
    rate_t(t = 0, params = c(1, 2), cov_funcs = function(t) t),
    "'cov_funcs' must be a list"
  )
})

# ---------- ExponentialRate ----------

test_that("ExponentialRate errors with mismatched lengths", {
  expect_error(
    ExponentialRate(matrix(c(1, 2), ncol = 1), c(0.5, 1, 2)),
    "length.*parameters.*must be one more"
  )
})

test_that("ExponentialRate returns numeric", {
  cov <- matrix(c(1, 2, 3, 4), ncol = 2)
  pars <- c(0.5, 1, -0.5)
  result <- ExponentialRate(cov, pars)
  expect_length(result, 1)
  expect_true(is.numeric(result))
})

test_that("ExponentialRate computes exp(bias + beta.x) per row", {
  cov <- matrix(c(1, 0), ncol = 1) # two rows: [1] and [0]
  pars <- c(0, 1) # bias=0, beta=1
  result <- ExponentialRate(cov, pars)
  # Row 1: exp(0 + 1*1) = e; Row 2: exp(0 + 1*0) = 1
  expect_equal(result, exp(1) + 1, tolerance = 1e-10)
})

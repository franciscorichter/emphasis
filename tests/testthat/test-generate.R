# Tests for rate / generation functions

test_that("nhExpRand generates the correct number of variates", {
  rate_func <- function(t) 1 + 0.5 * sin(t)
  result <- nhExpRand(n = 1, rate_func = rate_func, now = 0, tMax = 5)
  expect_length(result, 1)
  expect_true(result >= 0)
  expect_true(result <= 5)
})

test_that("nhExpRand errors on invalid inputs", {
  expect_error(nhExpRand(1, function(t) t, tMax = Inf),
               "Need a valid tMax")
  expect_error(nhExpRand(1, 42),
               "rate_func must be a function")
  expect_error(nhExpRand(1, function(t) t, now = -1, tMax = 5),
               "now must be greater")
})

test_that("generateNonHomogeneousExp errors when max_time <= start_time", {
  expect_error(
    generateNonHomogeneousExp(1, function(t) t + 1, start_time = 5, max_time = 3),
    "max_time.*must be greater"
  )
})

test_that("rate_t returns correct scalar", {
  cov_func1 <- function(t) sin(t)
  cov_func2 <- function(t) cos(t)
  params <- c(1, 0.5, -0.5)

  r <- rate_t(t = 0, params = params, cov_funcs = list(cov_func1, cov_func2))
  # At t=0: sin(0)=0, cos(0)=1 => 1*1 + 0.5*0 + (-0.5)*1 = 0.5
  expect_equal(r, 0.5)

  r_exp <- rate_t(t = 0, params = params,
                  cov_funcs = list(cov_func1, cov_func2),
                  use_exponential = TRUE)
  expect_equal(r_exp, exp(0.5))
})

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

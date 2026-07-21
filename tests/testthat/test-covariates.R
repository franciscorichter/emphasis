# Tests for the orthogonal covariate basis {N, M = P/N, D = E - M}.
# These lock the R-layer resolution, naming, and labelling of the basis.

test_that("model shortcuts resolve to the correct covariate slots", {
  expect_equal(emphasis:::.resolve_model("cr"), c(0L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model("dd"), c(1L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model("ma"), c(0L, 1L, 0L))  # mean age M
  expect_equal(emphasis:::.resolve_model("rd"), c(0L, 0L, 1L))  # deviation D
})

test_that("legacy shortcuts pd/ep remain aliases for ma/rd", {
  expect_equal(emphasis:::.resolve_model("pd"), emphasis:::.resolve_model("ma"))
  expect_equal(emphasis:::.resolve_model("ep"), emphasis:::.resolve_model("rd"))
})

test_that("model formulas parse N/M/D and legacy aliases", {
  expect_equal(emphasis:::.resolve_model(~ N + M), c(1L, 1L, 0L))
  expect_equal(emphasis:::.resolve_model(~ N + D), c(1L, 0L, 1L))
  expect_equal(emphasis:::.resolve_model(~ M + D), c(0L, 1L, 1L))
  # legacy aliases still accepted
  expect_equal(emphasis:::.resolve_model(~ N + PD), c(1L, 1L, 0L))
  expect_equal(emphasis:::.resolve_model(~ N + EP), c(1L, 0L, 1L))
})

test_that("unknown covariates in a formula are rejected", {
  expect_error(emphasis:::.resolve_model(~ N + Q), "Unknown covariate")
})

test_that("parameter names follow the {N, M, D} basis", {
  expect_equal(emphasis:::.par_names(c(0L, 0L, 0L)),
               c("beta_0", "gamma_0"))
  expect_equal(emphasis:::.par_names(c(0L, 1L, 0L)),
               c("beta_0", "beta_M", "gamma_0", "gamma_M"))
  expect_equal(emphasis:::.par_names(c(0L, 0L, 1L)),
               c("beta_0", "beta_D", "gamma_0", "gamma_D"))
  expect_equal(emphasis:::.par_names(c(1L, 1L, 1L)),
               c("beta_0", "beta_N", "beta_M", "beta_D",
                 "gamma_0", "gamma_N", "gamma_M", "gamma_D"))
})

test_that("model labels use the N/M/D covariate names", {
  expect_equal(emphasis:::.model_label(c(0L, 0L, 0L)), "CR")
  expect_equal(emphasis:::.model_label(c(1L, 0L, 0L)), "N")
  expect_equal(emphasis:::.model_label(c(0L, 1L, 0L)), "M")
  expect_equal(emphasis:::.model_label(c(0L, 0L, 1L)), "D")
  expect_equal(emphasis:::.model_label(c(1L, 1L, 1L)), "N + M + D")
})

test_that("expand/contract round-trips under the new basis", {
  mb <- c(1L, 1L, 1L)
  compact <- c(0.5, 0.02, 0.03, 0.01, 0.1, -0.01, 0.02, 0.005)
  full8 <- emphasis:::.expand_pars(compact, mb)
  expect_length(full8, 8L)
  expect_equal(emphasis:::.contract_pars(full8, mb), compact)
})

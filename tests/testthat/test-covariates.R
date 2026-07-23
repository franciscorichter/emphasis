# Tests for the canonical two-covariate basis {N, D} (D = E - M age-imbalance).
# M (slot 2) is retained internally only to centre D and is not user-selectable.

test_that("model shortcuts resolve to the correct covariate slots", {
  expect_equal(emphasis:::.resolve_model("cr"), c(0L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model("dd"), c(1L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model("d"),  c(0L, 0L, 1L))  # age-imbalance D
  expect_equal(emphasis:::.resolve_model("nd"), c(1L, 0L, 1L))  # N + D (recommended)
})

test_that("legacy shortcuts ep/rd remain aliases for d", {
  expect_equal(emphasis:::.resolve_model("ep"), emphasis:::.resolve_model("d"))
  expect_equal(emphasis:::.resolve_model("rd"), emphasis:::.resolve_model("d"))
})

test_that("model formulas parse N and D (and the legacy EP alias)", {
  expect_equal(emphasis:::.resolve_model(~ N),      c(1L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model(~ N + D),  c(1L, 0L, 1L))
  expect_equal(emphasis:::.resolve_model(~ D),      c(0L, 0L, 1L))
  expect_equal(emphasis:::.resolve_model(~ N + EP), c(1L, 0L, 1L))  # legacy alias
})

test_that("M is no longer a user-selectable covariate", {
  expect_error(emphasis:::.resolve_model("ma"))                    # removed shortcut
  expect_error(emphasis:::.resolve_model(~ N + M), "Unknown covariate")
  expect_error(emphasis:::.resolve_model(~ N + Q), "Unknown covariate")
})

test_that("parameter names follow the {N, D} basis", {
  expect_equal(emphasis:::.par_names(c(0L, 0L, 0L)),
               c("beta_0", "gamma_0"))
  expect_equal(emphasis:::.par_names(c(1L, 0L, 0L)),
               c("beta_0", "beta_N", "gamma_0", "gamma_N"))
  expect_equal(emphasis:::.par_names(c(0L, 0L, 1L)),
               c("beta_0", "beta_D", "gamma_0", "gamma_D"))
  expect_equal(emphasis:::.par_names(c(1L, 0L, 1L)),
               c("beta_0", "beta_N", "beta_D",
                 "gamma_0", "gamma_N", "gamma_D"))
  # internal length-3 slot machinery still labels the full vector (M included)
  expect_equal(emphasis:::.par_names(c(1L, 1L, 1L)),
               c("beta_0", "beta_N", "beta_M", "beta_D",
                 "gamma_0", "gamma_N", "gamma_M", "gamma_D"))
})

test_that("model labels use the N/D covariate names", {
  expect_equal(emphasis:::.model_label(c(0L, 0L, 0L)), "CR")
  expect_equal(emphasis:::.model_label(c(1L, 0L, 0L)), "N")
  expect_equal(emphasis:::.model_label(c(0L, 0L, 1L)), "D")
  expect_equal(emphasis:::.model_label(c(1L, 0L, 1L)), "N + D")
})

test_that("expand/contract round-trips (internal length-3 machinery intact)", {
  mb <- c(1L, 1L, 1L)
  compact <- c(0.5, 0.02, 0.03, 0.01, 0.1, -0.01, 0.02, 0.005)
  full8 <- emphasis:::.expand_pars(compact, mb)
  expect_length(full8, 8L)
  expect_equal(emphasis:::.contract_pars(full8, mb), compact)
})

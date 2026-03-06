# Tests for simulation functions

# ---------- Input validation for simulate_tree ----------

test_that("simulate_tree errors on invalid pars", {
  expect_error(
    simulate_tree(pars = "x", max_t = 5),
    "non-empty numeric"
  )
  expect_error(
    simulate_tree(pars = numeric(0), max_t = 5),
    "non-empty numeric"
  )
})

test_that("simulate_tree errors on invalid max_t", {
  expect_error(
    simulate_tree(pars = c(0.5, 0.1), max_t = -1),
    "positive number"
  )
  expect_error(
    simulate_tree(pars = c(0.5, 0.1), max_t = 0),
    "positive number"
  )
})

test_that("simulate_tree errors on invalid max_lin", {
  expect_error(
    simulate_tree(pars = c(0.5, 0.1), max_t = 5, max_lin = 0),
    "positive number"
  )
})

test_that("simulate_tree accepts max_tries = 0 (no retries)", {
  skip("C++ integration test: run locally with devtools::test(filter='simulate')")
  set.seed(1)
  # max_tries = 0 is valid: extinct runs are not retried (used for GAM training)
  result <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, max_tries = 0,
                          useDDD = FALSE)
  expect_true(result$status %in% c("done", "extinct", "too_large"))
})

test_that("simulate_tree errors on wrong pars length for model", {
  # "dd" = c(1,0,0) needs 4 pars, supplying 2 should error
  expect_error(
    simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "dd"),
    "requires 4 parameters"
  )
})

test_that("simulate_tree errors on negative max_t", {
  expect_error(
    simulate_tree(pars = c(0.5, 0.1), max_t = -1),
    "positive number"
  )
})

test_that(".resolve_model handles shortcuts and binary vectors", {
  expect_equal(emphasis:::.resolve_model("cr"), c(0L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model("dd"), c(1L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model("pd"), c(0L, 1L, 0L))
  expect_equal(emphasis:::.resolve_model("ep"), c(0L, 0L, 1L))
  expect_equal(emphasis:::.resolve_model(c(1, 1, 0)), c(1L, 1L, 0L))
  expect_error(
    emphasis:::.resolve_model(c(1, 2, 0)),
    "length-3 binary integer vector"
  )
})

test_that(".expand_pars errors on wrong length", {
  expect_error(
    emphasis:::.expand_pars(c(1, 2), c(1L, 0L, 0L)),
    "requires 4 parameters"
  )
})

# ---------- C++ simulation tests (run locally) ----------

test_that("simulate_tree cr model returns valid output", {
  skip("C++ integration test: run locally with devtools::test(filter='simulate')")

  set.seed(42)
  result <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")

  expect_type(result, "list")
  expect_named(result, c("tes", "tas", "L", "status", "survival_prob"))
  expect_true(result$status %in% c("done", "extinct", "too_large"))
})

test_that("simulate_tree dd model returns valid output", {
  skip("C++ integration test: run locally with devtools::test(filter='simulate')")

  set.seed(42)
  result <- simulate_tree(
    pars = c(0.8, -0.02, 0.2, -0.01),
    max_t = 5, model = "dd"
  )

  expect_type(result, "list")
  expect_true(result$status %in% c("done", "extinct", "too_large"))
})

test_that("simulate_tree returns survival_prob", {
  skip("C++ integration test: run locally with devtools::test(filter='simulate')")

  set.seed(42)
  result <- simulate_tree(pars = c(0.8, 0.1), max_t = 5, model = "cr")
  expect_true(is.numeric(result$survival_prob))
  expect_true(result$survival_prob >= 0 && result$survival_prob <= 1)
})

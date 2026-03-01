# Tests for Differential Evolution fitting

test_that("emphasis_de errors on mismatched bounds", {
  expect_error(
    emphasis_de(1:3, 1, 10, 100,
                sd_vec = c(0.1), lower_bound = c(0, 0), upper_bound = c(1),
                max_lambda = 10),
    "same length"
  )
})

test_that("emphasis_de integration test", {
  skip("C++ integration test: run locally with devtools::test(filter='de')")
  set.seed(42)
  brts <- c(0.9, 0.7, 0.5, 0.3, 0.1)

  result <- emphasis_de(
    brts           = brts,
    num_iterations = 2,
    num_points     = 10,
    max_missing    = 100,
    sd_vec         = c(0.1, 0.3, 0.02, 0.02),
    lower_bound    = c(0, 0, -0.1, -0.1),
    upper_bound    = c(0.5, 2, 0.1, 0.1),
    maxN           = 20,
    max_lambda     = 100,
    disc_prop      = 0.5,
    verbose        = FALSE,
    num_threads    = 1
  )

  expect_type(result, "list")
  expected_names <- c("parameters", "time", "fhatdiff", "minloglik",
                      "meanloglik", "min_pars", "mean_pars", "obtained_estim", "times")
  expect_true(all(expected_names %in% names(result)))
  expect_length(result$obtained_estim, 4)
})

test_that("emphasis_de_factorial integration test", {
  skip("C++ integration test: run locally with devtools::test(filter='de')")
  set.seed(1)
  brts <- c(0.9, 0.7, 0.5, 0.3, 0.1)

  result <- emphasis_de_factorial(
    brts           = brts,
    num_iterations = 2,
    num_points     = 10,
    max_missing    = 100,
    sd_vec         = c(0.1, 0.3, 0.02, 0.02),
    lower_bound    = c(0, 0, -0.1, -0.1),
    upper_bound    = c(0.5, 2, 0.1, 0.1),
    maxN           = 20,
    max_lambda     = 100,
    verbose        = FALSE,
    num_threads    = 1
  )

  expect_type(result, "list")
  expect_true("logf_grid" %in% names(result))
  expect_true("logg_grid" %in% names(result))
})

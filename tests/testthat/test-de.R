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

  # 8-param layout for CR model
  lb8 <- c(0, 0, 0, 0, 0, 0, 0, 0)
  ub8 <- c(2, 0, 0, 0, 0.5, 0, 0, 0)
  sd8 <- c(0.3, 0, 0, 0, 0.1, 0, 0, 0)
  result <- emphasis_de(
    brts           = brts,
    num_iterations = 2,
    num_points     = 10,
    max_missing    = 100,
    sd_vec         = sd8,
    lower_bound    = lb8,
    upper_bound    = ub8,
    maxN           = 20,
    max_lambda     = 100,
    disc_prop      = 0.5,
    verbose        = FALSE,
    num_threads    = 1,
    model          = c(0L, 0L, 0L)
  )

  expect_type(result, "list")
  expected_names <- c("parameters", "time", "fhatdiff", "minloglik",
                      "meanloglik", "min_pars", "mean_pars", "obtained_estim", "times")
  expect_true(all(expected_names %in% names(result)))
  expect_length(result$obtained_estim, 8)
})

test_that("emphasis_de_factorial integration test", {
  skip("C++ integration test: run locally with devtools::test(filter='de')")
  set.seed(1)
  brts <- c(0.9, 0.7, 0.5, 0.3, 0.1)

  lb8 <- c(0, 0, 0, 0, 0, 0, 0, 0)
  ub8 <- c(2, 0, 0, 0, 0.5, 0, 0, 0)
  sd8 <- c(0.3, 0, 0, 0, 0.1, 0, 0, 0)
  result <- emphasis_de_factorial(
    brts           = brts,
    num_iterations = 2,
    num_points     = 10,
    max_missing    = 100,
    sd_vec         = sd8,
    lower_bound    = lb8,
    upper_bound    = ub8,
    maxN           = 20,
    max_lambda     = 100,
    verbose        = FALSE,
    num_threads    = 1,
    model          = c(0L, 0L, 0L)
  )

  expect_type(result, "list")
  expect_true("logf_grid" %in% names(result))
  expect_true("logg_grid" %in% names(result))
})

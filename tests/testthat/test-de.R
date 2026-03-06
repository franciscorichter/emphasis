# Tests for Monte Carlo Cross-Entropy Method (CEM) optimiser

test_that("emphasis_cem errors on mismatched bounds", {
  expect_error(
    emphasis_cem(brts = 1:3, max_iter = 1, num_points = 10, max_missing = 100,
                 sd_vec = c(0.1), lower_bound = c(0, 0), upper_bound = c(1),
                 max_lambda = 10),
    "same length"
  )
})

test_that("emphasis_cem integration test", {
  skip("C++ integration test: run locally with devtools::test(filter='de')")
  set.seed(42)
  brts <- c(0.9, 0.7, 0.5, 0.3, 0.1)

  # 8-param layout for CR model
  lb8 <- c(0, 0, 0, 0, 0, 0, 0, 0)
  ub8 <- c(2, 0, 0, 0, 0.5, 0, 0, 0)
  sd8 <- c(0.3, 0, 0, 0, 0.1, 0, 0, 0)
  result <- emphasis_cem(
    brts        = brts,
    max_iter    = 2,
    num_points  = 10,
    max_missing = 100,
    sd_vec      = sd8,
    lower_bound = lb8,
    upper_bound = ub8,
    maxN        = 20,
    max_lambda  = 100,
    disc_prop   = 0.5,
    verbose     = FALSE,
    num_threads = 1,
    model       = c(0L, 0L, 0L)
  )

  expect_type(result, "list")
  expect_named(result, c("best_loglik", "best_pars", "obtained_estim",
                         "loglik_var", "converged", "time"))
  expect_length(result$obtained_estim, 8)
  expect_true(result$converged %in% c("annealing", "plateau", "max_iter"))
  expect_true(all(result$best_loglik <= 0 | is.na(result$best_loglik)))
})

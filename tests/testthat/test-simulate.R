# Tests for simulation functions
# C++ simulation tests run locally only (skip in all automated checks)

test_that("sim_tree_pd_cpp returns valid output", {
  skip("C++ integration test: run locally with devtools::test(filter='simulate')")

  set.seed(42)
  result <- sim_tree_pd_cpp(
    pars     = c(0.05, 0.8, 0, 0),
    max_t    = 1,
    max_lin  = 1e4,
    max_tries = 50,
    useDDD   = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("tes", "tas", "L", "brts", "status"))
  expect_true(result$status %in% c("done", "extinct", "too_large"))
})

test_that("sim_tree_is_extinct_pd returns a data frame with correct columns", {
  skip("C++ integration test: run locally with devtools::test(filter='simulate')")

  set.seed(7)
  result <- sim_tree_is_extinct_pd(
    pars     = c(0.05, 0.8, 0, 0),
    max_t    = 1,
    num_repl = 3,
    max_lin  = 1e4
  )

  expect_s3_class(result, "data.frame")
  expect_named(result, c("is_extinct", "t", "N", "P", "break_condition"))
  expect_equal(nrow(result), 3)
  expect_true(all(result$break_condition %in% c("none", "max_t", "extinction", "max_lin")))
})

test_that("sim_tree_pd_grid returns a data frame", {
  skip("C++ integration test: run locally with devtools::test(filter='simulate')")

  result <- sim_tree_pd_grid(
    mu_vec     = c(0.05, 0.1),
    lambda_vec = c(0.5, 0.8),
    b_n_vec    = c(0, 0),
    b_p_vec    = c(0, 0),
    max_t      = 1,
    num_repl   = 2,
    max_N      = 1e4
  )

  expect_s3_class(result, "data.frame")
  expect_named(result, c("is_extinct", "t", "N", "P"))
})

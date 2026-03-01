# Tests for simulation functions
# C++ backend tests only run in R CMD check (not devtools::test())
# because devtools::load_all() can cause C++ library loading issues locally.
in_r_cmd_check <- function() {
  nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
}

test_that("sim_tree_pd_cpp returns valid output", {
  skip_on_cran()
  skip_if(!in_r_cmd_check(), "C++ test: run via R CMD check")

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
  skip_on_cran()
  skip_if(!in_r_cmd_check(), "C++ test: run via R CMD check")

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
  skip_on_cran()
  skip_if(!in_r_cmd_check(), "C++ test: run via R CMD check")

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

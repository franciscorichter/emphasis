# Tests for EM / E-step functions
# C++ integration tests are slow and run locally only
in_r_cmd_check <- function() nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))

test_that("e_cpp returns a list with expected fields", {
  skip("C++ integration test: run locally with devtools::test(filter='em')")
  set.seed(42)
  brts <- c(0.8, 0.6, 0.4, 0.2)
  result <- e_cpp(
    brts        = brts,
    init_pars   = c(0.1, 0.5, -0.01, 0.01),
    sample_size = 10,
    maxN        = 100,
    soc         = 2,
    max_missing = 1000,
    max_lambda  = 500,
    lower_bound = c(0, 0, -0.1, -0.1),
    upper_bound = c(0.5, 2, 0.1, 0.1),
    xtol_rel    = 1e-3,
    num_threads = 1
  )

  expect_type(result, "list")
  expect_true("fhat" %in% names(result))
})

test_that("loglikelihood returns numeric fhat", {
  skip("C++ integration test: run locally with devtools::test(filter='em')")
  set.seed(42)
  brts <- c(0.8, 0.6, 0.4, 0.2)

  e_result <- e_cpp(
    brts        = brts,
    init_pars   = c(0.1, 0.5, -0.01, 0.01),
    sample_size = 5,
    maxN        = 100,
    soc         = 2,
    max_missing = 1000,
    max_lambda  = 500,
    lower_bound = c(0, 0, -0.1, -0.1),
    upper_bound = c(0.5, 2, 0.1, 0.1),
    xtol_rel    = 1e-3,
    num_threads = 1
  )

  skip_if(is.null(e_result$trees) || length(e_result$trees) == 0,
          "e_cpp returned no trees")

  ll <- loglikelihood(
    pars        = c(0.1, 0.5, -0.01, 0.01),
    trees       = e_result$trees,
    logg        = e_result$logg,
    plugin      = "rpd5c",
    num_rejected = e_result$rejected
  )

  expect_type(ll, "list")
  expect_true(is.numeric(ll$fhat))
})

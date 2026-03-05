# Tests for mc_loglik (Monte Carlo log-likelihood)

test_that("mc_loglik returns a list on a small tree", {
  skip_on_cran()
  skip("mc_loglik integration test skipped locally to avoid long runtime")
  set.seed(42)
  tree <- ape::rphylo(8, 0.5, 0)
  brts <- ape::branching.times(tree)

  result <- mc_loglik(
    brts        = brts,
    pars        = c(0.1, 0.5, -0.01, 0.01),
    sample_size = 1,
    maxN        = 50,

    max_missing = 500,
    max_lambda  = 50000,
    lower_bound = c(0, 0, -0.1, -0.1),
    upper_bound = c(0.5, 2, 0.1, 0.1),
    xtol_rel    = 1e-5,
    num_threads = 1
  )

  expect_type(result, "list")
  expect_true("fhat" %in% names(result))
})

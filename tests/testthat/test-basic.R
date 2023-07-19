context("basic utility")

test_that("usage", {
  # create reconstructed bd tree
  
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.4, fossils = FALSE)
  brts <- ape::branching.times(focal_tree)

  max_missing <- 5000
  disc_prop <- 0.5 # DE Mutation proportion
  num_it <- 20 # Number of iterations
  sd_vec <- c(0.1, 0.1, 0.05, 0.01) # Initial parameter sampling variation
  lower_bound <- c(0, 0, 0.0, 0.0)
  upper_bound <-  c(1, 2, 0.0, 0.0)
  max_lambda <- 1000
  ax <- emphasis::emphasis_de(brts = brts,
                              num_iterations = num_it,
                              num_points = 1000,
                              max_missing = 1000,
                              sd_vec = sd_vec,
                              lower_bound = lower_bound,
                              upper_bound = upper_bound,
                              maxN = 10,
                              max_lambda = 1000,
                              disc_prop = disc_prop,
                              verbose = TRUE,
                              num_threads = 4) # only one for CRAN/CI

  testthat::expect_equal(length(ax$parameters), num_it)
  emphasis_estim <- tail(ax$min_pars, 1)
  # next two parameters were not fitted:
  testthat::expect_equal(emphasis_estim[[3]], 0)
  testthat::expect_equal(emphasis_estim[[4]], 0)
  
  if (requireNamespace("DDD")) {
  
    testthat::expect_output(
      ml_estim <- DDD::bd_ML(brts = brts)
    )

    div_ml <- ml_estim$lambda0 - ml_estim$mu0
    div_em <- emphasis_estim[[2]] - emphasis_estim[[1]]
  
    testthat::expect_equal(div_ml, div_em, tolerance = 0.1)

    ll_ml <- ml_estim$loglik
    ll_em <- -tail(ax$minloglik, 1)
    testthat::expect_equal(ll_ml, ll_em, tolerance = 1)
  }
})

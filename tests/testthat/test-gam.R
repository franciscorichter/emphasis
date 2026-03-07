# Tests for GAM-based methods (train_GAM, train_likelihood_GAM, etc.)

test_that("train_GAM errors when all columns are constant", {
  sims <- replicate(10, list(status = "done"), simplify = FALSE)
  pars_mat <- matrix(0.5, nrow = 10, ncol = 2)
  colnames(pars_mat) <- c("beta_0", "gamma_0")
  expect_error(
    train_GAM(sims, pars_mat, model = "cr"),
    "All parameter columns are constant"
  )
})

test_that("train_GAM fits a GAM with varying predictors", {
  skip_if_not_installed("mgcv")
  set.seed(1)
  n <- 50
  pars_mat <- cbind(beta_0 = runif(n, 0.2, 1.0),
                    gamma_0 = runif(n, 0.05, 0.3))
  # Simulate survival: higher beta_0 -> more likely to survive
  survived <- rbinom(n, 1, plogis(2 * pars_mat[, "beta_0"] - 1))
  sims <- lapply(survived, function(s) {
    list(status = if (s == 1) "done" else "extinct")
  })

  gam_fit <- suppressMessages(train_GAM(sims, pars_mat, model = "cr"))
  expect_s3_class(gam_fit, "gam")
})

test_that("predict_survival returns correct length", {
  skip_if_not_installed("mgcv")
  set.seed(1)
  n <- 50
  pars_mat <- cbind(beta_0 = runif(n, 0.2, 1.0),
                    gamma_0 = runif(n, 0.05, 0.3))
  survived <- rbinom(n, 1, 0.5)
  sims <- lapply(survived, function(s) {
    list(status = if (s == 1) "done" else "extinct")
  })

  gam_fit <- suppressMessages(train_GAM(sims, pars_mat, model = "cr"))
  newpars <- cbind(beta_0 = c(0.3, 0.6, 0.9), gamma_0 = c(0.1, 0.15, 0.2))
  pred <- predict_survival(gam_fit, newpars)
  expect_length(pred, 3L)
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("train_likelihood_GAM errors with too few valid points", {
  skip_if_not_installed("mgcv")
  surface <- data.frame(beta_0 = 1:5, gamma_0 = 1:5,
                        fhat = c(NA, NA, NA, NA, -10),
                        n_trees = rep(50, 5))
  expect_error(
    suppressMessages(emphasis:::train_likelihood_GAM(surface)),
    "Too few valid"
  )
})

test_that("train_likelihood_GAM fits with enough valid points", {
  skip_if_not_installed("mgcv")
  set.seed(1)
  n <- 20
  surface <- data.frame(
    beta_0 = seq(0.1, 1.5, length.out = n),
    gamma_0 = seq(0.05, 0.5, length.out = n),
    fhat = rnorm(n, -10, 2),
    n_trees = rep(50L, n)
  )
  gam_fit <- suppressMessages(
    emphasis:::train_likelihood_GAM(surface)
  )
  expect_s3_class(gam_fit, "gam")
})

test_that("find_MLE returns MLE within bounds", {
  skip_if_not_installed("mgcv")
  set.seed(1)
  # Create a known quadratic surface: peak at beta_0=0.5
  n <- 30
  x <- seq(0.1, 1.0, length.out = n)
  surface <- data.frame(
    beta_0 = x,
    fhat = -5 * (x - 0.5)^2 - 3,
    n_trees = rep(50L, n)
  )
  gam_fit <- suppressMessages(
    emphasis:::train_likelihood_GAM(surface, par_names = "beta_0")
  )
  mle <- emphasis:::find_MLE(gam_fit,
    lower_bound = 0.1, upper_bound = 1.0, par_names = "beta_0")
  expect_true(mle$pars[["beta_0"]] > 0.3 && mle$pars[["beta_0"]] < 0.7)
  expect_equal(mle$convergence, 0L)
})

test_that("diagnose_gam runs without error on a GAM fit", {
  skip_if_not_installed("mgcv")
  set.seed(1)
  n <- 20
  x <- seq(0.1, 1.0, length.out = n)
  surface <- data.frame(
    beta_0 = x,
    fhat = -5 * (x - 0.5)^2 - 3,
    n_trees = rep(50L, n)
  )
  gam_fit <- suppressMessages(
    emphasis:::train_likelihood_GAM(surface, par_names = "beta_0")
  )

  # Create a fake emphasis_fit with GAM details
  fake_fit <- structure(list(
    pars = c(beta_0 = 0.5),
    loglik = -3,
    loglik_var = NA_real_,
    n_pars = 1L,
    AIC = 8,
    method = "gam",
    model = c(0L, 0L, 0L),
    details = list(
      surface = surface,
      gam_fit = gam_fit,
      par_names = "beta_0"
    )
  ), class = "emphasis_fit")

  d <- diagnose_gam(fake_fit, plot = FALSE)
  expect_s3_class(d, "gam_diagnostics")
  expect_true(d$surface_summary$n_valid == n)
  expect_true(d$gam_summary$deviance_explained > 0)

  out <- capture.output(print(d))
  expect_true(any(grepl("GAM-based MLE", out)))
})

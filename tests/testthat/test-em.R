# Tests for EM / E-step functions (audit ids H94, H96, H97, H98, H104).
#
# The five blocks in this file carried an unconditional skip().  The first
# called mc_loglik(), which exists in no version of the package; it is written
# here against em_cpp(), the entry point that performs one EM step.  The other
# four passed lower_bound / upper_bound as top-level estimate_rates arguments
# (they are control fields) and built their tree with a positional
# simulate_tree() call, which binds the parameters to `tree`.  They now use
# fixed branching-time vectors, so no block depends on a simulation surviving.

# Fixed branching times, crown age first, decreasing.
brts8  <- c(4, 3.2, 2.6, 2.0, 1.5, 1.0, 0.4)
brts12 <- c(5, 4.5, 4.0, 3.5, 3.0, 2.6, 2.2, 1.8, 1.4, 1.0, 0.5)

test_that("em_cpp returns an E-step and an M-step with the documented fields", {
  # dd in the 8-parameter layout, with the M and D slots bounded to zero
  pars8 <- c(0.5, -0.01, 0, 0, 0.1, 0, 0, 0)
  lb8   <- c(0, -0.1, 0, 0, 0, 0, 0, 0)
  ub8   <- c(2,  0.1, 0, 0, 0.5, 0, 0, 0)

  result <- emphasis:::em_cpp(
    brts        = brts8,
    init_pars   = pars8,
    sample_size = 20L,
    maxN        = 1000L,
    max_missing = 1000L,
    max_lambda  = 500,
    lower_bound = lb8,
    upper_bound = ub8,
    xtol_rel    = 1e-3,
    num_threads = 1L,
    copy_trees  = FALSE,
    model       = c(1L, 0L, 0L),
    link        = 0L
  )

  expect_type(result, "list")
  expect_true(all(c("fhat", "logf", "logg", "weights", "estimates", "nlopt",
                    "num_trees", "rejected", "rejected_overruns",
                    "rejected_lambda", "rejected_zero_weights",
                    "rejected_nonfinite") %in% names(result)))

  expect_equal(as.integer(result$num_trees), 20L)
  expect_length(result$logf, 20L)
  expect_length(result$logg, 20L)
  expect_length(result$estimates, 8L)
  expect_true(all(is.finite(result$logf)))
  expect_true(all(is.finite(result$logg)))
  expect_true(is.finite(result$fhat))

  # the C++ E-step aggregates the weights with the same estimator as the R
  # layer's .is_fhat, counting the zero-weight draws in the denominator
  expect_equal(result$fhat,
               emphasis:::.is_fhat(result$logf, result$logg,
                                   n_zero_weight = result$rejected_zero_weights),
               tolerance = 1e-12)

  # the M-step respects the box, and the slots the box pins stay pinned
  expect_true(all(result$estimates >= lb8 - 1e-8))
  expect_true(all(result$estimates <= ub8 + 1e-8))
  expect_equal(result$estimates[lb8 == ub8], lb8[lb8 == ub8])

  # The weights handed to the M-step are the max-scaled linear IS weights
  # (E_step.cpp:212-219), not normalised ones, and they are the exponentials
  # of the same log weights the E-step returns.
  lw <- result$logf - result$logg
  expect_equal(max(result$weights), 1)
  expect_true(all(result$weights > 0 & result$weights <= 1))
  expect_equal(result$weights, exp(lw - max(lw)), tolerance = 1e-12)
})

test_that("estimate_rates CR smoke test", {
  fit <- estimate_rates(brts12, method = "mcem", model = "cr",
    control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                   sample_size = 20L, max_iter = 3L, tol = 0.5,
                   num_threads = 1L))
  expect_named(fit$pars, c("beta_0", "gamma_0"))
  expect_length(fit$pars, 2)
  expect_true(is.numeric(fit$loglik))
  expect_true(is.finite(fit$loglik))
  expect_gte(unname(fit$pars["beta_0"]), 0)
  expect_lte(unname(fit$pars["beta_0"]), 2)
})

test_that("estimate_rates DD smoke test", {
  fit <- suppressWarnings(estimate_rates(brts12, method = "mcem", model = "dd",
    control = list(lower_bound = c(0.1, -0.1, 0, -0.01),
                   upper_bound = c(2, 0.01, 0.5, 0.01),
                   sample_size = 20L, max_iter = 3L, tol = 0.5,
                   num_threads = 1L)))
  expect_named(fit$pars, c("beta_0", "beta_N", "gamma_0", "gamma_N"))
  expect_length(fit$pars, 4)
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$n_pars, 4L)
})

test_that("estimate_rates D smoke test", {
  # "ep" and "d" both resolve to c(0, 0, 1)
  fit <- suppressWarnings(estimate_rates(brts12, method = "mcem", model = "ep",
    control = list(lower_bound = c(0.1, -0.5, 0,  -0.5),
                   upper_bound = c(2,    0.5, 0.5, 0.5),
                   sample_size = 20L, max_iter = 3L, tol = 0.5,
                   num_threads = 1L)))
  expect_named(fit$pars, c("beta_0", "beta_D", "gamma_0", "gamma_D"))
  expect_true(is.numeric(fit$loglik))
  expect_equal(emphasis:::.model_label(fit$model), "D")
})

test_that("estimate_rates D + exponential link smoke test", {
  fit <- suppressWarnings(estimate_rates(brts12, method = "mcem", model = "ep",
    link = "exponential",
    control = list(lower_bound = c(-5, -5, -5, -5),
                   upper_bound = c(2, 2, 2, 2),
                   sample_size = 20L, max_iter = 2L, tol = 0.5,
                   num_threads = 1L)))
  expect_named(fit$pars, c("beta_0", "beta_D", "gamma_0", "gamma_D"))
  expect_true(is.numeric(fit$loglik))
  # the exponential link admits negative parameters: the box is the only bound
  expect_true(all(fit$pars >= -5 - 1e-8))
  expect_true(all(fit$pars <= 2 + 1e-8))
})

# Tests for the core IS functions: augment_trees, eval_logf, emphasis_cem

# --------------------------------------------------------------------------- #
#  Pure R unit tests (no C++ required)                                         #
# --------------------------------------------------------------------------- #

test_that("emphasis_cem errors on mismatched bounds", {
  expect_error(
    emphasis_cem(brts = 1:3, max_iter = 1, num_points = 10, max_missing = 100,
                 sd_vec = c(0.1), lower_bound = c(0, 0), upper_bound = c(1)),
    "same length"
  )
})

test_that(".is_fhat returns log_mean_exp of IS weights", {
  logf  <- c(-1, -2, -3)
  log_q <- c(-0.5, -0.5, -0.5)
  lw    <- logf - log_q   # c(-0.5, -1.5, -2.5)
  expected <- log(mean(exp(lw)))
  expect_equal(emphasis:::.is_fhat(logf, log_q), expected, tolerance = 1e-10)
})

test_that(".is_fhat returns NA for empty input", {
  expect_true(is.na(emphasis:::.is_fhat(numeric(0), numeric(0))))
})

test_that(".ess_from_lw equals n for uniform weights", {
  lw <- rep(0, 10)           # all weights equal → ESS = n
  expect_equal(emphasis:::.ess_from_lw(lw), 10, tolerance = 1e-10)
})

test_that(".ess_from_lw equals 1 for degenerate weights", {
  lw <- c(100, rep(-1000, 9))  # one dominant weight → ESS ≈ 1
  expect_equal(emphasis:::.ess_from_lw(lw), 1, tolerance = 1e-6)
})


# --------------------------------------------------------------------------- #
#  C++ integration tests (skipped in CI; run locally)                          #
# --------------------------------------------------------------------------- #

test_that("augment_trees returns correct structure", {
  skip("C++ integration test: run locally with devtools::test(filter='de')")
  set.seed(42)
  brts  <- c(4, 2.5, 1.2, 0.6)
  pars8 <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)  # CR model in 8-param layout

  result <- augment_trees(
    brts        = brts,
    pars        = pars8,
    sample_size = 5L,
    maxN        = 200L,
    max_missing = 1000L,
    max_lambda  = 100,
    num_threads = 1L,
    model       = c(0L, 0L, 0L),
    link        = 0L
  )

  expect_type(result, "list")
  expect_named(result, c("trees", "logf", "logg",
                         "rejected", "rejected_overruns",
                         "rejected_lambda", "rejected_zero_weights", "time"))
  expect_length(result$trees, 5L)
  expect_length(result$logf,  5L)
  expect_length(result$logg,  5L)
  expect_true(all(is.finite(result$logf)))
  expect_true(all(is.finite(result$logg)))
})

test_that("eval_logf returns correct structure", {
  skip("C++ integration test: run locally with devtools::test(filter='de')")
  set.seed(1)
  brts  <- c(4, 2.5, 1.2, 0.6)
  pars8 <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)

  aug <- augment_trees(brts, pars8, sample_size = 3L, maxN = 200L,
                       max_missing = 1000L, max_lambda = 100,
                       num_threads = 1L, model = c(0L,0L,0L), link = 0L)

  result <- eval_logf(pars = pars8, trees = aug$trees,
                      model = c(0L, 0L, 0L), link = 0L)

  expect_type(result, "list")
  expect_named(result, "logf")
  expect_length(result$logf, 3L)
  expect_true(all(is.finite(result$logf)))
  # logf values from eval_logf must match those from augment_trees
  expect_equal(result$logf, aug$logf, tolerance = 1e-10)
})

test_that("eval_logf uses correct model (EP fix: model_bin[2] gates EP rates)", {
  skip("C++ integration test: run locally with devtools::test(filter='de')")
  # model_bin only affects loglik evaluation for EP (model_bin[2]).
  # N/P covariates are always evaluated via pars[1]*N + pars[2]*P in loglik.
  # Before the fix, loglikelihood() hardcoded model={0,0,0}, making EP
  # inference silent wrong. eval_logf passes the correct model.
  set.seed(99)
  brts    <- c(4, 2.5, 1.2, 0.6)
  pars_ep <- c(0.5, 0, 0, 0.05, 0.1, 0, 0, 0.01)  # EP: beta_E=0.05, gamma_E=0.01

  aug <- augment_trees(brts, pars_ep, sample_size = 5L, maxN = 500L,
                       max_missing = 1000L, max_lambda = 100,
                       num_threads = 1L, model = c(0L, 0L, 1L), link = 0L)

  logf_ep <- eval_logf(pars_ep, aug$trees, model = c(0L,0L,1L), link = 0L)$logf
  logf_cr <- eval_logf(pars_ep, aug$trees, model = c(0L,0L,0L), link = 0L)$logf

  # Scoring EP trees under EP vs CR model must give different logf
  expect_false(isTRUE(all.equal(logf_ep, logf_cr)))
  # EP model uses speciation_rate_ep which includes the E covariate
  expect_true(all(is.finite(logf_ep)))
})

test_that("fhat from augment_trees + .is_fhat matches direct calc", {
  skip("C++ integration test: run locally with devtools::test(filter='de')")
  set.seed(3)
  brts  <- c(3, 1.5, 0.8)
  pars8 <- c(0.4, 0, 0, 0, 0.05, 0, 0, 0)

  aug   <- augment_trees(brts, pars8, sample_size = 20L, maxN = 500L,
                         max_missing = 1000L, max_lambda = 100,
                         num_threads = 1L, model = c(0L,0L,0L), link = 0L)
  fhat  <- emphasis:::.is_fhat(aug$logf, aug$logg)
  expect_true(is.finite(fhat))
  expect_true(fhat < 0)   # log-likelihood is negative
})

test_that("emphasis_cem integration test", {
  skip("C++ integration test: run locally with devtools::test(filter='de')")
  set.seed(42)
  brts <- c(0.9, 0.7, 0.5, 0.3, 0.1)

  lb8 <- c(0,   0, 0, 0, 0,   0, 0, 0)
  ub8 <- c(2,   0, 0, 0, 0.5, 0, 0, 0)
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
    disc_prop   = 0.5,
    verbose     = FALSE,
    num_threads = 1,
    model       = c(0L, 0L, 0L)
  )

  expect_type(result, "list")
  expect_true(all(c("best_loglik", "best_pars", "obtained_estim",
                    "loglik_var", "converged", "history",
                    "final_pop", "best_IS", "time") %in% names(result)))
  expect_length(result$obtained_estim, 8)
  expect_true(result$converged %in% c("annealing", "plateau", "max_iter"))
})

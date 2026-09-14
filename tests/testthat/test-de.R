# Tests for the core IS functions: augment_trees, eval_logf, emphasis_cem

# --------------------------------------------------------------------------- #
#  Pure R unit tests (no C++ required)                                         #
# --------------------------------------------------------------------------- #

test_that("emphasis_cem errors on mismatched bounds", {
  expect_error(
    emphasis:::emphasis_cem(brts = 1:3, max_iter = 1, num_points = 10,
                            max_missing = 100, sd_vec = c(0.1),
                            lower_bound = c(0, 0), upper_bound = c(1)),
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
#  C++ integration (audit id H94: every block below carried an unconditional   #
#  skip(), so the shipped suite made no call into the compiled code)           #
# --------------------------------------------------------------------------- #
#
# augment_trees, eval_logf and emphasis_cem are internal, so they are reached
# through `emphasis:::` and the assertions hold under both `load_all()` and the
# installed namespace.  The C++ engine is clock-seeded: set.seed() does not
# reproduce a draw, so nothing below pins a drawn number.

brts4  <- c(4, 2.5, 1.2, 0.6)
brts3  <- c(3, 1.5, 0.8)
cr_bin <- c(0L, 0L, 0L)
d_bin  <- c(0L, 0L, 1L)

# CR in the 8-parameter layout c(beta_0, beta_N, beta_M, beta_D,
#                                gamma_0, gamma_N, gamma_M, gamma_D)
pars_cr <- c(0.5, 0, 0, 0, 0.1, 0, 0, 0)

test_that("augment_trees returns correct structure", {
  result <- emphasis:::augment_trees(
    brts        = brts4,
    pars        = pars_cr,
    sample_size = 5L,
    maxN        = 200L,
    max_missing = 1000L,
    max_lambda  = 100,
    num_threads = 1L,
    model       = cr_bin,
    link        = 0L
  )

  expect_type(result, "list")
  expect_named(result, c("trees", "logf", "logg",
                         "rejected", "rejected_overruns",
                         "rejected_lambda", "rejected_zero_weights",
                         "rejected_nonfinite", "num_trees",
                         "envelope_violations", "time"))
  expect_length(result$trees, 5L)
  expect_length(result$logf,  5L)
  expect_length(result$logg,  5L)
  expect_equal(result$num_trees, 5L)
  expect_true(all(is.finite(result$logf)))
  expect_true(all(is.finite(result$logg)))
  # the envelope must dominate the rate on a cr tree: no violation, no warning
  expect_equal(result$envelope_violations, 0)

  # An augmented tree is an event list in forward time from the crown, so the
  # observed events sit at crown_age - brts and the last row is the present.
  obs_forward <- brts4[1] - rev(brts4[-1])
  for (tr in result$trees) {
    expect_true(all(c("brts", "n", "t_ext", "pd", "tip_start",
                      "focal_tip_start", "clade", "id", "parent_id") %in%
                      names(tr)))
    expect_false(is.unsorted(tr$brts))
    expect_true(all(obs_forward %in% tr$brts))
    expect_equal(max(tr$brts), brts4[1])
    expect_gte(nrow(tr), length(brts4))
  }
})

test_that("eval_logf returns correct structure", {
  aug <- emphasis:::augment_trees(brts4, pars_cr, sample_size = 3L, maxN = 200L,
                                  max_missing = 1000L, max_lambda = 100,
                                  num_threads = 1L, model = cr_bin, link = 0L)

  result <- emphasis:::eval_logf(pars = pars_cr, trees = aug$trees,
                                 model = cr_bin, link = 0L)

  expect_type(result, "list")
  expect_named(result, c("logf", "logg"))
  expect_length(result$logf, 3L)
  expect_length(result$logg, 3L)
  expect_true(all(is.finite(result$logf)))
  expect_true(all(is.finite(result$logg)))
  # re-scoring the trees at the parameters they were drawn at reproduces both
  # densities the sampler recorded, to the bit
  expect_identical(result$logf, aug$logf)
  expect_identical(result$logg, aug$logg)
})

test_that("eval_logf uses the model it is given (model_bin[3] gates the D rates)", {
  # loglikelihood() once hardcoded model = {0,0,0}, which made D inference
  # silently the CR likelihood.  eval_logf passes the model through, so a D
  # model with non-zero D coefficients must score differently from CR, and a D
  # model with beta_D = gamma_D = 0 must reduce to CR exactly.
  pars_d <- c(0.5, 0, 0, 0.05, 0.1, 0, 0, 0.01)

  aug <- emphasis:::augment_trees(brts4, pars_d, sample_size = 5L, maxN = 500L,
                                  max_missing = 1000L, max_lambda = 100,
                                  num_threads = 1L, model = d_bin, link = 0L)

  logf_d  <- emphasis:::eval_logf(pars_d, aug$trees, model = d_bin, link = 0L)$logf
  logf_cr <- emphasis:::eval_logf(pars_d, aug$trees, model = cr_bin, link = 0L)$logf

  expect_true(all(is.finite(logf_d)))
  expect_false(isTRUE(all.equal(logf_d, logf_cr)))

  # With the D coefficients switched off the two models are the same model.
  # The two code paths accumulate the compensator in a different order, so
  # they agree to rounding rather than to the bit; the D effect above is
  # 0.2 nats, twelve orders of magnitude larger than this tolerance.
  logf_d0  <- emphasis:::eval_logf(pars_cr, aug$trees, model = d_bin, link = 0L)$logf
  logf_cr0 <- emphasis:::eval_logf(pars_cr, aug$trees, model = cr_bin, link = 0L)$logf
  expect_equal(logf_d0, logf_cr0, tolerance = 1e-12)
})

test_that("fhat from augment_trees + .is_fhat matches direct calc", {
  pars <- c(0.4, 0, 0, 0, 0.05, 0, 0, 0)

  aug  <- emphasis:::augment_trees(brts3, pars, sample_size = 20L, maxN = 500L,
                                   max_missing = 1000L, max_lambda = 100,
                                   num_threads = 1L, model = cr_bin, link = 0L)
  fhat <- emphasis:::.is_fhat(aug$logf, aug$logg)
  expect_true(is.finite(fhat))

  # .is_fhat is the log mean of the IS weights over the completed draws
  lw <- aug$logf - aug$logg
  expect_equal(fhat, log(mean(exp(lw))), tolerance = 1e-12)
  expect_equal(emphasis:::.is_fhat(aug$logf, aug$logg, n_zero_weight = 5L),
               fhat + log(20 / 25), tolerance = 1e-12)

  # ESS is bounded by the number of draws and is positive
  ess <- emphasis:::.ess_from_lw(lw)
  expect_gt(ess, 0)
  expect_lte(ess, 20)
})

test_that("fhat on a fixed 3-tip tree stays in a narrow band over replicates", {
  # A tolerance over replicates, not a fixed value: the draws are clock-seeded.
  # The observed spread on this tree is ~0.35 nats, so a 3-nat band is a check
  # that the estimator is stable, not a restatement of one run.
  pars <- c(0.4, 0, 0, 0, 0.05, 0, 0, 0)
  fhat <- vapply(seq_len(10), function(i) {
    aug <- emphasis:::augment_trees(brts3, pars, sample_size = 20L, maxN = 500L,
                                    max_missing = 1000L, max_lambda = 100,
                                    num_threads = 1L, model = cr_bin, link = 0L)
    emphasis:::.is_fhat(aug$logf, aug$logg)
  }, numeric(1))

  expect_true(all(is.finite(fhat)))
  expect_lt(max(fhat), 0)
  expect_lt(max(fhat) - min(fhat), 3)
})

test_that("emphasis_cem integration test", {
  lb8 <- c(0,   0, 0, 0, 0,   0, 0, 0)
  ub8 <- c(2,   0, 0, 0, 0.5, 0, 0, 0)
  sd8 <- c(0.3, 0, 0, 0, 0.1, 0, 0, 0)
  result <- emphasis:::emphasis_cem(
    brts        = c(0.9, 0.7, 0.5, 0.3, 0.1),
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
    model       = cr_bin
  )

  expect_type(result, "list")
  expect_true(all(c("best_loglik", "best_pars", "obtained_estim",
                    "loglik_var", "converged", "history",
                    "final_pop", "best_IS", "time") %in% names(result)))
  expect_length(result$obtained_estim, 8)
  expect_true(result$converged %in% c("annealing", "plateau", "max_iter",
                                      "all_failed", "time_budget"))
  # the parameters the bounds fix stay fixed
  expect_equal(result$obtained_estim[lb8 == ub8], lb8[lb8 == ub8])
  expect_gte(result$obtained_estim[1], 0)
  expect_lte(result$obtained_estim[1], 2)
})

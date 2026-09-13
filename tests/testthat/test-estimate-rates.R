# Pins for estimate_rates() bookkeeping (audit ids H35, H21/H12 init,
# H22 maxN check, H31 print, H20 loglik lag): n_pars excludes fixed
# parameters, a symmetric box is a legal start for both samplers, the maxN
# check applies to the sampler that reads maxN, the reported loglik is the
# driver's own, and stop_reason / iterations / n_failed are surfaced at the
# top level and in print().  diagnose_mcem counts iterations, not trace rows.

# Fixed 20-tip branching times (crown age 5), decreasing order.
brts20 <- c(5, 4.4, 3.9, 3.6, 3.1, 2.8, 2.5, 2.2, 2.0, 1.7, 1.5, 1.3, 1.1,
            0.9, 0.75, 0.6, 0.45, 0.3, 0.15)

# ---------- H35: n_pars counts only free parameters ----------

test_that("n_pars excludes parameters fixed by lower_bound == upper_bound (H35)", {
  set.seed(5)
  lb <- c(0, -0.1, 0, 0)
  ub <- c(1, 0, 0.5, 0)   # gamma_N fixed at 0
  fit <- suppressWarnings(estimate_rates(
    brts20, model = "dd", link = "linear", method = "mcem",
    control = list(lower_bound = lb, upper_bound = ub,
                   sample_size = 20L, max_iter = 2L, max_time = 60,
                   num_threads = 1L)))
  expect_equal(fit$n_pars, 3L)
  expect_equal(fit$n_pars, sum(lb != ub))
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$AIC, -2 * fit$loglik + 2 * 3)
  expect_equal(unname(fit$pars["gamma_N"]), 0)
})

# ---------- H21 / H12: a symmetric box is a legal starting point ----------
#
# The default init is the midpoint of the bounds, so a symmetric box starts at
# lambda == mu.  The BDI transition probability there is the lambda = mu limit
# of the CR kernel (pinned in test-bdi-cr.R), so the midpoint needs no shift
# and no guard: these fits run.

test_that("symmetric box [0,1]^2 with default init runs on the BDI sampler (H12, H21)", {
  set.seed(121)
  fit <- estimate_rates(
    brts20, model = "cr", method = "mcem",
    control = list(lower_bound = c(0, 0), upper_bound = c(1, 1),
                   sampling = "bdi", sample_size = 20L, max_iter = 5L,
                   num_threads = 1L))
  expect_false(identical(fit$stop_reason, "e_step_failure"))
  expect_gte(fit$iterations, 1L)
  expect_true(is.finite(fit$loglik))
  # the midpoint (0.5, 0.5) is not returned unchanged
  expect_false(isTRUE(all.equal(unname(fit$pars), c(0.5, 0.5))))
})

test_that("symmetric box [0,2]^2 with default init runs on the BDI sampler (H21)", {
  set.seed(7)
  fit <- estimate_rates(
    brts20, model = "cr", method = "mcem",
    control = list(lower_bound = c(0, 0), upper_bound = c(2, 2),
                   sampling = "bdi", sample_size = 20L, max_iter = 5L,
                   num_threads = 1L))
  expect_false(identical(fit$stop_reason, "e_step_failure"))
  expect_gte(fit$iterations, 1L)
  expect_true(is.finite(fit$loglik))
})

test_that("symmetric box on the exponential link starts at lambda == mu and runs", {
  set.seed(8)
  fit <- estimate_rates(
    brts20, model = "cr", method = "mcem", link = "exponential",
    control = list(lower_bound = c(-2, -2), upper_bound = c(1, 1),
                   sampling = "bdi", sample_size = 20L, max_iter = 2L,
                   num_threads = 1L))
  expect_false(identical(fit$stop_reason, "e_step_failure"))
  expect_true(is.finite(fit$loglik))
})

test_that("mu fixed at the lambda midpoint runs on the BDI sampler (H12)", {
  set.seed(11)
  # lambda box [0, 1] with mu fixed at 0.5 = the lambda midpoint: the default
  # init has lambda == mu and mu cannot move.  This is a profile fit over
  # lambda, and it is a supported configuration.
  fit <- estimate_rates(
    brts20, model = "cr", method = "mcem",
    control = list(lower_bound = c(0, 0.5), upper_bound = c(1, 0.5),
                   sampling = "bdi", sample_size = 20L, max_iter = 2L,
                   num_threads = 1L))
  expect_s3_class(fit, "emphasis_fit")
  expect_true(is.finite(fit$loglik))
  expect_equal(unname(fit$pars["gamma_0"]), 0.5)
  expect_false(identical(fit$stop_reason, "e_step_failure"))
})

test_that("the thinning sampler also accepts a default init with lambda == mu", {
  set.seed(9)
  fit <- suppressWarnings(estimate_rates(
    brts20, model = "cr", method = "mcem",
    control = list(lower_bound = c(0, 0.5), upper_bound = c(1, 0.5),
                   sampling = "dynamic_fresh", sample_size = 20L,
                   max_iter = 2L, max_time = 60, num_threads = 1L)))
  expect_s3_class(fit, "emphasis_fit")
  expect_equal(unname(fit$pars["gamma_0"]), 0.5)
})

# ---------- H22 (part A): maxN must cover num_trees, under thinning ----------

test_that("estimate_rates errors when maxN < num_trees under thinning (H22 A)", {
  expect_error(
    estimate_rates(
      brts20, model = "cr", method = "mcem", init_pars = c(1, 0.3),
      control = list(sampling = "dynamic_fresh", num_trees = 300L, maxN = 200L,
                     max_iter = 2L, lower_bound = c(0, 0), upper_bound = c(4, 4),
                     num_threads = 1L)),
    "maxN")
})

test_that("the BDI sampler ignores maxN, including maxN < num_trees", {
  # .mcem_bdi has no maxN argument: it draws every tree to completion.
  expect_false("maxN" %in% names(formals(emphasis:::.mcem_bdi)))
  set.seed(12)
  small <- estimate_rates(
    brts20, model = "cr", method = "mcem", init_pars = c(1, 0.3),
    control = list(sampling = "bdi", num_trees = 30L, maxN = 5L,
                   max_iter = 2L, lower_bound = c(0, 0), upper_bound = c(4, 4),
                   num_threads = 1L))
  set.seed(12)
  big <- estimate_rates(
    brts20, model = "cr", method = "mcem", init_pars = c(1, 0.3),
    control = list(sampling = "bdi", num_trees = 30L, maxN = 5000L,
                   max_iter = 2L, lower_bound = c(0, 0), upper_bound = c(4, 4),
                   num_threads = 1L))
  expect_true(is.finite(small$loglik))
  expect_identical(small$pars, big$pars)
  expect_identical(small$loglik, big$loglik)
})

test_that("a non-finite maxN is treated like NULL", {
  set.seed(13)
  fit <- suppressWarnings(estimate_rates(
    brts20, model = "cr", method = "mcem", init_pars = c(1, 0.3),
    control = list(sampling = "dynamic_fresh", num_trees = 20L, maxN = NA,
                   max_iter = 2L, lower_bound = c(0, 0), upper_bound = c(4, 4),
                   max_time = 60, num_threads = 1L)))
  expect_true(is.finite(fit$loglik))
  expect_equal(fit$iterations, 2L)
  # the resolved cap is the NULL default, not NA
  expect_equal(unique(fit$details$mcem$maxN), 2000)
})

test_that("default mcem maxN is NULL so max(2000, 10 * num_trees) applies", {
  expect_null(estimate_rates_control("mcem")$maxN)
  # num_trees = 300 with the default maxN runs (maxN resolves to 3000 > 300)
  set.seed(22)
  fit <- suppressWarnings(estimate_rates(
    brts20, model = "cr", method = "mcem", init_pars = c(1, 0.3),
    control = list(sampling = "dynamic_fresh", num_trees = 300L, max_iter = 1L,
                   lower_bound = c(0, 0), upper_bound = c(4, 4),
                   max_time = 60, num_threads = 1L)))
  expect_equal(fit$iterations, 1L)
  expect_true(is.finite(fit$loglik))
})

# ---------- the stopping-rule default is shared by control and drivers ------

test_that("estimate_rates_control and both MCEM drivers carry the same tol", {
  expect_equal(estimate_rates_control("mcem")$tol, 1e-2)
  expect_equal(eval(formals(emphasis:::.mcem_bdi)$tol), 1e-2)
  expect_equal(eval(formals(emphasis:::.mcem_dynamic_fresh)$tol), 1e-2)
  expect_equal(estimate_rates_control("mcem")$patience, 3L)
})

# ---------- H20: the reported loglik is the driver's ----------
#
# The drivers evaluate a final E-step at the parameters they return and report
# its fhat.  The trace's rows pair theta_k with fhat(theta_{k-1}), so taking
# the trace tail reports the previous iterate's value whenever that final
# E-step failed.  The driver is mocked so the failure path is deterministic.

# One iteration row of a BDI-shaped trace, with fhat at the previous iterate.
lagged_trace <- function(fhat = -12.5) {
  cbind(
    as.data.frame(as.list(stats::setNames(
      c(1, 0, 0, 0, 0.3, 0, 0, 0), paste0("par", 1:8)))),
    data.frame(fhat = fhat, delta_max = 0.2, abs_step = 0.2, drift = 0.2,
               m_step = TRUE, m_moved = TRUE, rejected = 0L, n_nonfinite = 0L,
               num_trees = 10L, ESS = 8, time = 1))
}

run_mocked_mcem <- function(driver) {
  testthat::local_mocked_bindings(.mcem_bdi = driver, .package = "emphasis",
                                  .env = parent.frame())
  estimate_rates(
    brts20, model = "cr", method = "mcem", init_pars = c(1, 0.3),
    control = list(lower_bound = c(0, 0), upper_bound = c(4, 4),
                   sampling = "bdi", sample_size = 10L, max_iter = 1L,
                   num_threads = 1L))
}

test_that("a failed final E-step gives loglik NA, not the previous iterate's fhat (H20)", {
  fit <- run_mocked_mcem(function(...) list(
    mcem = lagged_trace(-12.5), pars = c(1.1, 0, 0, 0, 0.3, 0, 0, 0),
    iterations = 1L, stop_reason = "e_step_failure", loglik = NA_real_,
    loglik_var = NA_real_, final_IS = NULL, n_failed = 4L))
  expect_true(is.na(fit$loglik))
  expect_true(is.na(fit$AIC))
  expect_equal(fit$n_failed, 4L)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Did not converge")
})

test_that("the driver's loglik is preferred over the trace tail", {
  fit <- run_mocked_mcem(function(...) list(
    mcem = lagged_trace(-12.5), pars = c(1.1, 0, 0, 0, 0.3, 0, 0, 0),
    iterations = 1L, stop_reason = "max_iter", loglik = -7.25,
    loglik_var = NA_real_, final_IS = NULL, n_failed = 0L))
  expect_equal(fit$loglik, -7.25)
  expect_equal(fit$AIC, -2 * (-7.25) + 2 * 2)
})

test_that("a return without loglik falls back to the trace, not to loglik_var", {
  # `$` partial-matches, so a driver that reports only loglik_var must not
  # have its variance read as a log-likelihood.
  fit <- run_mocked_mcem(function(...) list(
    mcem = lagged_trace(-12.5), pars = c(1.1, 0, 0, 0, 0.3, 0, 0, 0),
    iterations = 1L, stop_reason = "max_iter", loglik_var = 0.5,
    final_IS = NULL, n_failed = 0L))
  expect_equal(fit$loglik, -12.5)
  expect_equal(fit$loglik_var, 0.5)
})

# ---------- H31: stop_reason / iterations at top level and in print ----------

test_that("all-fail MCEM run reports stop_reason, 0 iterations and prints them (H31)", {
  set.seed(2)
  # lambda fixed at 0: every augmented tree has zero weight, every E-step fails
  fit <- suppressWarnings(estimate_rates(
    brts20, model = "cr", link = "linear", method = "mcem",
    control = list(lower_bound = c(0, 0.1), upper_bound = c(0, 0.1),
                   sampling = "dynamic_fresh", sample_size = 10L, maxN = 100L,
                   max_iter = 20L, max_time = 60, num_threads = 1L)))
  expect_true(is.na(fit$loglik))
  expect_true(is.na(fit$AIC))
  expect_equal(fit$stop_reason, "e_step_failure")
  expect_equal(fit$iterations, 0L)
  expect_gt(fit$n_failed, 0L)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "e_step_failure")
  expect_match(out, "Did not converge")
  expect_match(out, "Iterations:\\s+0")
})

test_that("print shows iterations and stop reason for a completed run", {
  set.seed(3)
  fit <- estimate_rates(
    brts20, model = "cr", method = "mcem", init_pars = c(1, 0.3),
    control = list(lower_bound = c(0, 0), upper_bound = c(4, 4),
                   sampling = "bdi", sample_size = 20L, max_iter = 3L,
                   num_threads = 1L))
  expect_true(fit$stop_reason %in% c("converged", "max_iter", "time_budget"))
  # the final E-step is a trace row but not an iteration
  expect_equal(fit$iterations, sum(emphasis:::.mcem_iter_rows(fit$details$mcem)))
  expect_equal(nrow(fit$details$mcem), fit$iterations + 1L)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, paste0("Stop reason:\\s+", fit$stop_reason))
  expect_match(out, paste0("Iterations:\\s+", fit$iterations))
  expect_false(grepl("Did not converge", out))
})

# ---------- diagnose_mcem reports iterations, not trace rows ----------

test_that(".mcem_iter_rows reads both drivers' flags for the final E-step row", {
  bdi_trace  <- data.frame(fhat = c(-2, -1, -1), m_step = c(TRUE, TRUE, FALSE))
  thin_trace <- data.frame(fhat = c(-2, -1, -1),
                           final_estep = c(FALSE, FALSE, TRUE))
  expect_equal(emphasis:::.mcem_iter_rows(bdi_trace), c(TRUE, TRUE, FALSE))
  expect_equal(emphasis:::.mcem_iter_rows(thin_trace), c(TRUE, TRUE, FALSE))
  # a trace with neither flag is all iterations
  expect_equal(emphasis:::.mcem_iter_rows(data.frame(fhat = c(-2, -1))),
               c(TRUE, TRUE))
  expect_equal(emphasis:::.mcem_iter_rows(NULL), logical(0))
})

test_that("diagnose_mcem excludes the final E-step row on both samplers", {
  set.seed(14)
  fb <- estimate_rates(
    brts20, model = "cr", method = "mcem", init_pars = c(1, 0.3),
    control = list(lower_bound = c(0, 0), upper_bound = c(4, 4),
                   sampling = "bdi", sample_size = 20L, max_iter = 3L,
                   num_threads = 1L))
  db <- diagnose_mcem(fb, plot = FALSE)
  expect_equal(nrow(db$convergence), fb$iterations)
  expect_equal(nrow(fb$details$mcem), fb$iterations + 1L)
  expect_false(any(is.na(db$convergence$delta_max)))
  # the fhat of the final E-step is reported separately and is the fit's loglik
  expect_equal(db$final_fhat, fb$loglik)
  expect_match(paste(capture.output(print(db)), collapse = "\n"),
               sprintf("Iterations:\\s+%d", fb$iterations))

  set.seed(15)
  ft <- suppressWarnings(estimate_rates(
    brts20, model = "cr", method = "mcem", init_pars = c(1, 0.3),
    control = list(lower_bound = c(0, 0), upper_bound = c(4, 4),
                   sampling = "dynamic_fresh", sample_size = 20L, max_iter = 3L,
                   max_time = 60, num_threads = 1L)))
  dt <- diagnose_mcem(ft, plot = FALSE)
  expect_equal(nrow(dt$convergence), ft$iterations)
  expect_equal(nrow(ft$details$mcem), ft$iterations + 1L)
  expect_false(any(is.na(dt$convergence$delta_max)))
  expect_equal(dt$final_fhat, ft$loglik)
})

test_that("diagnose_mcem's rejected column is the trace's n_rejected (H34)", {
  # The thinning trace splits rejections into channels; n_rejected is their
  # sum, and `rejected` alone counts unhandled exceptions.  max_missing = 2
  # forces augmentation overruns, so the two differ.
  # ape::rcoal(20) with set.seed(4), rescaled to crown age 5 (H34.R).
  brts4 <- c(5.000000, 0.824025, 0.722480, 0.699953, 0.461824, 0.410440,
             0.364074, 0.331263, 0.313102, 0.177920, 0.143375, 0.125860,
             0.124385, 0.113191, 0.097307, 0.085446, 0.072051, 0.059174,
             0.002050)
  cr_bin <- c(0L, 0L, 0L)
  ex <- function(p) emphasis:::.expand_pars(p, cr_bin)
  set.seed(4)
  d <- emphasis:::.mcem_dynamic_fresh(
    brts4, ex(c(0.6, 0.5)), sample_size = 20L, maxN = 5000L, max_missing = 2,
    lower_bound = ex(c(0, 0)), upper_bound = ex(c(2, 1)), max_iter = 2L,
    xtol = 1e-3, tol = 1e-2, patience = 3L, num_threads = 1L, verbose = FALSE,
    model = cr_bin, link = 0L, max_time = 120)
  dg <- diagnose_mcem(d, plot = FALSE)
  iter_rows <- emphasis:::.mcem_iter_rows(d$mcem)
  expect_gt(sum(d$mcem$n_rejected), 0L)
  expect_equal(sum(d$mcem$rejected), 0L)          # a different channel
  expect_equal(dg$convergence$rejected, d$mcem$n_rejected[iter_rows])
  expect_equal(nrow(dg$convergence), d$iterations)
})

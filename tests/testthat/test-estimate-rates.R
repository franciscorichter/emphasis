# Pins for estimate_rates() bookkeeping (audit ids H35, H21/H12 init,
# H22 maxN check, H31 print): n_pars excludes fixed parameters, the default
# init avoids lambda == mu, maxN < num_trees is rejected, and stop_reason /
# iterations are surfaced at the top level and in print().

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

# ---------- H21 / H12: default init avoids lambda == mu ----------

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

test_that("symmetric box on the exponential link also starts with lambda != mu", {
  set.seed(8)
  fit <- estimate_rates(
    brts20, model = "cr", method = "mcem", link = "exponential",
    control = list(lower_bound = c(-2, -2), upper_bound = c(1, 1),
                   sampling = "bdi", sample_size = 20L, max_iter = 2L,
                   num_threads = 1L))
  expect_false(identical(fit$stop_reason, "e_step_failure"))
  expect_true(is.finite(fit$loglik))
})

test_that("default init errors when mu is fixed at the lambda midpoint under BDI", {
  expect_error(
    estimate_rates(
      brts20, model = "cr", method = "mcem",
      control = list(lower_bound = c(0, 0.5), upper_bound = c(1, 0.5),
                     sampling = "bdi", sample_size = 20L, max_iter = 2L,
                     num_threads = 1L)),
    "lambda == mu")
})

test_that("the thinning sampler still accepts a default init with lambda == mu", {
  set.seed(9)
  fit <- suppressWarnings(estimate_rates(
    brts20, model = "cr", method = "mcem",
    control = list(lower_bound = c(0, 0.5), upper_bound = c(1, 0.5),
                   sampling = "dynamic_fresh", sample_size = 20L,
                   max_iter = 2L, max_time = 60, num_threads = 1L)))
  expect_s3_class(fit, "emphasis_fit")
  expect_equal(unname(fit$pars["gamma_0"]), 0.5)
})

# ---------- H22 (part A): maxN must cover num_trees ----------

test_that("estimate_rates errors when maxN < num_trees (H22 A)", {
  expect_error(
    estimate_rates(
      brts20, model = "cr", method = "mcem", init_pars = c(1, 0.3),
      control = list(sampling = "dynamic_fresh", num_trees = 300L, maxN = 200L,
                     max_iter = 2L, lower_bound = c(0, 0), upper_bound = c(4, 4),
                     num_threads = 1L)),
    "maxN")
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
  expect_equal(fit$iterations, nrow(fit$details$mcem))
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, paste0("Stop reason:\\s+", fit$stop_reason))
  expect_match(out, paste0("Iterations:\\s+", fit$iterations))
  expect_false(grepl("Did not converge", out))
})

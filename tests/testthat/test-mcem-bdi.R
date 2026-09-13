# Tests for .mcem_bdi: stopping rule, patience accounting, trace columns,
# final E-step at the returned iterate, and recovery from E-step failures.
# Audit ids H11(d), H20, H21, H12 (recovery).
#
# Fixed branching times throughout; num_threads = 1, rho = 1, cond = NULL.

ns <- asNamespace("emphasis")

# 12-tip CR tree (TreeSim::sim.bd.taxa(12, 1, 0.5, 0.2), seed 21).
brts12 <- c(6.385824, 2.063997, 1.19255, 0.923743, 0.884126, 0.822976,
            0.725585, 0.718801, 0.539214, 0.077287, 0.068123)

# 9-tip DD tree of dev/audit/checks/H11.R (DDD::dd_sim(c(1.5, 0.4, 1.1/0.12),
# age = 6, ddmodel = 1), seed 11).  lambda(N) = 1.5 - 0.12 N is zero at N = 12.5.
brts_dd <- c(6, 4.848493, 4.401821, 3.108164, 3.073914, 2.835023, 1.838828,
             0.50463)

cr8  <- function(lam, mu) c(lam, 0, 0, 0, mu, 0, 0, 0)
lb8  <- cr8(0, 0)
ub8  <- cr8(3, 3)

# Relative change of the stopping rule, recomputed from two parameter vectors.
rel_change <- function(new, old, eps = 1e-2)
  max(abs(new - old) / pmax(abs(old), eps))

run_bdi <- function(pars, max_iter, tol = 1e-2, patience = 3L,
                    sample_size = 10L, lower_bound = lb8, upper_bound = ub8) {
  ns$.mcem_bdi(brts12, pars = pars, sample_size = sample_size,
               max_missing = 1e4, lower_bound = lower_bound,
               upper_bound = upper_bound, max_iter = max_iter, xtol = 1e-3,
               tol = tol, patience = patience, num_threads = 1L,
               model = c(0L, 0L, 0L), link = 0L, rho = 1)
}

# --------------------------------------------------------------------------- #
#  Stopping rule (H21)                                                         #
# --------------------------------------------------------------------------- #

test_that("delta_max is the relative change with floor eps = 1e-2 and is scaled by the parameter, not the box", {
  steps <- c(0.004, 0.004, 0.05, 0.004, 0.004, 0.004)   # relative moves of lambda
  k <- 0L
  testthat::local_mocked_bindings(
    m_cpp = function(e_step, init_pars, ...) {
      k <<- k + 1L
      est <- init_pars
      est[1] <- est[1] * (1 + steps[k])
      if (k == 1L) est[5] <- est[5] + 5e-4     # mu sits below eps
      list(estimates = est, nlopt = 4L, time = 0)
    },
    .package = "emphasis")

  set.seed(1)
  r <- run_bdi(cr8(0.5, 1e-3), max_iter = 20L, tol = 1e-2, patience = 3L,
               lower_bound = cr8(0, 0), upper_bound = cr8(300, 300))

  # Three consecutive iterations below tol are needed; the large third step
  # resets the streak, so the run stops at iteration 6.
  expect_identical(r$stop_reason, "converged")
  expect_equal(r$iterations, 6L)

  em <- r$mcem[r$mcem$m_step, ]
  expect_equal(nrow(em), 6L)
  # Iteration 1: lambda 0.5 -> 0.502 (0.4 %), mu 1e-3 -> 1.5e-3, measured
  # against eps = 1e-2 (5 %), so delta_max = 0.05 whatever the box.
  expect_equal(em$delta_max[1], 5e-4 / 1e-2)
  expect_equal(em$abs_step[1], 0.5 * 0.004)
  # Iteration 3: lambda step of 5 % dominates.
  expect_equal(em$delta_max[3], 0.05, tolerance = 1e-8)
  # drift is NA for the first `patience` iterations, then the relative
  # displacement over the last `patience` iterations.
  expect_true(all(is.na(em$drift[1:3])))
  pc <- grep("^par[0-9]+$", names(em), value = TRUE)
  P  <- as.matrix(em[, pc])
  for (j in 4:6) {
    expect_equal(em$delta_max[j], rel_change(P[j, ], P[j - 1, ]))
    expect_equal(em$abs_step[j], max(abs(P[j, ] - P[j - 1, ])))
    expect_equal(em$drift[j], rel_change(P[j, ], P[j - 3, ]))
  }
  # The recorded delta never references the box.
  expect_false(any(abs(em$delta_max - em$abs_step / 300) < 1e-12))
})

test_that("the default tol of .mcem_bdi is 1e-2", {
  expect_equal(eval(formals(ns$.mcem_bdi)$tol), 1e-2)
})

test_that("boxes [0,3], [0,30] and [0,300] give estimates within tolerance of each other (H21_widebox)", {
  fit_box <- function(ub) {
    set.seed(21)
    estimate_rates(brts12, model = "cr", method = "mcem",
                   init_pars = c(1.2, 0.9),
                   control = list(lower_bound = c(0, 0), upper_bound = c(ub, ub),
                                  sampling = "bdi", sample_size = 200L,
                                  num_threads = 1L, max_iter = 25L,
                                  tol = 1e-2, patience = 3L))
  }
  fits  <- lapply(c(3, 30, 300), fit_box)
  pars  <- t(sapply(fits, function(f) unname(f$pars)))
  iters <- sapply(fits, function(f) f$details$iterations)

  # The range-scaled rule stopped the wide boxes after exactly `patience`
  # iterations while still 0.5 away from the narrow-box fit.
  expect_true(all(iters > 3L))
  expect_lt(diff(range(pars[, 1])), 0.15)
  expect_lt(diff(range(pars[, 2])), 0.15)

  # delta_max is a function of the parameter columns alone, in every box.
  for (f in fits) {
    m  <- f$details$mcem
    em <- m[m$m_step, ]
    pc <- grep("^par[0-9]+$", names(em), value = TRUE)
    P  <- rbind(ns$.expand_pars(c(1.2, 0.9), c(0L, 0L, 0L)), as.matrix(em[, pc]))
    for (k in seq_len(nrow(em)))
      expect_equal(em$delta_max[k], rel_change(P[k + 1, ], P[k, ]))
  }
})

# --------------------------------------------------------------------------- #
#  Patience accounting and trace columns (H11 d)                              #
# --------------------------------------------------------------------------- #

test_that("an M-step that returns its start unchanged does not count toward patience", {
  testthat::local_mocked_bindings(
    m_cpp = function(e_step, init_pars, ...)
      list(estimates = init_pars, nlopt = 4L, time = 0),
    .package = "emphasis")

  set.seed(2)
  r <- run_bdi(cr8(0.5, 0.2), max_iter = 5L, tol = Inf, patience = 3L)

  expect_identical(r$stop_reason, "max_iter")
  expect_equal(r$iterations, 5L)
  em <- r$mcem[r$mcem$m_step, ]
  expect_true(all(em$delta_max == 0))
  expect_false(any(em$m_moved))
})

test_that("non-finite draws and rejections are recorded from the E-step, and non-finite draws are excluded from the M-step set", {
  real_aug <- ns$.augment_tree_bdi
  n_trees_seen <- integer(0)
  testthat::local_mocked_bindings(
    .augment_tree_bdi = function(tree, pars, ...) {
      e <- real_aug(tree, pars, ...)
      e$logf[1]    <- -Inf
      e$weights[1] <- -Inf
      e$n_rejected <- 7L
      e
    },
    m_cpp = function(e_step, init_pars, ...) {
      n_trees_seen <<- c(n_trees_seen, length(e_step$trees))
      est <- init_pars
      est[1] <- est[1] * 1.001
      list(estimates = est, nlopt = 4L, time = 0)
    },
    .package = "emphasis")

  set.seed(3)
  r <- run_bdi(cr8(0.5, 0.2), max_iter = 2L, tol = 0, sample_size = 10L)

  m <- r$mcem
  expect_type(m$rejected, "integer")
  expect_type(m$n_nonfinite, "integer")
  expect_true(all(m$rejected == 7L))
  expect_true(all(m$n_nonfinite == 1L))
  expect_true(all(m$num_trees == 10L))
  expect_equal(n_trees_seen, c(9L, 9L))
  expect_equal(r$final_IS$n_rejected, 7L)
  expect_equal(r$final_IS$rejected_zero_weights, 1L)
})

test_that("the H11 part D configuration does not report convergence through zero-delta iterations", {
  lb <- c(0.01, -1, 0.001, 0); ub <- c(5, 0, 2, 0)
  set.seed(11)
  fit <- estimate_rates(brts_dd, method = "mcem", model = "dd",
                        init_pars = c(1.5, -0.12, 0.4, 0),
                        control = list(lower_bound = lb, upper_bound = ub,
                                       sample_size = 200L, max_iter = 4L,
                                       max_missing = 30L, num_threads = 1L,
                                       sampling = "bdi"))
  m  <- fit$details$mcem
  em <- m[m$m_step, ]
  expect_true(all(c("delta_max", "abs_step", "drift", "m_moved",
                    "n_nonfinite", "rejected") %in% names(m)))
  expect_false(identical(fit$details$stop_reason, "converged") &&
               any(utils::tail(em$delta_max, 3L) == 0))
  if (identical(fit$details$stop_reason, "converged")) {
    last <- utils::tail(em, 3L)
    expect_true(all(last$m_moved))
    expect_true(all(last$delta_max > 0))
  }
  expect_true(any(em$m_moved))
})

# --------------------------------------------------------------------------- #
#  Final E-step at the returned iterate (H20)                                 #
# --------------------------------------------------------------------------- #

test_that("fit$loglik is fhat at fit$pars, from an E-step recorded as the last trace row", {
  set.seed(20)
  fit <- estimate_rates(brts12, model = "cr", method = "mcem",
                        init_pars = c(1.2, 0.9),
                        control = list(lower_bound = c(0, 0), upper_bound = c(3, 3),
                                       sampling = "bdi", sample_size = 50L,
                                       num_threads = 1L, max_iter = 2L,
                                       tol = 1e-2))
  m  <- fit$details$mcem
  th <- as.numeric(fit$details$pars)

  expect_equal(fit$details$iterations, 2L)
  expect_equal(nrow(m), 3L)
  expect_equal(m$m_step, c(TRUE, TRUE, FALSE))
  expect_true(is.na(m$delta_max[3]))
  pc <- grep("^par[0-9]+$", names(m), value = TRUE)
  expect_equal(unname(unlist(m[3, pc])), th)

  # CR BDI is zero-variance IS: fhat(theta) does not depend on the draw.
  f_K <- ns$.augment_tree_bdi(brts12, th, model_bin = c(0L, 0L, 0L),
                              sample_size = 50L, max_missing = 1e4,
                              link = 0L, rho = 1)$fhat
  expect_equal(fit$loglik, f_K, tolerance = 1e-8)
  expect_equal(fit$details$loglik, fit$loglik)
  expect_equal(fit$details$final_IS$fhat, fit$loglik)
  expect_equal(fit$details$final_IS$ESS, 50)
  expect_equal(fit$AIC, -2 * f_K + 2 * fit$n_pars)

  # Two M-steps from a far start move the likelihood by more than the IS
  # noise, so the lagged value fhat(theta_{K-1}) is not the reported one.
  expect_gt(abs(m$fhat[2] - fit$loglik), 1e-3)
})

# --------------------------------------------------------------------------- #
#  Recovery from E-step failure (H12)                                          #
# --------------------------------------------------------------------------- #

test_that("an E-step failure restarts from the last successful iterate, not from the box centre", {
  real_aug <- ns$.augment_tree_bdi
  calls <- list()
  testthat::local_mocked_bindings(
    .augment_tree_bdi = function(tree, pars, ...) {
      calls[[length(calls) + 1L]] <<- pars
      if (length(calls) == 2L) stop("simulated E-step failure")
      real_aug(tree, pars, ...)
    },
    .package = "emphasis")

  set.seed(4)
  p0 <- cr8(0.5, 0.2)
  r  <- run_bdi(p0, max_iter = 3L, tol = 0, patience = 3L,
                lower_bound = cr8(0, 0), upper_bound = cr8(4, 4))

  # Loop pass 1: call 1 at p0 (ok) -> theta_1.  Pass 2: call 2 at theta_1
  # fails -> restart at p0.  Pass 3: call 3 at p0 (ok) -> theta_1'.
  # Call 4: final E-step at theta_1'.
  expect_length(calls, 4L)
  expect_identical(calls[[3]], calls[[1]])
  expect_false(isTRUE(all.equal(calls[[3]], 0.8 * calls[[2]] + 0.2 * cr8(2, 2))))
  expect_equal(r$iterations, 2L)
  expect_equal(r$n_failed, 1L)
  expect_identical(r$stop_reason, "max_iter")
  expect_equal(nrow(r$mcem), 3L)
  expect_equal(calls[[4]], r$pars)
})

test_that("eight consecutive E-step failures return the init unchanged with iterations 0 and loglik NA", {
  n_calls <- 0L
  testthat::local_mocked_bindings(
    .augment_tree_bdi = function(tree, pars, ...) {
      n_calls <<- n_calls + 1L
      stop("simulated E-step failure")
    },
    .package = "emphasis")

  p0 <- cr8(0.7, 0.3)
  r  <- run_bdi(p0, max_iter = 50L, tol = 0)

  expect_identical(r$stop_reason, "e_step_failure")
  expect_identical(r$pars, p0)
  expect_equal(r$iterations, 0L)
  expect_equal(r$n_failed, 8L)
  expect_null(r$mcem)
  expect_true(is.na(r$loglik))
  expect_true(is.na(r$loglik_var))
  expect_null(r$final_IS)
  expect_equal(n_calls, 9L)   # 8 attempts plus the final E-step
})

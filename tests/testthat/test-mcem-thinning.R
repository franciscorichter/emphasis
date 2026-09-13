# Pins the thinning MCEM driver (.mcem_dynamic_fresh) after the wave-1 fix for
# H20, H21, H22, H30, H31, H34:
#   - stopping rule max_j |dtheta_j| / max(|theta_{k-1,j}|, 1e-2) < tol for
#     `patience` consecutive iterations, independent of the bound box
#   - a final E-step at theta_K so fhat / final_IS describe the returned pars
#   - maxN ratchets on E-step failure and is never reset after a success
#   - a failure restarts from the last estimate; the box-centre perturbation
#     starts only from the second consecutive failure; the streak resets
#   - eight consecutive failures return the last estimate (init) flagged with
#     iterations = 0 and n_failed = 8, not a mixture with the box centre
#   - separate rejection columns whose sum is final_IS$n_rejected
#
# Branching-time vectors are fixed (from the audit scripts H22.R, H30.R, H34.R).
# num_threads = 1 throughout.

# ape::rphylo(20, 1, 0.3) with set.seed(22)  (H22.R)
brts22 <- c(2.984340, 2.456040, 2.041011, 1.919104, 1.703097, 1.313093,
            0.939284, 0.715800, 0.553944, 0.455806, 0.427175, 0.308270,
            0.174959, 0.141573, 0.112195, 0.069056, 0.058366, 0.044335,
            0.035010)
# ape::rcoal(20) with set.seed(1), rescaled to crown age 5  (H30.R)
brts1 <- c(5.000000, 3.822147, 3.059192, 1.966150, 1.604249, 1.358484,
           0.622046, 0.467531, 0.393533, 0.285495, 0.276148, 0.225482,
           0.201295, 0.154060, 0.057677, 0.044974, 0.041380, 0.038051,
           0.013895)
# ape::rcoal(20) with set.seed(4), rescaled to crown age 5  (H34.R)
brts4 <- c(5.000000, 0.824025, 0.722480, 0.699953, 0.461824, 0.410440,
           0.364074, 0.331263, 0.313102, 0.177920, 0.143375, 0.125860,
           0.124385, 0.113191, 0.097307, 0.085446, 0.072051, 0.059174,
           0.002050)

cr_bin <- c(0L, 0L, 0L)
ex <- function(p) emphasis:::.expand_pars(p, cr_bin)

rej_cols <- c("rejected", "rejected_overruns", "rejected_lambda")

test_that("sample_size > maxN: maxN ratchets once and the fit stays away from the box centre (H22 B)", {
  lb <- c(0, 0); ub <- c(4, 4); centre <- (lb + ub) / 2
  # estimate_rates() rejects maxN < num_trees up front (H22 part A), so the
  # ratchet is exercised on the driver directly.
  set.seed(2222)
  d <- suppressWarnings(emphasis:::.mcem_dynamic_fresh(
    brts22, ex(c(1, 0.3)), sample_size = 300L, maxN = 200L, max_missing = 1e4,
    lower_bound = ex(lb), upper_bound = ex(ub), max_iter = 40L, xtol = 1e-3,
    tol = 1e-2, patience = 3L, num_threads = 1L, verbose = FALSE,
    model = cr_bin, link = 0L, max_time = 300))
  fit <- list(pars = emphasis:::.contract_pars(d$pars, cr_bin))

  # The first E-step fails structurally (200 attempts < 300 trees); the
  # doubled cap then holds for the rest of the run.
  expect_equal(d$n_failed, 1L)
  expect_true(all(d$mcem$maxN >= 400L))
  expect_equal(d$maxN, 400L)
  expect_gt(d$iterations, 1L)
  expect_equal(nrow(d$mcem), d$iterations + 1L)   # + final E-step row

  # No pull toward the box centre: the clean fit on this tree is near
  # (1.0, 0.16); the old alternating fit settled near (1.9, 1.4).
  expect_gt(max(abs(as.numeric(fit$pars) - centre)), 1.0)
  expect_lt(as.numeric(fit$pars)[2], 1)

  # Stopping metric is the relative step, recomputed from the trace rows
  # (rows 2..K measure against the previous row; row 1 against init).
  m <- d$mcem[!d$mcem$final_estep, ]
  P <- as.matrix(m[, grep("^par[0-9]+$", names(m))])
  for (k in seq_len(nrow(P))[-1]) {
    ref <- max(abs(P[k, ] - P[k - 1, ]) / pmax(abs(P[k - 1, ]), 1e-2))
    expect_equal(m$delta_max[k], ref)
    expect_equal(m$abs_step[k], max(abs(P[k, ] - P[k - 1, ])))
  }
  p0 <- ex(c(1, 0.3))
  expect_equal(m$delta_max[1], max(abs(P[1, ] - p0) / pmax(abs(p0), 1e-2)))
})

test_that("final E-step: fit$loglik is fhat at fit$pars and final_IS matches the last row (H20)", {
  lb8 <- ex(c(0, 0)); ub8 <- ex(c(2, 1))
  set.seed(20)
  res <- emphasis:::.mcem_dynamic_fresh(
    brts1, ex(c(0.5, 0.1)), sample_size = 30L, maxN = 3000L, max_missing = 1e4,
    lower_bound = lb8, upper_bound = ub8, max_iter = 3L, xtol = 1e-3,
    tol = 1e-2, patience = 3L, num_threads = 1L, verbose = FALSE,
    model = cr_bin, link = 0L, max_time = 120)
  m <- res$mcem
  expect_equal(res$stop_reason, "max_iter")
  expect_equal(res$iterations, 3L)
  expect_true(res$final_estep)
  expect_equal(nrow(m), 4L)
  expect_equal(m$final_estep, c(FALSE, FALSE, FALSE, TRUE))

  last <- m[nrow(m), ]
  # The final row is evaluated at the returned pars ...
  expect_equal(as.numeric(last[, paste0("par", 1:8)]), as.numeric(res$pars))
  # ... and its fhat is what .run_mcem reports as loglik, and what final_IS holds
  expect_equal(res$final_IS$fhat, last$fhat)
  expect_equal(length(res$final_IS$logf), 30L)
  expect_true(is.na(last$delta_max))

  # Through estimate_rates: loglik is the final-row fhat, not the lagged one
  set.seed(21)
  fit <- estimate_rates(
    brts1, model = "cr", method = "mcem", init_pars = c(0.5, 0.1),
    control = list(sampling = "dynamic_fresh", num_trees = 30L, maxN = 3000L,
                   max_iter = 2L, tol = 1e-2, patience = 3L,
                   lower_bound = c(0, 0), upper_bound = c(2, 1), num_threads = 1L))
  mm <- fit$details$mcem
  expect_equal(fit$loglik, mm$fhat[nrow(mm)])
  expect_true(mm$final_estep[nrow(mm)])
  expect_equal(fit$details$iterations, 2L)
})

test_that("eight consecutive failures return the init flagged with iterations 0 (H22 C, H31)", {
  lb8 <- ex(c(0, 0)); ub8 <- ex(c(4, 4)); init8 <- ex(c(1, 0.3))
  boom <- function(p) stop("injected failure")
  set.seed(3)
  res <- suppressWarnings(emphasis:::.mcem_dynamic_fresh(
    brts22, init8, sample_size = 20L, maxN = 200L, max_missing = 1e4,
    lower_bound = lb8, upper_bound = ub8, max_iter = 20L, xtol = 1e-3,
    tol = 1e-2, patience = 3L, num_threads = 1L, verbose = FALSE,
    conditional = boom, model = cr_bin, link = 0L, max_time = 120))
  expect_equal(res$stop_reason, "e_step_failure")
  expect_equal(res$iterations, 0L)
  expect_equal(res$n_failed, 8L)
  expect_false(res$final_estep)
  expect_null(res$mcem)
  expect_null(res$final_IS)
  expect_true(is.na(res$loglik_var))
  # The returned pars is init, not 0.8^8 init + (1 - 0.8^8) centre
  expect_equal(as.numeric(res$pars), as.numeric(init8))
  centre8 <- (lb8 + ub8) / 2
  mixture <- 0.8^8 * init8 + (1 - 0.8^8) * centre8
  expect_gt(max(abs(as.numeric(res$pars) - mixture)), 0.1)
  # maxN escalated on every failure up to the cap
  expect_equal(res$maxN, 50000L)

  # Same through estimate_rates with a structural cause (lambda fixed at 0:
  # every augmented tree has zero weight), as in H31.R
  w <- character(0)
  fit <- withCallingHandlers(
    estimate_rates(brts1, model = "cr", link = "linear", method = "mcem",
                   control = list(lower_bound = c(0, 0.1), upper_bound = c(0, 0.1),
                                  sampling = "dynamic_fresh", sample_size = 10L,
                                  maxN = 100L, max_iter = 20L, max_time = 60,
                                  num_threads = 1L)),
    warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
  expect_equal(fit$details$stop_reason, "e_step_failure")
  expect_identical(fit$details$iterations, 0L)
  expect_equal(fit$details$n_failed, 8L)
  expect_true(is.na(fit$loglik))
  expect_equal(as.numeric(fit$pars), c(0, 0.1))
  expect_true(length(w) >= 1L)
})

test_that("an interleaved E-step failure resets the streak and restarts from the last estimate (H30)", {
  lb8 <- ex(c(0, 0)); ub8 <- ex(c(2, 1)); ip8 <- ex(c(0.5, 0.1))
  flag <- new.env(); flag$fail <- FALSE
  cond <- function(p) {
    if (flag$fail) { flag$fail <- FALSE; stop("injected E-step failure") }
    0
  }
  msgs <- character(0)
  set.seed(1)
  res <- withCallingHandlers(
    emphasis:::.mcem_dynamic_fresh(
      brts1, ip8, sample_size = 20L, maxN = 2000L, max_missing = 1e4,
      lower_bound = lb8, upper_bound = ub8, max_iter = 20L, xtol = 1e-3,
      tol = 1e6,            # every successful iteration counts as stable
      patience = 3L, num_threads = 1L, verbose = TRUE, conditional = cond,
      model = cr_bin, link = 0L, max_time = 120),
    message = function(m) {
      txt <- conditionMessage(m); msgs <<- c(msgs, txt)
      if (grepl("^Iteration 2:", txt)) flag$fail <- TRUE
      invokeRestart("muffleMessage")
    })
  m <- res$mcem[!res$mcem$final_estep, ]

  # Iterations 1, 2 succeed, 3 fails, then three more successes are needed:
  # the streak does not span the failure.
  expect_equal(res$stop_reason, "converged")
  expect_equal(res$n_failed, 1L)
  expect_equal(res$iterations, 5L)
  expect_equal(nrow(m), 5L)
  expect_true(any(grepl("^Iteration 3: E-step failed", msgs)))
  expect_true(any(grepl("restarting from the last estimate", msgs)))
  expect_false(any(grepl("perturbed toward center", msgs)))

  # The retry after the failure ran from the last estimate: the third
  # recorded step is measured against the second recorded iterate, not
  # against a box-centre mixture.
  P <- as.matrix(m[, paste0("par", 1:8)])
  ref <- max(abs(P[3, ] - P[2, ]) / pmax(abs(P[2, ]), 1e-2))
  expect_equal(m$delta_max[3], ref)
  centre8 <- (lb8 + ub8) / 2
  pert <- pmin(pmax(0.8 * P[2, ] + 0.2 * centre8, lb8), ub8)
  expect_false(isTRUE(all.equal(m$delta_max[3],
                                max(abs(P[3, ] - pert) / pmax(abs(pert), 1e-2)))))
  # maxN was raised by the failure and kept afterwards
  expect_equal(m$maxN, c(2000L, 2000L, 4000L, 4000L, 4000L))
  expect_equal(res$maxN, 4000L)
})

test_that("second consecutive failure perturbs toward the centre; pars stays the last estimate", {
  lb8 <- ex(c(0, 0)); ub8 <- ex(c(2, 1)); ip8 <- ex(c(0.5, 0.1))
  flag <- new.env(); flag$left <- 0L
  cond <- function(p) {
    if (flag$left > 0L) { flag$left <- flag$left - 1L; stop("injected E-step failure") }
    0
  }
  msgs <- character(0)
  set.seed(2)
  res <- withCallingHandlers(
    emphasis:::.mcem_dynamic_fresh(
      brts1, ip8, sample_size = 20L, maxN = 2000L, max_missing = 1e4,
      lower_bound = lb8, upper_bound = ub8, max_iter = 20L, xtol = 1e-3,
      tol = 1e6, patience = 2L, num_threads = 1L, verbose = TRUE,
      conditional = cond, model = cr_bin, link = 0L, max_time = 120),
    message = function(m) {
      txt <- conditionMessage(m); msgs <<- c(msgs, txt)
      if (grepl("^Iteration 1:", txt)) flag$left <- 2L
      invokeRestart("muffleMessage")
    })
  expect_equal(res$n_failed, 2L)
  expect_equal(res$iterations, 3L)      # 1 success, 2 failures, 2 successes
  expect_true(any(grepl("^Iteration 2: E-step failed \\(1 consecutive\\).*restarting", msgs)))
  expect_true(any(grepl("^Iteration 3: E-step failed \\(2 consecutive\\).*perturbed", msgs)))
  m <- res$mcem[!res$mcem$final_estep, ]
  P <- as.matrix(m[, paste0("par", 1:8)])
  # After the second failure the E-step ran from 0.8 * theta_1 + 0.2 * centre,
  # and the recorded step is measured against that point.
  centre8 <- (lb8 + ub8) / 2
  pert <- pmin(pmax(0.8 * P[1, ] + 0.2 * centre8, lb8), ub8)
  expect_equal(m$delta_max[2], max(abs(P[2, ] - pert) / pmax(abs(pert), 1e-2)))
  expect_equal(m$maxN, c(2000L, 8000L, 8000L))
})

test_that("rejection columns are separate and sum to final_IS$n_rejected (H34)", {
  lb8 <- ex(c(0, 0)); ub8 <- ex(c(2, 1)); ip8 <- ex(c(0.6, 0.5))
  set.seed(4)
  res <- emphasis:::.mcem_dynamic_fresh(
    brts4, ip8, sample_size = 20L, maxN = 5000L,
    max_missing = 2,          # tiny -> many augmentation overruns
    lower_bound = lb8, upper_bound = ub8, max_iter = 2L, xtol = 1e-3,
    tol = 1e-2, patience = 3L, num_threads = 1L, verbose = FALSE,
    model = cr_bin, link = 0L, max_time = 120)
  m <- res$mcem
  expect_true(all(c(rej_cols, "rejected_zero_weights", "n_rejected") %in% names(m)))
  expect_equal(m$n_rejected, rowSums(m[, rej_cols]))
  last <- m[nrow(m), ]
  expect_true(last$final_estep)
  expect_equal(sum(last[, rej_cols]), res$final_IS$n_rejected)
  expect_equal(last$n_rejected, res$final_IS$n_rejected)
  expect_equal(last$rejected_zero_weights, res$final_IS$rejected_zero_weights)
  # With max_missing = 2 overruns occur; they are no longer hidden behind
  # the unhandled-exception counter.
  expect_gt(sum(m$rejected_overruns), 0L)
  expect_gt(res$final_IS$n_rejected, 0L)
})

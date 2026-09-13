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


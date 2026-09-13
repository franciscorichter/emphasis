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


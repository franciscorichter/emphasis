# Tests for simulation functions

# ---------- Input validation for simulate_tree ----------

test_that("simulate_tree errors on invalid pars", {
  expect_error(
    simulate_tree(pars = "x", max_t = 5),
    "non-empty numeric"
  )
  expect_error(
    simulate_tree(pars = numeric(0), max_t = 5),
    "non-empty numeric"
  )
})

test_that("simulate_tree errors on invalid max_t", {
  expect_error(
    simulate_tree(pars = c(0.5, 0.1), max_t = -1),
    "positive number"
  )
  expect_error(
    simulate_tree(pars = c(0.5, 0.1), max_t = 0),
    "positive number"
  )
})

test_that("simulate_tree errors on invalid max_lin", {
  expect_error(
    simulate_tree(pars = c(0.5, 0.1), max_t = 5, max_lin = 0),
    "positive number"
  )
})

test_that("simulate_tree accepts max_tries = 0 (no retries)", {
  set.seed(1)
  # max_tries = 0 is valid: extinct runs are not retried (used for GAM training)
  result <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, max_tries = 0,
                          useDDD = FALSE)
  expect_true(result$status %in% c("done", "extinct", "too_large"))
})

test_that("simulate_tree errors on wrong pars length for model", {
  # "dd" = c(1,0,0) needs 4 pars, supplying 2 should error
  expect_error(
    simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "dd"),
    "requires 4 parameters"
  )
})

test_that("simulate_tree errors on negative max_t", {
  expect_error(
    simulate_tree(pars = c(0.5, 0.1), max_t = -1),
    "positive number"
  )
})

test_that(".resolve_model handles shortcuts and binary vectors", {
  # five slots {N, M, D, ED, K}; a shorter vector is an earlier layout, padded
  expect_equal(emphasis:::.resolve_model("cr"), c(0L, 0L, 0L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model("dd"), c(1L, 0L, 0L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model("d"), c(0L, 0L, 1L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model("ep"), c(0L, 0L, 1L, 0L, 0L))
  expect_equal(emphasis:::.resolve_model(c(1, 0, 1)), c(1L, 0L, 1L, 0L, 0L))
  expect_error(
    emphasis:::.resolve_model(c(1, 2, 0)),
    "binary integer vector of length 3 to 5"
  )
})

test_that(".expand_pars errors on wrong length", {
  expect_error(
    emphasis:::.expand_pars(c(1, 2), c(1L, 0L, 0L)),
    "requires 4 parameters"
  )
})

# ---------- C++ simulation (audit id H94: these blocks were skipped) ---------
#
# The C++ engine is clock-seeded, so set.seed() does not reproduce a draw.
# Every assertion below is therefore either structural or a tolerance over
# replicates, never a fixed value.

test_that("simulate_tree cr model returns valid output", {
  result <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")

  expect_type(result, "list")
  expect_named(result, c("tes", "tas", "L", "status", "survival_prob", "attempts"))
  expect_true(result$status %in% c("done", "extinct", "too_large"))
  if (result$status == "done") {
    # Ltable in DDD format: birth time, parent, id, death time (-1 = extant)
    expect_equal(ncol(result$L), 4L)
    expect_equal(result$L[1:2, 1], c(5, 5))       # both crown lineages at max_t
    expect_true(all(result$L[, 1] <= 5))
    expect_gte(sum(result$L[, 4] == -1), 2L)      # both crown clades survive
    expect_s3_class(result$tes, "phylo")
    expect_equal(ape::Ntip(result$tes), sum(result$L[, 4] == -1))
  } else {
    expect_null(result$L)
  }
})

test_that("simulate_tree dd model returns valid output", {
  result <- simulate_tree(
    pars = c(0.8, -0.02, 0.2, -0.01),
    max_t = 5, model = "dd"
  )

  expect_type(result, "list")
  expect_true(result$status %in% c("done", "extinct", "too_large"))
  if (result$status == "done") {
    expect_equal(ncol(result$L), 4L)
    # beta_N = -0.02 puts the zero of lambda at N = 40, so a dd tree of these
    # parameters cannot run away
    expect_lt(nrow(result$L), 1e5)
  }
})

test_that("simulate_tree returns survival_prob", {
  result <- simulate_tree(pars = c(0.8, 0.1), max_t = 5, model = "cr")
  expect_true(is.numeric(result$survival_prob))
  expect_true(result$survival_prob >= 0 && result$survival_prob <= 1)
  # survival_prob is 1/attempts on success and 0 on failure; max_tries = 1
  # allows two attempts
  expect_true(result$survival_prob %in% c(0, 0.5, 1))
  expect_equal(result$survival_prob > 0, result$status == "done")
})

# ---------- the forward simulator against the birth-death moments -----------
#
# The simulator starts from the two crown lineages and stops with status
# "extinct" as soon as either crown clade is empty (general_tree.hpp:290,
# N1 < 1 || N2 < 1), so "done" is the event that both clades are alive at
# max_t.  With u the probability that one clade is extinct at T,
#
#   P(done) = (1 - u)^2 ,   E[N_T | done] = 2 e^{(lambda - mu) T} / (1 - u) ,
#
# the second because E[N_T] = e^{rT} per clade with N = 0 on extinction, and
# the two clades are independent.  max_tries = 0 means one attempt, so the
# draws are unconditional and these are the moments to compare against.

bd_extinction_prob <- function(lambda, mu, tt) {
  ert <- exp((lambda - mu) * tt)
  mu * (ert - 1) / (lambda * ert - mu)
}

test_that("forward cr simulation matches the birth-death moments", {
  skip_on_cran()     # a 300-replicate Monte Carlo gate with 4-sigma tolerances
  lambda <- 0.5
  mu     <- 0.1
  tt     <- 5
  reps   <- 300L

  n_extant <- vapply(seq_len(reps), function(i) {
    x <- simulate_tree(pars = c(lambda, mu), max_t = tt, model = "cr",
                       max_tries = 0, useDDD = FALSE)
    if (is.null(x$L)) 0L else sum(x$L[, 4] == -1)
  }, integer(1))

  u        <- bd_extinction_prob(lambda, mu, tt)
  p_done   <- (1 - u)^2
  mean_cnd <- 2 * exp((lambda - mu) * tt) / (1 - u)

  done   <- n_extant > 0L
  p_hat  <- mean(done)
  se_p   <- sqrt(p_done * (1 - p_done) / reps)
  expect_lt(abs(p_hat - p_done), 4 * se_p)

  kept   <- n_extant[done]
  expect_gt(length(kept), 100L)
  se_n   <- stats::sd(kept) / sqrt(length(kept))
  expect_lt(abs(mean(kept) - mean_cnd), 4 * se_n)

  # every surviving draw carries at least the two crown lineages
  expect_true(all(kept >= 2L))
})

test_that("an overflowing clade is not retried, an extinct one is", {
  set.seed(5)
  # lambda far above what max_lin allows: the first draw overflows and the
  # retry budget is not spent on it
  big <- simulate_tree(pars = c(5, 0.0), max_t = 10, model = "cr", link = "linear",
                       max_lin = 20L, max_tries = 20L)
  expect_identical(big$status, "too_large")
  expect_identical(big$attempts, 1L)
  # extinction dominates: retries happen, within the budget
  set.seed(6)
  dead <- simulate_tree(pars = c(0.05, 2.0), max_t = 10, model = "cr", link = "linear",
                        max_lin = 1000L, max_tries = 3L)
  expect_true(dead$attempts >= 1L && dead$attempts <= 4L)
  if (dead$status == "extinct") expect_identical(dead$attempts, 4L)
})

test_that("a clade that churns under max_lin stops at ten times max_lin lineages in all", {
  # speciation and extinction balance near 100 alive lineages, both near 3:
  # the clade stays well under max_lin = 200 alive while thousands are born and die
  gN <- log(10) / 100   # extinction rises from 0.3 to 3 at 100 lineages, where it meets speciation
  set.seed(7)
  t0 <- proc.time()[["elapsed"]]
  s <- simulate_tree(pars = c(log(3), 0, log(0.3), gN), max_t = 60, model = "dd", link = "exponential",
                     max_lin = 200L, max_tries = 50L)
  expect_lt(proc.time()[["elapsed"]] - t0, 30)
  expect_identical(s$status, "too_large")
  # an ordinary clade is untouched by the table cap
  set.seed(8)
  s2 <- simulate_tree(pars = c(log(0.4), log(0.04)), max_t = 10, model = "cr", link = "exponential",
                      max_lin = 1000L, max_tries = 20L)
  expect_identical(s2$status, "done")
})

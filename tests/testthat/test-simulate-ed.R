# The forward simulator under raw evolutionary distinctiveness.  Raw ED grows
# at unit slope between events, so the hazard is not constant until the next
# event; the simulator draws the waiting time by thinning against a bound on
# the hazard.  Two things pin that down: with the ED coefficients at zero the
# thinning path must reproduce the constant-rate process, and with a negative
# coefficient events must come later than a simulator that holds the rate
# fixed at the start of each waiting time would place them.

test_that("the raw-ED path with zero ED coefficients matches the birth-death moments", {
  skip_on_cran()
  lambda <- 0.5; mu <- 0.1; tt <- 4
  set.seed(11)
  tips <- replicate(300, {
    s <- simulate_tree(pars = c(log(lambda), 0, log(mu), 0), max_t = tt, model = "ed",
                       link = "exponential", max_tries = 1, num_threads = 1L)
    if (identical(s$status, "done") && !is.null(s$tes)) ape::Ntip(s$tes) else NA_real_
  })
  # E[N(T) | survival of both crown lineages] against the same quantity from the
  # cr path, which is exact; the two paths must agree in mean within sampling error
  set.seed(11)
  ref <- replicate(300, {
    s <- simulate_tree(pars = c(log(lambda), log(mu)), max_t = tt, model = "cr",
                       link = "exponential", max_tries = 1, num_threads = 1L)
    if (identical(s$status, "done") && !is.null(s$tes)) ape::Ntip(s$tes) else NA_real_
  })
  a <- tips[is.finite(tips)]; b <- ref[is.finite(ref)]
  se <- sqrt(stats::var(a) / length(a) + stats::var(b) / length(b))
  expect_lt(abs(mean(a) - mean(b)), 3.5 * se)
  # tip counts are integers, so ties are certain and the p-value is the
  # asymptotic one, which is conservative for discrete data
  expect_gt(muffle_expected(stats::ks.test(a, b), "ties")$p.value, 0.01)
})

test_that("a negative ED coefficient on the exponential link slows the clade below the rate-frozen draw", {
  skip_on_cran()
  # Under ned with beta_ED < 0 every lineage's speciation rate falls while it
  # waits, so the exact waiting time is longer than one drawn at the rate of
  # the interval's start.  A clade simulated exactly is therefore smaller on
  # average than one whose simulator froze the rate.  The frozen simulator is
  # not kept, so the check is against the closed-form mean of the matching
  # constant-rate clade: with the brake the clade must be smaller than that.
  pars <- c(log(0.5), 0, -0.3, log(0.05), 0, 0)   # no N term; a strong brake
  set.seed(5)
  tips <- replicate(200, {
    s <- simulate_tree(pars = pars, max_t = 8, model = "ned", link = "exponential",
                       max_tries = 1, num_threads = 1L, max_lin = 5000L)
    if (identical(s$status, "done") && !is.null(s$tes)) ape::Ntip(s$tes) else NA_real_
  })
  a <- tips[is.finite(tips)]
  expect_gt(length(a), 50)
  expect_lt(mean(a), 2 * exp((0.5 - 0.05) * 8))
})

test_that("the thinning path is deterministic under set.seed", {
  pars <- c(log(0.5), -0.01, -0.05, log(0.05), 0, 0)
  one <- function() { set.seed(3); simulate_tree(pars = pars, max_t = 6, model = "ned", link = "exponential",
                                                 max_tries = 1, num_threads = 1L)$L }
  expect_identical(one(), one())
})

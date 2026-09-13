# Thinning envelope of do_augment_tree_cont (audit finding H7).
#
# The sampler proposes speciation times with a homogeneous rate lambda_max on
# each inter-event segment and accepts with probability nh(t) / lambda_max.
# lambda_max must dominate nh(t) on the segment; otherwise the sampler draws
# fewer missing lineages than sampling_prob charges and fhat is biased.
#
# Fixed 7-tip CR tree of dev/audit/checks/H7.R. The C++ engine is seeded from
# the clock, so every assertion is a tolerance on replicates: 3 se for the
# binomial check, 0.03 + 3 se per grid point plus a pooled chi-square for the
# fhat gaps.

brts7 <- c(6, 4.5, 3.0, 2.0, 1.2, 0.5)
T7 <- brts7[1]
s7 <- T7 - brts7[-1]                       # forward branching times

t_ext_tip <- 1e11
t_ext_extinct <- 0
t_ext_unsampled <- 5e10
is_missing_row <- function(df) {
  !(df$t_ext == t_ext_extinct | df$t_ext == t_ext_tip | df$t_ext == t_ext_unsampled)
}

# cumulative hazard int_0^T nh(u) du of the CR rate on the observed tree,
# n(u) = 2 + #{s_i < u}, nh(u) = n(u) lam (1 - exp(-mu (T - u)))
cr_hazard <- function(lam, mu) {
  knots <- c(0, s7, T7)
  nvec <- 2 + seq_along(knots) - 1
  h <- 0
  for (i in seq_len(length(knots) - 1)) {
    a <- knots[i]; b <- knots[i + 1]
    h <- h + nvec[i] * lam * ((b - a) - (1 / mu) * (exp(-mu * (T7 - b)) - exp(-mu * (T7 - a))))
  }
  h
}

draw_cr <- function(lam, mu, n, link = 0L) {
  pars <- if (link == 1L) c(log(lam), 0, 0, 0, log(mu), 0, 0, 0) else c(lam, 0, 0, 0, mu, 0, 0, 0)
  augment_trees(brts = brts7, pars = pars, sample_size = n, maxN = 50L * n,
                max_missing = 200L, max_lambda = 1e6, num_threads = 1L,
                model = c(0L, 0L, 0L), link = link, rho = 1)
}

is_fhat <- function(logf, logg, n_zero = 0) {
  lw <- logf - logg
  m <- max(lw)
  w <- exp(lw - m)
  c(fhat = log(mean(w)) + m - log(1 + n_zero / length(w)),
    se = sd(w) / (mean(w) * sqrt(length(w))))
}

test_that("P(no missing lineage) matches the closed form exp(-int nh) (H7.R part 1)", {
  skip_on_cran()
  thinning_envelope_violations(reset = TRUE)
  n <- 10000L
  raw <- draw_cr(0.4, 0.2, n)
  expect_equal(raw$rejected_overruns + raw$rejected_lambda + raw$rejected_zero_weights, 0)

  n_missing <- vapply(raw$trees, function(df) sum(is_missing_row(df)), numeric(1))
  p0 <- mean(n_missing == 0)
  se <- sqrt(p0 * (1 - p0) / n)
  p0_closed <- exp(-cr_hazard(0.4, 0.2))
  expect_lt(abs(p0 - p0_closed), 3 * se)

  # the density charged for an empty augmentation is the same closed form
  logg_empty <- raw$logg[n_missing == 0]
  expect_equal(logg_empty, rep(-cr_hazard(0.4, 0.2), length(logg_empty)), tolerance = 1e-8)

  expect_equal(thinning_envelope_violations(), 0)
})

test_that("fhat reproduces DDD::bd_loglik across the theta grid (H7.R part 4)", {
  skip_on_cran()
  thinning_envelope_violations(reset = TRUE)
  grid <- rbind(c(0.3, 0.1), c(0.4, 0.2), c(0.5, 0.3), c(0.6, 0.45), c(0.35, 0.3))
  n <- 8000L
  z <- numeric(nrow(grid))
  for (k in seq_len(nrow(grid))) {
    lam <- grid[k, 1]; mu <- grid[k, 2]
    raw <- draw_cr(lam, mu, n)
    f <- is_fhat(raw$logf, raw$logg, raw$rejected_zero_weights)
    # installed DDD: pars2 = c(tdmodel, cond, btorph, verbose, soc)
    ref <- DDD::bd_loglik(pars1 = c(lam, mu, 0, 0), pars2 = c(0, 0, 1, 0, 2),
                          brts = brts7, missnumspec = 0)
    gap <- unname(f["fhat"] - ref)
    z[k] <- gap / unname(f["se"])
    expect_lt(abs(gap), 0.03 + 3 * unname(f["se"]),
              label = sprintf("|fhat - bd_loglik| at lambda = %.2f, mu = %.2f", lam, mu))
  }
  # pooled: five gaps of zero mean give sum(z^2) ~ chi^2(5)
  expect_lt(sum(z^2), qchisq(0.999, df = nrow(grid)))
  expect_equal(thinning_envelope_violations(), 0)
})

test_that("acceptance probability never exceeds 1 on the other rate paths", {
  skip_on_cran()
  thinning_envelope_violations(reset = TRUE)
  # exponential link, constant rates within a segment
  raw <- draw_cr(0.4, 0.2, 2000L, link = 1L)
  expect_equal(raw$rejected_lambda, 0)
  # dd: lambda depends on n only, still constant within a segment
  raw <- augment_trees(brts = brts7, pars = c(0.6, -0.04, 0, 0, 0.2, 0, 0, 0),
                       sample_size = 2000L, maxN = 100000L, max_missing = 200L,
                       max_lambda = 1e6, num_threads = 1L, model = c(1L, 0L, 0L),
                       link = 0L, rho = 1)
  expect_equal(raw$rejected_lambda, 0)
  # M active: lambda varies within a segment, safety-factor envelope
  raw <- augment_trees(brts = brts7, pars = c(0.4, 0, 0.05, 0, 0.2, 0, 0.02, 0),
                       sample_size = 2000L, maxN = 100000L, max_missing = 200L,
                       max_lambda = 1e6, num_threads = 1L, model = c(0L, 1L, 0L),
                       link = 0L, rho = 1)
  expect_equal(raw$rejected_lambda, 0)
  expect_equal(thinning_envelope_violations(), 0)
})

# Thinning envelope of do_augment_tree_cont (audit finding H7).
#
# The sampler proposes speciation times with a homogeneous rate lambda_max on
# each inter-event segment and accepts with probability nh(t) / lambda_max.
# lambda_max must dominate nh(t) on the segment; otherwise the sampler draws
# fewer missing lineages than sampling_prob charges and fhat is biased.
#
# Fixed 7-tip CR tree of dev/audit/checks/H7.R. The C++ engine is seeded from
# the clock, so every assertion is a tolerance on replicates: 3 se for the
# binomial check, 4 se per grid point plus a pooled chi-square for the fhat
# gaps. Measured separation of those assertions on the build before the fix,
# where the envelope collapsed to nh(next_bt): P(no missing) 0.0570 against a
# closed form of 0.0442, 5.5 se; per-point |z| 4.4 / 4.4 / 5.0 / 4.4 / 6.3 and
# 5.3 / 4.6 / 2.1 / 4.4 / 6.0 over two replicates; pooled sum(z^2) 122 and 109
# against a critical value of 20.5. After the fix, over three replicates:
# max |z| 1.86, pooled 4.6 / 5.4 / 1.6.
#
# The envelope is nh at the start of the segment. Since Model::nh_rate holds
# lambda and mu at the segment's own covariates, the only thing carrying t is
# the survival factor, which falls with t: nh is non-increasing on every
# segment and the start value dominates it, for every model and every link.
# The last block holds that directly.

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

  expect_equal(raw$envelope_violations, 0)
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
    # in standard errors, not in log-likelihood units: the earlier
    # |gap| < 0.03 + 3 se allowed 0.05 to 0.10 at n = 8000 and passed on the
    # pre-fix build at all five points.
    expect_lt(abs(z[k]), 4,
              label = sprintf("|fhat - bd_loglik| / se at lambda = %.2f, mu = %.2f", lam, mu))
    expect_equal(raw$envelope_violations, 0)
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
  expect_equal(raw$envelope_violations, 0)
  # dd: lambda depends on n only, still constant within a segment
  raw <- augment_trees(brts = brts7, pars = c(0.6, -0.04, 0, 0, 0.2, 0, 0, 0),
                       sample_size = 2000L, maxN = 100000L, max_missing = 200L,
                       max_lambda = 1e6, num_threads = 1L, model = c(1L, 0L, 0L),
                       link = 0L, rho = 1)
  expect_equal(raw$rejected_lambda, 0)
  expect_equal(raw$envelope_violations, 0)
  # M active: lambda and mu read the mean pendant age, held at the segment's
  # own value, so nh still falls within the segment.  M needs the topology --
  # without it the proposal's P is not reproducible by the scorer and the call
  # is refused -- so supply a parent_tip_start for these branching times.
  raw <- augment_trees(brts = brts7, pars = c(0.4, 0, 0.05, 0, 0.2, 0, 0.02, 0),
                       sample_size = 2000L, maxN = 100000L, max_missing = 200L,
                       max_lambda = 1e6, num_threads = 1L, model = c(0L, 1L, 0L),
                       link = 0L, rho = 1,
                       parent_tip_start = c(rep(0, length(brts7) - 1L), -1))
  expect_equal(raw$rejected_lambda, 0)
  expect_equal(raw$envelope_violations, 0)
  expect_equal(thinning_envelope_violations(), 0)
})

test_that("augment_trees rejects a pars vector shorter than 8", {
  skip_on_cran()
  # a length-7 vector returned trees before the guard in rcpp_mce.
  expect_error(
    augment_trees(brts = brts7, pars = c(0.4, 0, 0, 0, 0.2, 0, 0),
                  sample_size = 5L, maxN = 500L, max_missing = 200L,
                  max_lambda = 1e6, num_threads = 1L, model = c(0L, 0L, 0L),
                  link = 0L, rho = 1),
    "must have length 8")
})

test_that("max_lambda bounds the rate the sampler thins against", {
  skip_on_cran()
  # nd/linear, beta_0 = 1.5, beta_D = -0.5. The envelope is nh at the start of
  # the segment and nothing is inflated on top of it, so max_lambda and the
  # envelope test the same number. A bound of 8 returns every tree; the draw
  # used to be tested against an envelope twice the endpoint rate, where the
  # same bound left 15 trees out of maxN = 20000, i.e. an error, not trees.
  raw <- augment_trees(brts = brts7, pars = c(1.5, 0, 0, -0.5, 0.2, 0, 0, 0),
                       sample_size = 200L, maxN = 20000L, max_missing = 200L,
                       max_lambda = 8, num_threads = 1L, model = c(1L, 0L, 1L),
                       link = 0L, rho = 1)
  expect_equal(raw$num_trees, 200L)
  expect_length(raw$trees, 200L)
})

test_that("the envelope dominates on the D-dependent models too", {
  skip_on_cran()
  # This block used to record a defect. With beta_D active the speciation rate
  # rose within a segment while the survival factor fell, so nh(t) peaked
  # between the endpoints and the safety factor of 2 on the larger endpoint was
  # not a bound: over 200 trees at beta_D = -0.5, 480 to 630 of about 2500
  # segments carried an interior nh above twice the larger endpoint, the worst
  # by a factor of 10^2 to 10^5, while the sampler's own violation count stayed
  # at 0 because on a segment whose endpoint rates are clipped to zero the
  # envelope is zero and no candidate is proposed at all.
  #
  # Model::nh_rate now holds lambda and mu at the covariates of the segment and
  # does not read D, so the only thing carrying t is the survival factor
  # 1 - rho * exp(-mu * (T - t)), which falls with t. nh is therefore
  # non-increasing on every segment, the rate at its start is a dominating
  # envelope, and the safety factor is gone.
  p <- c(0.4, -0.01, 0, -0.5, 0.2, 0, 0, 0)
  thinning_envelope_violations(reset = TRUE)
  raw <- suppressWarnings(
    augment_trees(brts = brts7, pars = p, sample_size = 200L, maxN = 200000L,
                  max_missing = 400L, max_lambda = 1e6, num_threads = 1L,
                  model = c(1L, 0L, 1L), link = 0L, rho = 1))
  expect_equal(raw$num_trees, 200L)
  # the per-call count is the increment of the process-wide counter
  expect_equal(raw$envelope_violations, thinning_envelope_violations(reset = TRUE))
  expect_equal(raw$rejected_lambda, 0)
  expect_equal(raw$envelope_violations, 0)

  over <- 0L; segs <- 0L; worst <- 0
  for (df in raw$trees) {
    for (i in seq_len(nrow(df))) {
      prev <- if (i == 1L) 0 else df$brts[i - 1L]
      if (df$brts[i] <= prev) next
      # the half-open segment (prev, brts]: at t = prev the node lower_bound
      # finds is the one at prev, which governs the segment before this one
      g  <- seq(prev, df$brts[i], length.out = 201L)[-1L]
      nh <- pmax(0, eval_nh_rate(p, df, g, model = c(1L, 0L, 1L),
                                 link = 0L, rho = 1)$nh)
      start <- nh[1L]                         # what the envelope reads
      segs <- segs + 1L
      if (max(nh) > start * (1 + 1e-12)) over <- over + 1L
      if (start > 1e-12) worst <- max(worst, max(nh) / start)
    }
  }
  expect_gt(segs, 1000L)
  expect_identical(over, 0L)
  expect_lt(worst, 1 + 1e-12)
})

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
# The last block records a defect that is still open, not desired behaviour:
# the envelope does not dominate on the D-dependent models.

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
  # M active: lambda varies within a segment, safety-factor envelope
  raw <- augment_trees(brts = brts7, pars = c(0.4, 0, 0.05, 0, 0.2, 0, 0.02, 0),
                       sample_size = 2000L, maxN = 100000L, max_missing = 200L,
                       max_lambda = 1e6, num_threads = 1L, model = c(0L, 1L, 0L),
                       link = 0L, rho = 1)
  expect_equal(raw$rejected_lambda, 0)
  expect_equal(raw$envelope_violations, 0)
  expect_equal(thinning_envelope_violations(), 0)
})

test_that("augment_trees rejects a pars vector shorter than 8", {
  skip_on_cran()
  # the envelope reads pars[7] to decide whether the rates are constant
  # within a segment; a length-7 vector returned trees before the guard.
  expect_error(
    augment_trees(brts = brts7, pars = c(0.4, 0, 0, 0, 0.2, 0, 0),
                  sample_size = 5L, maxN = 500L, max_missing = 200L,
                  max_lambda = 1e6, num_threads = 1L, model = c(0L, 0L, 0L),
                  link = 0L, rho = 1),
    "must have length 8")
})

test_that("max_lambda bounds the rate, not the inflated envelope", {
  skip_on_cran()
  # nd/linear, beta_0 = 1.5, beta_D = -0.5: lambda and mu vary within a
  # segment, so the envelope is the larger endpoint rate times the safety
  # factor 2. Tested against the endpoint rate, a bound of 8 returns every
  # tree (1527 lambda rejections). Tested against the inflated envelope it is
  # the bound 4 that each draw meets, and every draw throws: 19985 lambda
  # rejections and 15 trees against maxN = 20000, i.e. an error, not trees.
  raw <- augment_trees(brts = brts7, pars = c(1.5, 0, 0, -0.5, 0.2, 0, 0, 0),
                       sample_size = 200L, maxN = 20000L, max_missing = 200L,
                       max_lambda = 8, num_threads = 1L, model = c(1L, 0L, 1L),
                       link = 0L, rho = 1)
  expect_equal(raw$num_trees, 200L)
  expect_length(raw$trees, 200L)
})

test_that("the envelope does not dominate on the D-dependent models (recorded, deferred)", {
  skip_on_cran()
  # RECORDED DEFECT, not an accepted behaviour. With beta_D active the
  # speciation rate rises within a segment while the survival factor falls, so
  # nh(t) peaks between the endpoints and the safety factor of 2 on the larger
  # endpoint is not a bound. A dominating envelope needs max(lambda) and
  # max(mu) over the segment separately (see the comment on envelope_safety in
  # src/augment_tree.cpp); deferred to wave 3 with H6, the larger error on the
  # same path.
  #
  # Measured here by scanning nh over each segment of trees the sampler itself
  # drew, rather than by counting its own rejections. The two do not see the
  # same thing: on a segment whose endpoint rates are both clipped to zero the
  # envelope is zero, no candidate is drawn at all, and nothing is counted even
  # though nh is positive inside. Over 200 trees at beta_D = -0.5 and at -2.0,
  # 480 to 630 of about 2500 segments carry an interior nh above twice the
  # larger endpoint, the worst by a factor of 10^2 to 10^5, while the sampler's
  # own violation count is 0. Before Model::nh_rate was given the pendant PD
  # the scorer uses, that count was 8 to 15 per 2000 draws at beta_D = -0.5 and
  # 24 to 46 at -2.0: the corrected P raises the level of lambda far more than
  # its slope, so the candidates the sampler does draw now fall under the
  # inflated envelope and the ones it should have drawn are never proposed.
  # When the envelope is repaired both numbers become 0.
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

  over <- 0L; segs <- 0L; worst <- 0
  for (df in raw$trees) {
    for (i in seq_len(nrow(df))) {
      prev <- if (i == 1L) 0 else df$brts[i - 1L]
      if (df$brts[i] <= prev) next
      g  <- seq(prev, df$brts[i], length.out = 201L)
      nh <- pmax(0, eval_nh_rate(p, df, g, model = c(1L, 0L, 1L),
                                 link = 0L, rho = 1)$nh)
      ends <- max(nh[2L], nh[length(nh)])          # the two the envelope reads
      segs <- segs + 1L
      if (max(nh) > 2 * ends) over <- over + 1L
      if (ends > 1e-12) worst <- max(worst, max(nh) / ends)
    }
  }
  expect_gt(segs, 1000L)
  expect_gt(over, 100L)      # becomes expect_identical(over, 0L) when repaired
  expect_gt(worst, 2)
})

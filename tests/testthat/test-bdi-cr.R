# CR BDI sampler: survival probability, hazard integral and exactness on
# both sides of lam0 = mu0 (audit ids H12, H54).

p_cr   <- emphasis:::.bdi_p_cr
int_cr <- emphasis:::.bdi_integral_cr

tp <- 5

# Numerical integral of the total BDI rate through .bdi_p_cr itself, so an
# error shared by p_cr and the closed form would cancel out of the
# comparison; p_cr is pinned separately against Nee's formula and against
# the hand-written critical limit below.
num_int <- function(t1, t2, n, k, lam, mu, rel.tol = 1e-12) {
  rate <- function(t) {
    p <- vapply(t, function(s) p_cr(s, lam, mu, tp), numeric(1))
    (n + 2 * k) * lam * (1 - p) + n * mu / (1 - p)
  }
  stats::integrate(rate, t1, t2, rel.tol = rel.tol)$value
}

# Critical-case integral written independently of the package code.
crit_int <- function(t1, t2, n, k, lam) {
  a1 <- 1 + lam * (tp - t1)
  a2 <- 1 + lam * (tp - t2)
  I_lam <- lam * (t2 - t1) - log(a1 / a2)
  I_mu  <- lam * (t2 - t1) + log((tp - t1) / (tp - t2))
  (n + 2 * k) * I_lam + n * I_mu
}

# Fixed 11-tip branching-time vector (crown age 5, decreasing).
brts11 <- c(5, 4.745201, 4.530461, 4.067871, 3.622029, 2.929002,
            1.468698, 1.386875, 1.302139, 0.044729)

# --------------------------------------------------------------------------- #
#  .bdi_p_cr                                                                   #
# --------------------------------------------------------------------------- #

test_that(".bdi_p_cr matches Nee's formula for mu < lam and mu > lam", {
  nee <- function(t, lam, mu) {
    E0 <- exp(-(lam - mu) * (tp - t))
    (lam - mu) / (lam - mu * E0)
  }
  for (t in c(0, 1, 2.5, 4.9)) {
    expect_equal(p_cr(t, 0.5, 0.3, tp), nee(t, 0.5, 0.3), tolerance = 1e-14)
    expect_equal(p_cr(t, 0.3, 0.5, tp), nee(t, 0.3, 0.5), tolerance = 1e-14)
  }
  p_hi <- p_cr(1, 0.3, 0.5, tp)
  expect_true(is.finite(p_hi) && p_hi > 0 && p_hi < 1)
})

test_that(".bdi_p_cr at lam == mu is the limit 1/(1 + lam*(tp - t))", {
  for (t in c(0, 1, 2.5, 4.9)) {
    expect_equal(p_cr(t, 0.4, 0.4, tp), 1 / (1 + 0.4 * (tp - t)), tolerance = 1e-14)
  }
  expect_false(is.nan(p_cr(1, 0.4, 0.4, tp)))
  # Continuity across the switch: d = +/-1e-10 agrees with the limit.
  lim <- 1 / (1 + 0.5 * (tp - 1))
  expect_equal(p_cr(1, 0.5, 0.5 - 1e-10, tp), lim, tolerance = 1e-8)
  expect_equal(p_cr(1, 0.5, 0.5 + 1e-10, tp), lim, tolerance = 1e-8)
  # Inside the relative tolerance the limit is used.
  expect_equal(p_cr(1, 0.5, 0.5 - 1e-14, tp), lim, tolerance = 1e-14)
})

# --------------------------------------------------------------------------- #
#  .bdi_integral_cr                                                            #
# --------------------------------------------------------------------------- #

test_that("closed-form CR integral equals integrate() on both sides of lam = mu", {
  cases <- list(c(0.5, 0.3), c(0.3, 0.5), c(0.4, 0.4))
  segs  <- list(c(1, 2), c(0, 4.9))
  for (cs in cases) for (n in c(0L, 1L, 3L)) for (sg in segs) {
    cf <- int_cr(sg[1], sg[2], n, 2L, cs[1], cs[2], tp)
    nm <- num_int(sg[1], sg[2], n, 2L, cs[1], cs[2])
    expect_true(is.finite(cf),
                info = sprintf("lam=%g mu=%g n=%d", cs[1], cs[2], n))
    expect_equal(cf, nm, tolerance = 1e-10,
                 info = sprintf("lam=%g mu=%g n=%d [%g,%g]", cs[1], cs[2], n, sg[1], sg[2]))
  }
})

test_that("critical branch equals the lam -> mu limit (H54 pin)", {
  # t1 = 1, t2 = 3, n = 2, k = 3, lam = 0.5: limit value 8.1425734963
  expect_equal(int_cr(1, 3, 2L, 3L, 0.5, 0.5, tp), crit_int(1, 3, 2L, 3L, 0.5),
               tolerance = 1e-12)
  expect_equal(int_cr(1, 3, 2L, 3L, 0.5, 0.5, tp), 8.1425734963, tolerance = 1e-9)
  # n = 0 uses only the lam-integral.
  expect_equal(int_cr(1, 2, 0L, 2L, 0.4, 0.4, tp), crit_int(1, 2, 0L, 2L, 0.4),
               tolerance = 1e-12)
})

test_that("general branch is continuous with the critical branch at d = +/-1e-10", {
  for (d in c(1e-10, -1e-10)) for (n in c(0L, 1L, 3L)) {
    lam <- 0.5; mu <- lam - d
    cf   <- int_cr(1, 3, n, 2L, lam, mu, tp)
    crit <- crit_int(1, 3, n, 2L, lam)
    expect_true(is.finite(cf))
    expect_equal(cf, crit, tolerance = 1e-8)
    # integrate() through .bdi_p_cr carries cancellation noise at this d.
    expect_equal(cf, num_int(1, 3, n, 2L, lam, mu, rel.tol = 1e-8), tolerance = 1e-7)
  }
})

test_that("integral edge cases: zero-length segment and t2 == tp", {
  expect_identical(int_cr(1, 1, 1L, 2L, 0.5, 0.3, tp), 0)
  # n > 0 and t2 at tp diverges on both sides and at criticality.
  expect_identical(int_cr(1, tp, 1L, 2L, 0.5, 0.3, tp), Inf)
  expect_identical(int_cr(1, tp, 1L, 2L, 0.3, 0.5, tp), Inf)
  expect_identical(int_cr(1, tp, 1L, 2L, 0.4, 0.4, tp), Inf)
  # n = 0 at t2 == tp stays finite.
  expect_true(is.finite(int_cr(1, tp, 0L, 2L, 0.5, 0.3, tp)))
  expect_true(is.finite(int_cr(1, tp, 0L, 2L, 0.3, 0.5, tp)))
  expect_true(is.finite(int_cr(1, tp, 0L, 2L, 0.4, 0.4, tp)))
})

test_that("the integral stays finite when |lam - mu| * (tp - t) passes the exp range", {
  # exp(-(lam0-mu0)*(tp-t)) overflows for mu0 > lam0 once the product passes
  # ~709.  The factored form keeps the exponent non-positive.
  expect_true(is.finite(int_cr(0, 3, 0L, 2L, 0.1, 200, tp)))
  expect_true(is.finite(int_cr(0, 3, 1L, 2L, 0.1, 200, tp)))
  expect_true(is.finite(int_cr(0, 3, 1L, 2L, 0.1, 1e4, tp)))
  # Agreement with the numerical integral, and with the value at mu - lam
  # = 100 where the unfactored form is still in range.
  expect_equal(int_cr(0, 3, 1L, 2L, 0.1, 200, tp),
               num_int(0, 3, 1, 2, 0.1, 200), tolerance = 1e-10)
  expect_equal(int_cr(0, 3, 1L, 2L, 0.1, 100.1, tp), 301.8, tolerance = 1e-9)
  expect_true(is.finite(p_cr(1, 0.1, 200, tp)))

  # And the sampler runs there: lam = 1, mu = exp(5.5) = 244.7 gives
  # (mu - lam) * tp = 1218, well past the overflow threshold.
  set.seed(11)
  a <- emphasis:::.augment_tree_bdi(brts11, pars = c(0, 5.5),
                                    model_bin = c(0L, 0L, 0L),
                                    sample_size = 5L, link = 1L, rho = 1)
  expect_length(a$trees, 5L)
  expect_true(is.finite(a$fhat))
  expect_lt(diff(range(a$weights)), 1e-8)
})

test_that(".bdi_find_event_time_cr runs with mu > lam and n > 0", {
  find_ev <- emphasis:::.bdi_find_event_time_cr
  set.seed(3)
  for (i in 1:20) {
    U <- stats::rexp(1)
    r <- find_ev(1, U, 1L, 2L, 0.3, 0.5, tp, 2)
    expect_true(is.finite(r))
    expect_true(r > 1)
    if (r <= 2) {
      expect_equal(int_cr(1, r, 1L, 2L, 0.3, 0.5, tp), U, tolerance = 1e-8)
    }
  }
})

# --------------------------------------------------------------------------- #
#  Sampler exactness on both sides of lam = mu                                 #
# --------------------------------------------------------------------------- #

test_that("CR BDI sampler gives zero-variance weights for mu < lam, mu > lam, mu = lam", {
  aug <- emphasis:::.augment_tree_bdi
  for (pars in list(c(0.5, 0.3), c(0.3, 0.5), c(0.4, 0.4), c(0.4, 0.41))) {
    set.seed(11)
    a <- aug(brts11, pars = pars, model_bin = c(0L, 0L, 0L),
             sample_size = 30L, link = 0L, rho = 1)
    expect_length(a$trees, 30L)
    expect_true(all(is.finite(a$weights)))
    expect_lt(diff(range(a$weights)), 1e-8)
    expect_true(is.finite(a$fhat))
  }
})

test_that("CR BDI fhat equals the BD likelihood up to -log((n-1)!) on both sides of lam = mu", {
  skip_if_not_installed("DDD")
  aug <- emphasis:::.augment_tree_bdi
  n_tips <- length(brts11) + 1L
  for (pars in list(c(0.5, 0.3), c(0.3, 0.5), c(0.4, 0.4), c(1, 1.3))) {
    set.seed(11)
    a <- aug(brts11, pars = pars, model_bin = c(0L, 0L, 0L),
             sample_size = 30L, link = 0L, rho = 1)
    # pars2 = (tdmodel, cond, btorph, printing, soc)
    bl_brts <- DDD::bd_loglik(pars1 = c(pars, 0, 0), pars2 = c(0, 0, 0, 0, 2),
                              brts = brts11, missnumspec = 0)
    bl_phy  <- DDD::bd_loglik(pars1 = c(pars, 0, 0), pars2 = c(0, 0, 1, 0, 2),
                              brts = brts11, missnumspec = 0)
    expect_equal(a$fhat - bl_brts, -lgamma(n_tips), tolerance = 1e-8,
                 info = sprintf("lam=%g mu=%g", pars[1], pars[2]))
    expect_equal(a$fhat, bl_phy, tolerance = 1e-8,
                 info = sprintf("lam=%g mu=%g", pars[1], pars[2]))
  }
})

test_that("CR BDI sampler is exact under the exponential link with mu > lam", {
  skip_if_not_installed("DDD")
  aug <- emphasis:::.augment_tree_bdi
  set.seed(12)
  a <- aug(brts11, pars = log(c(0.3, 0.5)), model_bin = c(0L, 0L, 0L),
           sample_size = 30L, link = 1L, rho = 1)
  bl <- DDD::bd_loglik(pars1 = c(0.3, 0.5, 0, 0), pars2 = c(0, 0, 0, 0, 2),
                       brts = brts11, missnumspec = 0)
  expect_lt(diff(range(a$weights)), 1e-8)
  expect_equal(a$fhat - bl, -lgamma(length(brts11) + 1L), tolerance = 1e-8)
})

# The widened BDI proposal: the gaussian link and incomplete sampling.
#
# What the widening claims, and what each block here checks:
#
#   1  the rho-sampled survival probability p_rho solves the same ODE with the
#      terminal condition p(tp) = rho, and matches the Stadler reparameterisation
#   2  the cumulative-hazard integrals match quadrature at rho < 1
#   3  EXACTNESS.  Under constant rates, on every supported link and at every
#      rho, the log-weights are constant and fhat is the closed-form likelihood
#   4  the augmentation emits unsampled extant lineages at rho < 1 and none at
#      rho = 1, and the tree data frame counts them in N(t)
#   5  diversity dependence: the DDD gate at rho = 1 is untouched, the mean-field
#      iteration converges at rho < 1 on the monotone links, and does not on the
#      gaussian one -- which is why the gate refuses it
#   6  the gate: what is newly accepted, what is still refused, and the reason
#
# All sampling in .augment_tree_bdi is pure R (stats::rexp / runif), so set.seed
# reproduces every draw.

## Fixed trees.  brts_t11 is the 10-branching tree used across the BDI tests;
## the other two are seeded rphylo draws, pinned here so the file needs no
## simulation at test time.
brts_t11 <- c(5, 4.745201, 4.530461, 4.067871, 3.622029, 2.929002,
              1.468698, 1.386875, 1.302139, 0.044729)

brts_dd20 <- c(6, 4.27428208320388, 4.1258975168475, 3.72207978339226,
               3.64993282835291, 3.17502914508657, 2.81163149722479,
               2.48700670286076, 2.38274529774947, 1.94421117175661,
               1.86664932095232, 1.81164887493753, 1.38503498765192,
               1.06976980095552, 0.460758774034712, 0.417936203054209,
               0.363178692045526, 0.253829586012209, 0.0517782331486725)

cr_bin <- c(0L, 0L, 0L)
dd_bin <- c(1L, 0L, 0L)

# The rho-sampled crown likelihood in the package's labelled-history
# convention, the reference the exactness gate is measured against.  It is the
# rho = 1 form of dev/validation/R/00-common.R::ll_cr_nee with p1 carrying the
# extra factor rho, i.e. Nee's p1 for a Bernoulli-sampled process:
#
#   p1(t) = rho * r^2 * exp(-r t) / [rho*lam + (lam*(1-rho) - mu)*exp(-r t)]^2
#
# log|.| on both r and the denominator so that mu > lam is in range: with
# r < 0 the denominator's zero lies at t < 0, so it keeps one sign over the
# tree, both factors change sign together, and the expression is squared.
#
# At lam = mu the expression is 0/0.  Expanding the denominator in r gives
# r*(1 + rho*lam*t) + O(r^2), so the limit is p1(t) = rho / (1 + rho*lam*t)^2 --
# the complete-sampling limit 1/(1 + lam*t)^2 of dev/validation/00-design.md
# S4.2 with rho carried through.
ll_cr_nee_rho <- function(lam, mu, rho, brts) {
  r <- lam - mu
  logp1 <- if (abs(r) <= 1e-10 * max(lam, mu)) {
    log(rho) - 2 * log1p(rho * lam * brts)
  } else {
    den <- rho * lam + (lam * (1 - rho) - mu) * exp(-r * brts)
    log(rho) + 2 * log(abs(r)) - r * brts - 2 * log(abs(den))
  }
  2 * logp1[1] + sum(log(lam) + logp1[-1])
}

# pars for a constant-rates model reaching (lam, mu) on each link.  On the
# gaussian link the covariate part of eta is zero, so the rate is the constant
# beta_0 * exp(-1/2) and beta_0 = lam * exp(1/2).
cr_pars <- function(lam, mu, link) switch(as.character(link),
  "0" = c(lam, mu), "1" = log(c(lam, mu)), "2" = c(lam, mu) * exp(0.5))


# ---------------------------------------------------------------------------
# 1. the rho-sampled survival probability
# ---------------------------------------------------------------------------

test_that("p_rho solves dp/dt = lam p^2 - (lam-mu) p with p(tp) = rho", {
  tp <- 7
  for (par in list(c(0.5, 0.3), c(0.3, 0.5), c(0.4, 0.4), c(1.2, 0.05))) {
    lam <- par[1]; mu <- par[2]
    for (rho in c(1, 0.9, 0.5, 0.1)) {
      expect_equal(.bdi_p_cr(tp, lam, mu, tp, rho), rho, tolerance = 1e-12)
      ts <- seq(0.1, tp - 0.1, length.out = 60)
      h  <- 1e-5
      num <- sapply(ts, function(t)
        (.bdi_p_cr(t + h, lam, mu, tp, rho) - .bdi_p_cr(t - h, lam, mu, tp, rho)) / (2 * h))
      p   <- sapply(ts, .bdi_p_cr, lam, mu, tp, rho)
      expect_lt(max(abs(num - (lam * p^2 - (lam - mu) * p))), 1e-8)
      # a probability, and monotone toward rho at the present
      expect_true(all(p > 0 & p <= 1))
    }
  }
})

test_that("p_rho equals rho times the Stadler-reparameterised p_1", {
  # p_rho(t; lam, mu) = rho * p_1(t; rho*lam, mu - lam*(1 - rho)).  The factor
  # rho is not cosmetic: p_1 is 1 at tp and p_rho is rho there.  This is the
  # survival-probability face of the likelihood identity
  # nee(lam, mu, rho) = nee(rho*lam, mu - lam*(1-rho), 1) + 2*log(rho), which
  # dev/audit/checks/H2.R verified to 2e-16.
  tp <- 7
  for (par in list(c(0.5, 0.3), c(0.3, 0.5), c(1.0, 0.2))) {
    lam <- par[1]; mu <- par[2]
    for (rho in c(0.9, 0.5, 0.1)) {
      ts <- seq(0, tp, length.out = 80)
      a <- sapply(ts, .bdi_p_cr, lam, mu, tp, rho)
      b <- rho * sapply(ts, .bdi_p_cr, rho * lam, mu - lam * (1 - rho), tp, 1)
      expect_equal(a, b, tolerance = 1e-13)
    }
  }
})

test_that("the rho = 1 path of p and of the hazard integral is untouched", {
  # The complete-sampling branch is the code the validation study measured.
  # Calling it with the default rho and with an explicit 1 must be identical.
  tp <- 7
  for (par in list(c(0.5, 0.3), c(0.3, 0.5), c(0.4, 0.4))) {
    for (t in c(0, 1.5, 6.9)) {
      expect_identical(.bdi_p_cr(t, par[1], par[2], tp),
                       .bdi_p_cr(t, par[1], par[2], tp, 1))
    }
    expect_identical(.bdi_integral_cr(1, 5, 3L, 4L, par[1], par[2], tp),
                     .bdi_integral_cr(1, 5, 3L, 4L, par[1], par[2], tp, 1))
  }
})


# ---------------------------------------------------------------------------
# 2. the cumulative hazard at rho < 1
# ---------------------------------------------------------------------------

test_that(".bdi_integral_cr matches quadrature at rho < 1, both signs of lam-mu", {
  tp <- 7
  for (par in list(c(0.5, 0.3), c(0.3, 0.5), c(0.4, 0.4), c(1.2, 0.05))) {
    lam <- par[1]; mu <- par[2]
    for (rho in c(0.9, 0.5, 0.1)) {
      for (nk in list(c(0L, 2L), c(3L, 5L), c(7L, 2L))) {
        n <- nk[1]; k <- nk[2]
        f <- function(s) {
          omp <- 1 - sapply(s, .bdi_p_cr, lam, mu, tp, rho)
          (n + 2 * k) * lam * omp + n * mu / omp
        }
        # the whole way to tp: at rho < 1 the extinction integral is finite
        # there, which is what lets a doomed lineage survive unsampled
        q <- stats::integrate(f, 1.3, tp, rel.tol = 1e-12,
                              subdivisions = 2000L)$value
        expect_equal(.bdi_integral_cr(1.3, tp, n, k, lam, mu, tp, rho), q,
                     tolerance = 1e-9)
      }
    }
  }
})

test_that("the extinction integral diverges at tp only at rho = 1", {
  # B(tp) = d*(1 - rho): zero at complete sampling (a doomed lineage must die
  # before the present), non-zero below it.
  expect_identical(.bdi_integral_cr(1, 7, 3L, 2L, 0.5, 0.3, 7, 1), Inf)
  expect_true(is.finite(.bdi_integral_cr(1, 7, 3L, 2L, 0.5, 0.3, 7, 0.9)))
  expect_true(is.finite(.bdi_integral_cr(1, 7, 3L, 2L, 0.5, 0.3, 7, 0.5)))
  # with no missing lineage alive there is no extinction term either way
  expect_true(is.finite(.bdi_integral_cr(1, 7, 0L, 2L, 0.5, 0.3, 7, 1)))
})


# ---------------------------------------------------------------------------
# 3. EXACTNESS -- the decisive gate
# ---------------------------------------------------------------------------

test_that("cr: zero-variance weights and fhat = the closed form, every link, every rho", {
  # The gate that holds today for links 0 and 1 at rho = 1 (test-bdi-cr.R,
  # dev/validation/R/00-selfcheck.R assertions c1/c1b) must hold for the
  # gaussian link and at rho < 1 as well.  Measured here: sd(lw) <= 2.4e-12 and
  # |fhat - closed form| <= 2.2e-12 over 3 trees x 3 links x 4 rho x 3 theta.
  trees <- list(t11 = brts_t11, dd20 = brts_dd20)
  for (brts in trees) {
    for (link in c(0L, 1L, 2L)) {
      for (rho in c(1, 0.9, 0.5, 0.2)) {
        for (tgt in list(c(0.5, 0.3), c(0.9, 0.1))) {
          lam <- tgt[1]; mu <- tgt[2]
          set.seed(99)
          e <- .augment_tree_bdi(brts, cr_pars(lam, mu, link), cr_bin,
                                 sample_size = 60L, link = link, rho = rho)
          expect_equal(e$n_valid, 60L)
          expect_equal(e$n_rejected, 0L)        # cr never rejects a survivor
          expect_equal(e$acc, 1)
          expect_lt(stats::sd(e$weights), 1e-9)
          expect_lt(abs(e$fhat - ll_cr_nee_rho(lam, mu, rho, brts)), 1e-8)
        }
      }
    }
  }
})

test_that("cr at mu > lam and at lam = mu is exact too, at every rho", {
  # mu >= lam is the branch the closed forms are written to keep in range
  # (H12); the critical case lam = mu is the 0/0 limit of both.
  for (tgt in list(c(0.3, 0.5), c(0.4, 0.4))) {
    lam <- tgt[1]; mu <- tgt[2]
    for (link in c(0L, 1L, 2L)) for (rho in c(1, 0.5)) {
      set.seed(3)
      e <- .augment_tree_bdi(brts_t11, cr_pars(lam, mu, link), cr_bin,
                             sample_size = 40L, link = link, rho = rho)
      expect_lt(stats::sd(e$weights), 1e-9)
      expect_lt(abs(e$fhat - ll_cr_nee_rho(lam, mu, rho, brts_t11)), 1e-8)
    }
  }
})

test_that("the gaussian link is the other two under a reparameterisation", {
  # beta_0 * exp(-1/2) = lam, so the same (lam, mu) reached through any link
  # gives the same likelihood.  This is what makes cr on the gaussian link
  # exact rather than merely available.
  for (rho in c(1, 0.5)) {
    f <- sapply(c(0L, 1L, 2L), function(link) {
      set.seed(21)
      .augment_tree_bdi(brts_dd20, cr_pars(0.7, 0.25, link), cr_bin,
                        sample_size = 40L, link = link, rho = rho)$fhat
    })
    expect_equal(diff(range(f)), 0, tolerance = 1e-9)
  }
})


# ---------------------------------------------------------------------------
# 4. unsampled extant lineages
# ---------------------------------------------------------------------------

test_that("the augmentation emits unsampled extant lineages at rho < 1 and none at rho = 1", {
  # H2 measured 0 unsampled-extant sentinels over 30 draws at rho = 0.5 while
  # the thinning sampler wrote up to 17.  That is the defect being closed.
  n_uns <- function(e) sapply(e$trees, function(tr) sum(tr$t_ext == 5e10))
  set.seed(8)
  e1 <- .augment_tree_bdi(brts_dd20, c(0.5, 0.3), cr_bin, sample_size = 60L,
                          link = 0L, rho = 1)
  expect_equal(max(n_uns(e1)), 0L)
  prev <- 0
  for (rho in c(0.9, 0.5, 0.2)) {
    set.seed(8)
    e <- .augment_tree_bdi(brts_dd20, c(0.5, 0.3), cr_bin, sample_size = 60L,
                           link = 0L, rho = rho)
    m <- mean(n_uns(e))
    expect_gt(m, prev)          # fewer sampled tips, more unsampled lineages
    prev <- m
  }
  expect_gt(prev, 5)
})

test_that("an unsampled extant lineage raises N(t) and never lowers it", {
  set.seed(4)
  e <- .augment_tree_bdi(brts_t11, c(0.6, 0.2), cr_bin, sample_size = 30L,
                         link = 0L, rho = 0.4)
  for (tr in e$trees) {
    # n is rebuilt from the sentinels: +1 on every non-extinction node
    expect_equal(tr$n, c(2, utils::head(2 + cumsum(ifelse(tr$t_ext == 0, -1, 1)),
                                        -1L)))
    uns <- tr[tr$t_ext == 5e10, , drop = FALSE]
    if (nrow(uns) > 0L) {
      # an unsampled extant has a birth node and no extinction node: its id
      # appears exactly once
      expect_true(all(table(tr$id[tr$id %in% uns$id]) == 1L))
      expect_true(all(uns$brts < brts_t11[1]))
      expect_true(all(uns$tip_start == uns$brts))
    }
  }
  # the closing node still marks the present
  expect_true(all(sapply(e$trees, function(tr) utils::tail(tr$brts, 1L)) ==
                    brts_t11[1]))
})

test_that("rho < 1 removes the survivor rejection channel entirely", {
  # At rho = 1 a surviving missing lineage has f = 0 and is rejected, and fhat
  # carries log(acc) for it.  At rho < 1 it is a legitimate configuration, so
  # acc is 1 and the correction is inert.
  set.seed(12)
  e1 <- .augment_tree_bdi(brts_dd20, c(0.8, 0.3, -0.02, 0)[c(1, 3, 2, 4)],
                          dd_bin, sample_size = 120L, link = 0L, rho = 1)
  expect_gt(e1$n_rejected, 0L)
  expect_lt(e1$acc, 1)
  set.seed(12)
  e2 <- .augment_tree_bdi(brts_dd20, c(0.8, 0.3, -0.02, 0)[c(1, 3, 2, 4)],
                          dd_bin, sample_size = 120L, link = 0L, rho = 0.5)
  expect_equal(e2$n_rejected, 0L)
  expect_equal(e2$acc, 1)
})


# ---------------------------------------------------------------------------
# 5. diversity dependence
# ---------------------------------------------------------------------------

test_that("the mean-field iteration converges for linear and exponential dd, rho included", {
  bt <- sort(brts_dd20[1] - brts_dd20[-1]); tp <- brts_dd20[1]
  for (rho in c(1, 0.9, 0.5, 0.2)) {
    for (K in c(20, 50, 1e4)) {
      s <- .bdi_iterate(.expand_pars(c(0.8, -(0.8 - 0.3) / K, 0.3, 0), dd_bin),
                        dd_bin, 0L, bt, tp, rho = rho)
      expect_true(s$converged)
      expect_lt(s$delta, 1e-4)
    }
    for (bN in c(-0.02, -0.05, -0.15)) {
      s <- .bdi_iterate(.expand_pars(c(log(0.8), bN, log(0.2), 0), dd_bin),
                        dd_bin, 1L, bt, tp, rho = rho)
      expect_true(s$converged)
    }
  }
})

test_that("the sweep budget is what small rho needs, not what rho = 1 needs", {
  # The iteration slows as rho falls, because the unsampled extant lineages
  # raise Nhat, which feeds back into the rates that produced them.  At
  # K = 20 on this tree it takes 6 sweeps at rho = 1 and 27 at rho = 0.2, so
  # the old cap of 20 returned a non-converged mean field there.
  bt <- sort(brts_dd20[1] - brts_dd20[-1]); tp <- brts_dd20[1]
  p  <- .expand_pars(c(0.8, -(0.8 - 0.3) / 20, 0.3, 0), dd_bin)
  its <- vapply(c(1, 0.9, 0.5, 0.2), function(rho)
    .bdi_iterate(p, dd_bin, 0L, bt, tp, rho = rho)$iterations, integer(1))
  expect_true(all(diff(its) > 0))       # monotone in decreasing rho
  expect_gt(max(its), 20L)              # past the old cap
  expect_lt(max(its), 60L)              # inside the new one
  # the cap only ever truncates: a run that converged is unaffected by it
  a <- .bdi_iterate(p, dd_bin, 0L, bt, tp, rho = 1)
  b <- .bdi_iterate(p, dd_bin, 0L, bt, tp, rho = 1, max_iter = 200L)
  expect_identical(a$iterations, b$iterations)
  expect_identical(a$delta, b$delta)
})

test_that("the mean-field iteration DIVERGES for gaussian dd, and not merely slowly", {
  # lambda(N) = beta_0*exp(-(beta_N*N - 1)^2/2) peaks at N = 1/beta_N and falls
  # after it, so the map N -> rate -> N can oscillate past the peak.  This is
  # the measurement behind .bdi_supported refusing dd on the gaussian link.
  #
  # The budget is deliberately 200 here -- more than three times what
  # .bdi_iterate allows -- so that the failure cannot be read as "the cap is
  # too low".  delta grows into the tens: these are divergent, not slow.
  bt <- sort(brts_dd20[1] - brts_dd20[-1]); tp <- brts_dd20[1]
  gp <- function(b0, bN) .expand_pars(c(b0 * exp(0.5), bN, 0.3 * exp(0.5), 0),
                                      dd_bin)
  # `auto` is the relaxed sweep.  It rescues the first two cases and only
  # shrinks the third, which is why the gate on this link stays shut: a
  # relaxation that fixes most of a region is not a licence to admit the region.
  for (cse in list(list(b0 = 4, bN = 0.20, rho = 1,   auto_ok = TRUE),
                   list(b0 = 2, bN = 0.10, rho = 0.5, auto_ok = TRUE),
                   list(b0 = 4, bN = 0.10, rho = 0.2, auto_ok = FALSE))) {
    # damping = 1 is the plain Picard iteration, which is what this measures.
    s <- .bdi_iterate(gp(cse$b0, cse$bN), dd_bin, 2L, bt, tp, rho = cse$rho,
                      max_iter = 200L, damping = 1)
    expect_false(s$converged)
    expect_gt(s$delta, 1)
    expect_identical(s$iterations, 200L)
    sd <- .bdi_iterate(gp(cse$b0, cse$bN), dd_bin, 2L, bt, tp, rho = cse$rho,
                       max_iter = 200L, damping = "auto")
    expect_identical(sd$converged, cse$auto_ok)
    # and it never leaves the residual larger than the plain iteration did
    expect_lte(sd$delta, s$delta)
  }
  # And the failing region is not an interval that could be carved out: at
  # rho = 0.5, beta_0 = 1 the iteration converges at beta_N = 0.08, fails at
  # 0.1, and converges again at 0.2.  A static gate on the link is the only
  # honest one.
  conv <- vapply(c(0.08, 0.1, 0.2), function(bN)
    .bdi_iterate(gp(1, bN), dd_bin, 2L, bt, tp, rho = 0.5,
                 max_iter = 200L, damping = 1)$converged, logical(1))
  expect_identical(conv, c(TRUE, FALSE, TRUE))
  # Relaxed, the hole closes: all three converge.  So "the failing region is
  # not an interval" is a statement about the plain iteration only.
  convd <- vapply(c(0.08, 0.1, 0.2), function(bN)
    .bdi_iterate(gp(1, bN), dd_bin, 2L, bt, tp, rho = 0.5,
                 max_iter = 200L, damping = "auto")$converged, logical(1))
  expect_identical(convd, c(TRUE, TRUE, TRUE))

  # A non-converged mean field is announced rather than used silently -- the
  # gate keeps dd/gaussian out of estimate_rates, but .augment_tree_bdi is
  # callable directly and must not stay quiet about it.  With the relaxation
  # off, so that there is a non-convergence left to announce.
  expect_warning(
    .augment_tree_bdi(brts_dd20, c(4 * exp(0.5), 0.1, 0.3 * exp(0.5), 0), dd_bin,
                      sample_size = 5L, link = 2L, rho = 0.2, damping = 1),
    "mean-field iteration did not converge")
})

test_that("dd at rho = 1 on the linear link still matches DDD::dd_loglik", {
  # The existing gate, re-run here: the widening must not move it.
  skip_if_not_installed("DDD")
  dd_pars <- function(l0, m0, K) c(l0, -(l0 - m0) / K, m0, 0)
  grid <- expand.grid(K = c(20, 50, 1e4), i = 1:3)
  lm <- rbind(c(0.8, 0.3), c(0.6, 0.1), c(1.0, 0.8))
  set.seed(101)
  gap <- vapply(seq_len(nrow(grid)), function(j) {
    l0 <- lm[grid$i[j], 1]; m0 <- lm[grid$i[j], 2]; K <- grid$K[j]
    ref <- DDD::dd_loglik(pars1 = c(l0, m0, K), pars2 = c(300, 1, 0, 1, 0, 2),
                          brts = brts_dd20, missnumspec = 0)
    e <- .augment_tree_bdi(brts_dd20, dd_pars(l0, m0, K), dd_bin,
                           sample_size = 500L, link = 0L, rho = 1)
    expect_true(e$mf_converged)
    e$fhat - ref
  }, numeric(1))
  expect_lt(diff(range(gap)), 0.25)
  expect_lt(max(abs(gap)), 0.3)
})

test_that("dd at rho < 1 draws unsampled lineages and keeps a usable ESS", {
  # No DDD reference exists at rho < 1 (dd_loglik takes missnumspec, a fixed
  # count, not a Bernoulli fraction), so this pins behaviour rather than a
  # likelihood value: the draws carry unsampled extant lineages, the weights
  # are finite, and the proposal does not collapse.
  #
  # Measured over seeds 11-16 at N = 200: ESS 40-71 (20-35%) at rho = 0.5
  # against 143-155 (72-78%) at rho = 1.  The mean-field proposal is worse at
  # rho < 1 than at rho = 1 -- the augmented N is larger and further from its
  # mean-field value -- and still far above the thinning proposal's few
  # per cent in these regimes.
  ess <- nonfin <- numeric(0)
  for (sd in 11:16) {
    set.seed(sd)
    e <- .augment_tree_bdi(brts_dd20, c(0.8, -0.02, 0.3, 0), dd_bin,
                           sample_size = 200L, link = 0L, rho = 0.5)
    expect_true(e$mf_converged)
    expect_equal(e$n_valid, 200L)
    expect_equal(e$n_rejected, 0L)
    expect_gt(mean(sapply(e$trees, function(tr) sum(tr$t_ext == 5e10))), 0)
    ess    <- c(ess, .ess_from_lw(e$weights))
    nonfin <- c(nonfin, e$n_nonfinite)
  }
  expect_gt(min(ess), 0.1 * 200)
  expect_gt(stats::median(ess), 0.2 * 200)
  # lambda(N) = 0.8 - 0.02 N hits the linear link's clamp at N = 40, and the
  # extra unsampled lineages take the augmented N there: a handful of draws
  # per E-step score -Inf and leave the M-step set.  Reaching the clamp is a
  # property of a linear rate at rho < 1, not a defect of the sampler, but it
  # must stay rare.
  expect_lt(max(nonfin), 0.05 * 200)
})


# ---------------------------------------------------------------------------
# 6. the gate
# ---------------------------------------------------------------------------

test_that(".bdi_supported accepts the widened scope", {
  for (rho in c(1, 0.9, 0.5, 1e-3)) {
    expect_true(.bdi_supported(c(0L, 0L, 0L), 0L, rho))   # cr linear
    expect_true(.bdi_supported(c(0L, 0L, 0L), 1L, rho))   # cr exponential
    expect_true(.bdi_supported(c(0L, 0L, 0L), 2L, rho))   # cr gaussian  (new)
    expect_true(.bdi_supported(c(1L, 0L, 0L), 0L, rho))   # dd linear
    expect_true(.bdi_supported(c(1L, 0L, 0L), 1L, rho))   # dd exponential
  }
})

test_that(".bdi_supported still refuses what is unsupported, and names the reason", {
  # D-dependent: out of reach by construction, not by omission.
  for (link in c(0L, 1L, 2L)) for (rho in c(1, 0.5)) {
    expect_false(.bdi_supported(c(0L, 0L, 1L), link, rho))
    expect_false(.bdi_supported(c(1L, 0L, 1L), link, rho))
    expect_false(.bdi_supported(c(0L, 1L, 0L), link, rho))   # M-dependent
  }
  # dd on the gaussian link: the Picard iteration does not converge there.
  expect_false(.bdi_supported(c(1L, 0L, 0L), 2L, 1))
  expect_false(.bdi_supported(c(1L, 0L, 0L), 2L, 0.5))
  # rho out of range
  expect_false(.bdi_supported(c(0L, 0L, 0L), 0L, 0))
  expect_false(.bdi_supported(c(0L, 0L, 0L), 0L, 1.5))
  expect_false(.bdi_supported(c(0L, 0L, 0L), 0L, NA_real_))

  expect_match(.bdi_unsupported_reason(c(0L, 0L, 1L), 0L, 1), "D-dependent")
  expect_match(.bdi_unsupported_reason(c(0L, 0L, 1L), 0L, 1), "pendant age")
  expect_match(.bdi_unsupported_reason(c(0L, 1L, 0L), 0L, 1), "M-dependent")
  expect_match(.bdi_unsupported_reason(c(1L, 0L, 0L), 2L, 1), "gaussian")
  expect_match(.bdi_unsupported_reason(c(1L, 0L, 0L), 2L, 1), "not monotone")
  expect_match(.bdi_unsupported_reason(c(0L, 0L, 0L), 0L, 2), "outside")
  expect_null(.bdi_unsupported_reason(c(0L, 0L, 0L), 2L, 0.5))
  expect_null(.bdi_unsupported_reason(c(1L, 0L, 0L), 1L, 0.5))
})

test_that("estimate_rates routes to BDI where the gate now allows it", {
  # rho < 1 used to announce a fallback; it must not any more.
  expect_message(
    fit <- estimate_rates(brts_t11, model = "cr", init_pars = c(0.5, 0.2),
                          method = "mcem",
                          control = list(sampling = "bdi", num_trees = 20L,
                                         max_iter = 3L, rho = 0.5,
                                         lower_bound = c(1e-3, 0),
                                         upper_bound = c(3, 3))),
    NA)
  expect_equal(fit$rho, 0.5)
  # the gaussian link for cr, likewise
  expect_message(
    estimate_rates(brts_t11, model = "cr", link = "gaussian",
                   init_pars = c(0.5, 0.2) * exp(0.5), method = "mcem",
                   control = list(sampling = "bdi", num_trees = 20L,
                                  max_iter = 3L, rho = 1,
                                  lower_bound = c(1e-3, 0) * exp(0.5),
                                  upper_bound = c(3, 3) * exp(0.5))),
    NA)
})

test_that("estimate_rates still falls back for dd on the gaussian link, with the reason", {
  expect_message(
    estimate_rates(brts_t11, model = "dd", link = "gaussian",
                   init_pars = c(0.8, 0.02, 0.2, 0), method = "mcem",
                   control = list(sampling = "bdi", num_trees = 10L,
                                  max_iter = 2L, maxN = 500L,
                                  lower_bound = c(1e-3, 0, 0, 0),
                                  upper_bound = c(3, 0.5, 3, 0))),
    "gaussian link")
})


# ---------------------------------------------------------------------------
# 7. the complete-sampling path the validation study measured is unchanged
# ---------------------------------------------------------------------------

test_that("cr and dd at rho = 1 on links 0 and 1 take the untouched branch", {
  # The bit-for-bit comparison against the pre-change build is not a fixture in
  # the repository: it was run at development time by capturing
  # .augment_tree_bdi's output and its eval_logf rescoring over 18 model x tree
  # cells x 12 seeded draws on the pristine parent commit and requiring
  # identical() afterwards (it held).  What can be checked from inside the
  # package is the invariant that makes that comparison hold at all: rho = 1
  # must take the complete-sampling branch, so the draw cannot depend on
  # whether rho was passed explicitly, and no unsampled sentinel is written.
  grid <- list(list(mb = cr_bin, link = 0L, pars = c(0.5, 0.3)),
               list(mb = cr_bin, link = 1L, pars = log(c(0.6, 0.2))),
               list(mb = dd_bin, link = 0L, pars = c(0.8, -0.01, 0.1, 0)),
               list(mb = dd_bin, link = 1L, pars = c(-0.3, -0.01, -1.5, 0)))
  for (g in grid) {
    set.seed(4242)
    a <- suppressWarnings(.augment_tree_bdi(brts_t11, g$pars, g$mb,
                                            sample_size = 10L, link = g$link))
    set.seed(4242)
    b <- suppressWarnings(.augment_tree_bdi(brts_t11, g$pars, g$mb,
                                            sample_size = 10L, link = g$link,
                                            rho = 1))
    expect_identical(a$trees, b$trees)
    expect_identical(a$logg, b$logg)
    expect_identical(a$fhat, b$fhat)
    # and no unsampled sentinel is ever written at complete sampling
    expect_equal(sum(sapply(a$trees, function(tr) sum(tr$t_ext == 5e10))), 0)
  }
})

# ---------------------------------------------------------------------------
# The mean-field sweep budget is part of the proposal, not a free parameter:
# raising it changes the mean field of any run that had not reached tolerance
# within it, and with it the draws and the estimate.  At rho = 1 the budget is
# the one the validation study measured and must not move.
# ---------------------------------------------------------------------------

test_that("the sweep budget at complete sampling is 20 and pins the proposal", {
  skip_on_cran()
  expect_equal(eval(formals(emphasis:::.bdi_iterate)$max_iter), NULL)
  # resolved value: 20 at rho = 1, larger below it
  bt <- sort(c(0.4, 1.1, 1.9, 2.6, 3.4, 4.0))
  p8 <- c(1.5, -0.06, 0, 0, 0.3, 0, 0, 0)          # dd, linear: a slow cell
  sol20  <- emphasis:::.bdi_iterate(p8, c(1L,0L,0L), 0L, bt, 5, rho = 1)
  sol200 <- emphasis:::.bdi_iterate(p8, c(1L,0L,0L), 0L, bt, 5, max_iter = 200L,
                                    rho = 1)
  expect_lte(sol20$iterations, 20L)
  # if this cell converges inside 20 the two agree; if it does not, they differ,
  # and that difference is exactly why the budget is pinned
  if (!isTRUE(sol20$converged)) {
    expect_gt(abs(sol20$delta - sol200$delta), 0)
  }
  # below rho = 1 the budget is larger, because the fixed point is slower there
  solr <- emphasis:::.bdi_iterate(p8, c(1L,0L,0L), 0L, bt, 5, rho = 0.2)
  expect_gt(solr$iterations, 0L)
  expect_lte(solr$iterations, 200L)
})

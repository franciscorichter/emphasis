# log q is the log density of what the thinning sampler drew (audit finding H6).
#
# The proposal, as it now stands:
#
#   births     an inhomogeneous Poisson process on (0, T) of intensity
#                nh(t) = N * lambda_seg * (1 - rho * exp(-mu_seg * (T - t)))
#              where N is the number of lineages alive, and lambda_seg, mu_seg
#              are the speciation and extinction rates of the model read on the
#              covariates (N, M) of the segment t falls in, with the mean pendant
#              age M held at the value it has where the segment begins.  The
#              proposal does not read D.  Only the survival factor carries t, and
#              it falls with t, so nh is decreasing on every segment: the rate at
#              the segment's start is a dominating thinning envelope, and
#              int nh dt is the closed form sampling_prob already used.
#   lifetime   exponential with rate mu_seg, truncated to (0, T - t), or (with
#              rho < 1) unsampled-extant with the complementary probability.
#   parent     uniform over the labelled attachments, 2 * tips + Ne of them (H5).
#
# Four things are held here.
#
#   2  the lifetimes are drawn from the density they are charged with (a KS test
#      of their probability integral transform).
#   3  E_q[f/q] is the marginal likelihood: on a 2-tip tree, against a
#      brute-force integral over augmentations with one missing lineage.
#   4  a d/nd model with its D coefficients at zero charges exactly what the cr
#      or dd model it extends charges.
#   +  the envelope never fails to dominate on a D model.
#
# The cr/dd bit-for-bit gate is dev/crdd_invariance.R: it has to carry a corpus
# of augmented trees across two builds, which a test cannot.

T_TIP <- 10e10; T_EXT <- 0; T_UNS <- 5e10
is_mis_row <- function(df) !(df$t_ext == T_EXT | df$t_ext == T_TIP | df$t_ext == T_UNS)
.lnk <- function(link, eta) if (link == 1L) exp(eta) else pmax(0, eta)

# The extinction rate sampling_prob charges a lifetime with: the segment's mu.
# Written here from the definition above, not from the C++.
.mu_charged <- function(df, i, p, link) {
  prev <- if (i == 1L) 0 else df$brts[i - 1L]
  pd <- df$pd[i] + df$n[i] * (prev - df$brts[i])     # P where the segment begins
  M  <- if (df$n[i] > 0) pd / df$n[i] else 0
  max(.lnk(link, p[5] + p[6] * df$n[i] + p[7] * M), 1e-10)
}

# A tree deep enough that mu * (T - t) is of order 1: on a shallow tree the
# truncated exponential is close to uniform on (0, T - t) whatever its rate, and
# the transform cannot see a wrong rate at all.
.deep_tree <- function() {
  set.seed(11)
  phy <- ape::rphylo(10L, 0.8, 0.25)
  phy$edge.length <- phy$edge.length * 4
  phy
}

.d_settings <- list(
  list(name = "d  linear", mb = c(0L, 0L, 1L), link = 0L,
       p = c( 0.12,  0.000, 0, 0.05,  0.45, 0, 0, 0.40)),
  list(name = "nd linear", mb = c(1L, 0L, 1L), link = 0L,
       p = c( 0.14, -0.004, 0, 0.05,  0.45, 0, 0, 0.40)),
  list(name = "d  exp",    mb = c(0L, 0L, 1L), link = 1L,
       p = c(-2.10,  0.000, 0, 0.05, -0.80, 0, 0, 0.60)),
  list(name = "nd exp",    mb = c(1L, 0L, 1L), link = 1L,
       p = c(-2.05, -0.004, 0, 0.05, -0.80, 0, 0, 0.60))
)


# --------------------------------------------------------------------------- #
#  Gate 2: the lifetimes come from the density they are charged with           #
# --------------------------------------------------------------------------- #

test_that("the drawn lifetimes pass a KS test against the density logg charges", {
  skip_on_cran()
  phy  <- .deep_tree()
  brts <- emphasis:::.extract_brts(phy)
  pts  <- emphasis:::.pts(brts)
  TT   <- brts[1]
  # On the build this test was written against, the same four settings gave
  # p = 3e-55 (d linear), 5e-96 (nd linear), 9e-05 (d exp), 1e-04 (nd exp) at
  # 2000 draws: the lifetime was drawn with the segment's D-free rate and
  # charged with the D-aware rate of the inserted node.
  for (s in .d_settings) {
    a <- augment_trees(as.numeric(brts), s$p, 3000L, 600000L, 800L, 1e6, 1L,
                       model = s$mb, link = s$link, rho = 1, parent_tip_start = pts)
    expect_equal(a$envelope_violations, 0)
    mu <- l <- r <- numeric(0)
    for (df in a$trees) {
      for (i in which(is_mis_row(df))) {
        mu <- c(mu, .mu_charged(df, i, s$p, s$link))
        l  <- c(l, df$t_ext[i] - df$brts[i])
        r  <- c(r, TT - df$brts[i])
      }
    }
    expect_gt(length(l), 3000)
    expect_gt(median(mu * r), 0.5)          # the transform can see the rate
    u <- (1 - exp(-mu * l)) / (1 - exp(-mu * r))
    expect_gt(suppressWarnings(stats::ks.test(u, "punif"))$p.value, 1e-4,
              label = sprintf("KS p of the lifetime transform, %s", s$name))
    # and a statistic the KS tail cannot hide: the transform has mean 1/2
    expect_lt(abs(mean(u) - 0.5) / (stats::sd(u) / sqrt(length(u))), 4,
              label = sprintf("|mean(PIT) - 1/2| / se, %s", s$name))
  }
})


# --------------------------------------------------------------------------- #
#  Gate 3: E_q[f/q] is the marginal likelihood                                 #
# --------------------------------------------------------------------------- #
#
# On a 2-tip tree an augmentation with one missing lineage is a pair (t, d).
# The observed tree has no internal node, so the sampler's candidate-parent list
# is empty at every t and the parent is a crown lineage -- which is also every
# labelled attachment the -log(2 tips + Ne) term counts: 2 per observed lineage,
# and here both observed lineages are crown lineages at tip_start 0, so f is the
# same for all four.  The stratum is therefore exactly
#
#   L1 = 4 * int_0^T int_t^T f(y, z(t, d)) dd dt.
#
# On a tree with an observed split it would not be: the sampler draws the parent
# uniformly over the nodes alive, which never include the crown lineages, so some
# labelled attachments have proposal probability zero (H45, open).  That is why
# the brute force is done on a 2-tip tree.

.gauss_legendre <- function(n) {                   # Golub-Welsch, mapped to [0,1]
  k <- 1:(n - 1); b <- k / sqrt(4 * k^2 - 1)
  J <- diag(0, n); J[cbind(k, k + 1)] <- b; J[cbind(k + 1, k)] <- b
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = (e$values[o] + 1) / 2, w = e$vectors[1, o]^2)
}

# The augmented 2-tip tree the sampler builds for a birth at t dying at d.
# P is the pendant PD of the alive set: {0, 0} before t, {0, t, t} on (t, d],
# {0, t} after -- the split resets the parent's tip start as well as the
# daughter's.  Checked against a drawn tree below.
.aug1 <- function(t, d, TT) data.frame(
  brts = c(t, d, TT), n = c(2, 3, 2), t_ext = c(d, 0, T_TIP),
  pd = c(2 * t, 3 * d - 2 * t, 2 * TT - t), tip_start = c(t, t, TT),
  focal_tip_start = c(0, t, -1), id = c(1L, 1L, 0L), parent_id = c(-1L, -1L, -1L))

.aug0 <- function(TT) data.frame(
  brts = TT, n = 2, t_ext = T_TIP, pd = 2 * TT,
  tip_start = TT, focal_tip_start = -1, id = 0L, parent_id = -1L)

.L1 <- function(TT, p, mb, lk, n = 60L) {
  gq <- .gauss_legendre(n); tt <- TT * gq$x
  trees <- vector("list", n * n); wgt <- numeric(n * n); k <- 0L
  for (ix in seq_len(n)) {
    t <- tt[ix]; dd <- t + (TT - t) * gq$x
    for (iy in seq_len(n)) {
      k <- k + 1L
      trees[[k]] <- .aug1(t, dd[iy], TT)
      wgt[k] <- gq$w[ix] * gq$w[iy] * TT * (TT - t)
    }
  }
  4 * sum(wgt * exp(eval_logf(p, trees, model = mb, link = lk, rho = 1)$logf))
}

test_that("the 2-tip augmented tree is the one the sampler builds", {
  skip_on_cran()
  TT <- 3
  p  <- c(0.30, 0, 0, 0.15, 0.25, 0, 0, 0.10)
  a <- augment_trees(TT, p, 400L, 40000L, 50L, 1e6, 1L, model = c(0L, 0L, 1L),
                     link = 0L, rho = 1, parent_tip_start = c(-1))
  nm <- vapply(a$trees, function(d) sum(is_mis_row(d)), 0)
  skip_if(!any(nm == 1L), "no draw with exactly one missing lineage")
  df <- a$trees[[which(nm == 1L)[1L]]]
  ref <- .aug1(df$brts[1], df$brts[2], TT)
  for (col in names(ref)) expect_equal(as.numeric(df[[col]]), as.numeric(ref[[col]]),
                                       tolerance = 1e-12, info = col)
  df0 <- a$trees[[which(nm == 0L)[1L]]]
  ref0 <- .aug0(TT)
  for (col in names(ref0)) expect_equal(as.numeric(df0[[col]]), as.numeric(ref0[[col]]),
                                        tolerance = 1e-12, info = col)
  # and its logg is the compensator alone: -int_0^T nh dt, in closed form
  lam <- 0.30; mu <- 0.25
  expect_equal(a$logg[which(nm == 0L)[1L]],
               -2 * lam * (TT - (1 / mu) * (1 - exp(-mu * TT))), tolerance = 1e-10)
})

test_that("E_q[f/q] is the brute-force marginal likelihood on a 2-tip tree", {
  skip_on_cran()
  TT <- 3
  for (cfg in list(list(name = "d linear", mb = c(0L, 0L, 1L), link = 0L,
                        p = c(0.30, 0, 0, 0.15, 0.25, 0, 0, 0.10)),
                   list(name = "d exp", mb = c(0L, 0L, 1L), link = 1L,
                        p = c(log(0.30), 0, 0, 0.15, log(0.25), 0, 0, 0.10)))) {
    l0 <- exp(eval_logf(cfg$p, list(.aug0(TT)), model = cfg$mb, link = cfg$link,
                        rho = 1)$logf)
    l1 <- .L1(TT, cfg$p, cfg$mb, cfg$link, 60L)
    # the quadrature is converged: halving the order moves it by < 1e-5 relative
    expect_lt(abs(l1 - .L1(TT, cfg$p, cfg$mb, cfg$link, 30L)) / l1, 1e-5)

    tot1 <- ss1 <- totle <- ssle <- w2 <- 0; N <- 0L; nz <- 0L
    for (b in 1:3) {
      a <- augment_trees(TT, cfg$p, 20000L, 2000000L, 200L, 1e6, 1L, model = cfg$mb,
                         link = cfg$link, rho = 1, parent_tip_start = c(-1))
      expect_equal(a$rejected_overruns, 0)
      expect_equal(a$rejected_lambda, 0)
      expect_equal(a$rejected_nonfinite, 0)
      expect_equal(a$envelope_violations, 0)
      # a draw with f = 0 carries weight zero but is still a draw from q, so it
      # belongs in the denominator
      nz <- nz + a$rejected_zero_weights
      w  <- exp(eval_logf(cfg$p, a$trees, model = cfg$mb, link = cfg$link,
                          rho = 1)$logf - a$logg)
      nm <- vapply(a$trees, function(d) sum(is_mis_row(d)), 0)
      wle <- ifelse(nm <= 1, w, 0)
      N <- N + length(w)
      tot1 <- tot1 + sum(w[nm == 1]); ss1 <- ss1 + sum(w[nm == 1]^2)
      totle <- totle + sum(wle); ssle <- ssle + sum(wle^2)
      w2 <- w2 + sum(w[nm >= 2])
    }
    M <- N + nz
    e1  <- tot1 / M;  s1  <- sqrt(ss1 / M - e1^2) / sqrt(M)
    ele <- totle / M; sle <- sqrt(ssle / M - ele^2) / sqrt(M)
    expect_lt(abs(e1 - l1) / s1, 4,
              label = sprintf("|IS - L1| / se over the one-missing draws, %s", cfg$name))
    expect_lt(abs(ele - l0 - l1) / sle, 4,
              label = sprintf("|IS - (L0 + L1)| / se over the draws with at most one, %s",
                              cfg$name))
    # the truncation, measured rather than assumed: the draws left out
    expect_gt(w2 / (tot1 + totle), 0)       # there are some
    expect_lt(w2 / (tot1 + totle), 0.20)    # and they are a stated fraction
  }
})


# --------------------------------------------------------------------------- #
#  Gate 4: the D model charges what the model it extends charges               #
# --------------------------------------------------------------------------- #

test_that("beta_D = gamma_D = beta_M = gamma_M = 0 gives the cr / dd log q", {
  skip_on_cran()
  set.seed(5)
  phy  <- ape::rphylo(12L, 0.8, 0.2)
  brts <- emphasis:::.extract_brts(phy)
  pts  <- emphasis:::.pts(brts)
  trees <- list()
  for (spec in list(list(c(0L, 0L, 0L), c(0.9, 0, 0, 0, 0.35, 0, 0, 0)),
                    list(c(0L, 0L, 1L), c(0.9, 0, 0, 0.1, 0.35, 0, 0, 0.03)))) {
    a <- augment_trees(as.numeric(brts), spec[[2L]], 8L, 20000L, 300L, 1e6, 1L,
                       model = spec[[1L]], link = 0L, rho = 1, parent_tip_start = pts)
    trees <- c(trees, a$trees)
  }
  skip_if(length(trees) == 0L, "augmentation drew no tree")

  for (link in c(0L, 1L, 2L)) {
    for (rho in c(1, 0.8)) {
      for (g in list(list(c(0L, 0L, 0L), c(0.90,  0.00, 0, 0, 0.35, 0, 0, 0)),
                     list(c(1L, 0L, 0L), c(0.90, -0.02, 0, 0, 0.35, 0, 0, 0)))) {
        base <- eval_logf(g[[2L]], trees, model = g[[1L]], link = link, rho = rho)$logg
        dmod <- eval_logf(g[[2L]], trees, model = c(g[[1L]][1L], 0L, 1L),
                          link = link, rho = rho)$logg
        expect_true(all(is.finite(base)))
        expect_identical(dmod, base)
      }
    }
  }
})

test_that("log q does not read the D coefficients at all", {
  skip_on_cran()
  # The proposal is D-free by construction, so moving beta_D or gamma_D may not
  # move log q by a bit -- under any model flag, any link, any rho.  That is
  # what makes the reduction above exact rather than approximate.
  set.seed(23)
  phy  <- ape::rphylo(9L, 0.9, 0.3)
  brts <- emphasis:::.extract_brts(phy)
  pts  <- emphasis:::.pts(brts)
  a <- augment_trees(as.numeric(brts), c(1.0, 0, 0, 0.12, 0.45, 0, 0, 0.05),
                     10L, 20000L, 300L, 1e6, 1L, model = c(0L, 0L, 1L),
                     link = 0L, rho = 1, parent_tip_start = pts)
  skip_if(length(a$trees) == 0L, "augmentation drew no tree")
  p0 <- c(0.9, -0.02, 0,  0.00, 0.35, 0, 0,  0.00)
  p1 <- c(0.9, -0.02, 0,  0.77, 0.35, 0, 0, -0.41)
  for (link in c(0L, 1L, 2L)) {
    for (rho in c(1, 0.8)) {
      for (mb in list(c(0L, 0L, 0L), c(1L, 0L, 0L), c(0L, 0L, 1L), c(1L, 0L, 1L))) {
        expect_identical(
          eval_logf(p0, a$trees, model = mb, link = link, rho = rho)$logg,
          eval_logf(p1, a$trees, model = mb, link = link, rho = rho)$logg)
      }
    }
  }
})


# --------------------------------------------------------------------------- #
#  The envelope now dominates, on every model                                  #
# --------------------------------------------------------------------------- #

test_that("nh falls within a segment, so the start-of-segment envelope dominates", {
  skip_on_cran()
  set.seed(5)
  phy  <- ape::rphylo(14L, 0.8, 0.2)
  brts <- emphasis:::.extract_brts(phy)
  pts  <- emphasis:::.pts(brts)
  p    <- c(0.9, -0.01, 0, -0.50, 0.35, 0, 0, 0.30)
  thinning_envelope_violations(reset = TRUE)
  a <- augment_trees(as.numeric(brts), p, 200L, 200000L, 400L, 1e6, 1L,
                     model = c(1L, 0L, 1L), link = 0L, rho = 1, parent_tip_start = pts)
  expect_equal(a$envelope_violations, 0)
  expect_equal(thinning_envelope_violations(reset = TRUE), 0)

  # scanned directly, rather than through the sampler's own rejections: nh is
  # non-increasing on every segment of every tree it drew
  segs <- 0L; worst <- 0
  for (df in a$trees) {
    for (i in seq_len(nrow(df))) {
      prev <- if (i == 1L) 0 else df$brts[i - 1L]
      if (df$brts[i] <= prev) next
      # the half-open segment (prev, brts]: at t = prev, lower_bound picks the
      # node at prev, i.e. the segment before this one
      g  <- seq(prev, df$brts[i], length.out = 65L)[-1L]
      nh <- eval_nh_rate(p, df, g, model = c(1L, 0L, 1L), link = 0L, rho = 1)$nh
      segs <- segs + 1L
      worst <- max(worst, max(diff(nh)) / max(1e-12, max(nh)))
    }
  }
  expect_gt(segs, 1000L)
  expect_lt(worst, 1e-12)
})

# ---------------------------------------------------------------------------
# The proposal reads M at the segment start, so P must extrapolate to the same
# value whether or not a later birth has split the segment.  Deterministic:
# no draws, no tolerance on a Monte Carlo quantity.
# ---------------------------------------------------------------------------

test_that("the proposal's P stays inside the bounds a pendant PD must satisfy", {
  skip_on_cran()
  # P(t) = sum over alive lineages of (t - tip_start) is bounded below by 0 and
  # above by N*t (attained only when every alive lineage dates from the crown).
  # The proposal reads P by extrapolating node.pd along the segment, so if that
  # extrapolation used the wrong anchor -- the defect this commit removes -- the
  # bounds break before any distributional test would notice.
  set.seed(9)
  phy <- ape::rphylo(10, 0.7, 0.25)
  b   <- emphasis:::.extract_brts(phy)
  pts <- emphasis:::.pts(b)
  a <- augment_trees(b, c(0.7, 0, 0, 0, 0.25, 0, 0, 0), sample_size = 40L,
                     maxN = 4000L, max_missing = 1e4, max_lambda = 1e6,
                     num_threads = 1L, model = c(0L, 0L, 0L), link = 0L,
                     rho = 1, parent_tip_start = pts)
  worst_lo <- 0; worst_hi <- 0; checked <- 0L
  for (tr in a$trees) {
    tr <- tr[order(tr$brts), ]
    t0 <- 0
    for (k in seq_len(nrow(tr))) {
      # probe inside the segment this node governs
      for (u in c(0.1, 0.5, 0.9)) {
        t <- t0 + u * (tr$brts[k] - t0)
        if (t <= t0) next
        P <- tr$pd[k] + tr$n[k] * (t - tr$brts[k])
        worst_lo <- min(worst_lo, P)
        worst_hi <- max(worst_hi, P - tr$n[k] * t)
        checked <- checked + 1L
      }
      t0 <- tr$brts[k]
    }
  }
  expect_gt(checked, 500L)
  expect_gte(worst_lo, -1e-9)   # P >= 0
  expect_lte(worst_hi,  1e-9)   # P <= N*t
})

test_that("an M-dependent model without the topology is refused, not mis-scored", {
  brts <- c(3, 1.4, 0.6)
  p <- c(0.4, 0, 0.05, 0, 0.1, 0, 0.02, 0)
  # M active (slot 2 of model), no parent_tip_start: log q would not be the
  # density the sampler drew from, so the call must error rather than return.
  expect_error(
    augment_trees(brts, p, sample_size = 1L, maxN = 50L, max_missing = 1e4,
                  max_lambda = 1e6, num_threads = 1L, model = c(0L, 1L, 0L),
                  link = 0L, rho = 1),
    "topology")
  # the same model with the topology supplied runs
  expect_silent({
    a <- augment_trees(brts, p, sample_size = 1L, maxN = 50L, max_missing = 1e4,
                       max_lambda = 1e6, num_threads = 1L, model = c(0L, 1L, 0L),
                       link = 0L, rho = 1, parent_tip_start = c(0, 0, -1))
  })
  # models without M are unaffected by the guard
  expect_silent({
    a <- augment_trees(brts, p, sample_size = 1L, maxN = 50L, max_missing = 1e4,
                       max_lambda = 1e6, num_threads = 1L, model = c(1L, 0L, 1L),
                       link = 0L, rho = 1)
  })
})

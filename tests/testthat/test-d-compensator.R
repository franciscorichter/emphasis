# The D-model compensator under the linear link, and the pendant PD the
# thinning sampler reads.
#
# Two things are held here:
#
#   * Model::loglik's ep_linear branch integrates its own model.  Under the
#     linear link the per-lineage rate on a segment is max(0, c_s + b*t), a line
#     clipped at zero, and the compensator is the sum of those areas — not
#     dt * n * (rate at the segment's end node), which is what every link but
#     exponential and gaussian used to get.
#
#   * Model::nh_rate reads P off the node that governs the candidate's segment,
#     P(t) = node.pd + node.n * (t - node.brts), so the density the sampler draws
#     from and the density the scorer charges for are built on the same P.
#
# The bookkeeping the compensator uses for "the lineages alive on this segment"
# is the exponential branch's: the two crown lineages at tip_start 0, plus every
# lineage whose birth node lies at or before the segment start and which has not
# died before it ends.  A split adds the daughter and leaves the parent's
# tip_start where it was, which is not what the alive multiset behind node.pd
# does (H99, open); the reference below uses the same set, so what is tested
# here is the integral, not that bookkeeping.

T_TIP <- 10e10; T_EXT <- 0

# --------------------------------------------------------------------------- #
#  References, independent of the C++                                          #
# --------------------------------------------------------------------------- #

# The lineages the ep_linear / ep_exp branches count as alive on the segment
# that ends at node i.  Returns their tip_starts.
.alive_ts <- function(df, i) {
  n <- nrow(df)
  prev <- if (i == 1L) 0 else df$brts[i - 1L]
  j <- which(df$t_ext != T_EXT & seq_len(n) != n &
             df$brts <= prev & df$t_ext >= df$brts[i])
  c(0, 0, df$tip_start[j])
}

# Numerical integral of max(0, c + b*t) over [t1, t2].  The integrand is a line
# with one kink; integrate() is called on each side of it, where the integrand
# is smooth.  Independent of the closed form under test, which does the case
# split algebraically.
.int_relu_num <- function(cs, b, t1, t2) {
  vapply(cs, function(c0) {
    if (b == 0) return(max(0, c0) * (t2 - t1))
    root <- -c0 / b
    cuts <- sort(unique(c(t1, t2, if (root > t1 && root < t2) root)))
    s <- 0
    for (k in seq_len(length(cuts) - 1L)) {
      s <- s + stats::integrate(function(t) pmax(0, c0 + b * t),
                                cuts[k], cuts[k + 1L],
                                rel.tol = 1e-13, abs.tol = 1e-15)$value
    }
    s
  }, 0)
}

# The whole segment in one adaptive call, told nothing about the kinks: a second
# opinion on the splitting above.
.int_segment_blind <- function(p, df, i) {
  prev <- if (i == 1L) 0 else df$brts[i - 1L]
  if (df$brts[i] <= prev) return(0)
  M  <- if (df$n[i] > 0) df$pd[i] / df$n[i] else 0
  ts <- .alive_ts(df, i)
  tot <- function(t) vapply(t, function(u) {
    D <- (u - ts) - M
    sum(pmax(0, p[1] + p[2] * df$n[i] + p[3] * M + p[4] * D)) +
    sum(pmax(0, p[5] + p[6] * df$n[i] + p[7] * M + p[8] * D))
  }, 0)
  stats::integrate(tot, prev, df$brts[i], rel.tol = 1e-9,
                   subdivisions = 2000L)$value
}

# log f for a D model under the linear link: the event terms exactly as
# speciation_rate_ep / extinction_rate_ep compute them, and a compensator that
# is integrated numerically rather than in closed form.
.logf_ref <- function(p, df) {
  n <- nrow(df)
  M <- ifelse(df$n > 0, df$pd / df$n, 0)
  E <- ifelse(df$t_ext == T_EXT, df$brts - df$tip_start,
       ifelse(df$focal_tip_start >= 0, df$brts - df$focal_tip_start, M))
  D <- E - M
  lambda <- pmax(0, p[1] + p[2] * df$n + p[3] * M + p[4] * D)
  mu     <- pmax(0, p[5] + p[6] * df$n + p[7] * M + p[8] * D)

  ev <- 0
  for (i in seq_len(n)) {
    if (df$t_ext[i] == T_EXT) ev <- ev + log(max(mu[i], 1e-300))
    else if (i != n)          ev <- ev + log(lambda[i])
  }

  inte <- 0
  for (i in seq_len(n)) {
    prev <- if (i == 1L) 0 else df$brts[i - 1L]
    if (df$brts[i] <= prev) next
    ts <- .alive_ts(df, i)
    A_l <- p[1] + p[2] * df$n[i] + (p[3] - p[4]) * M[i] - p[4] * ts
    A_m <- p[5] + p[6] * df$n[i] + (p[7] - p[8]) * M[i] - p[8] * ts
    inte <- inte + sum(.int_relu_num(A_l, p[4], prev, df$brts[i])) +
                   sum(.int_relu_num(A_m, p[8], prev, df$brts[i]))
  }
  ev - inte
}

# The alive multiset node.pd is built on, replayed from the augmented tree's own
# columns: an extinction node carries the tip_start of the lineage that died and
# a split node the tip_start of the lineage that split, so the set before every
# event can be rebuilt without knowing how it was produced.
.alive_before <- function(df) {
  out   <- vector("list", nrow(df))
  alive <- c(0, 0)                       # the two crown lineages
  for (i in seq_len(nrow(df))) {
    out[[i]] <- alive
    if (df$t_ext[i] == T_EXT) {
      alive <- alive[-which.min(abs(alive - df$tip_start[i]))]
    } else if (i != nrow(df)) {
      want  <- max(df$focal_tip_start[i], 0)
      alive <- c(alive[-which.min(abs(alive - want))], df$brts[i], df$brts[i])
    }
  }
  out
}

# --------------------------------------------------------------------------- #
#  Hand-built trees with known tip starts                                      #
# --------------------------------------------------------------------------- #

# Crown at 0 with two lineages at tip_start 0; a missing lineage born at 1 off a
# crown lineage and dead at 4; an observed split at 5 of the lineage the missing
# one came from; the present at 10.  The alive multiset runs
#   (0,1] {0,0}  (1,4] {0,1,1}  (4,5] {0,1}  (5,10] {0,5,5}
# and pd = N*t - sum(ts) is read straight off it.
.hand_tree <- function() {
  data.frame(
    brts            = c(1,   4,     5,     10),
    n               = c(2,   3,     2,     3),
    t_ext           = c(4,   T_EXT, T_TIP, T_TIP),
    pd              = c(2,   10,    9,     20),
    tip_start       = c(1,   1,     5,     10),
    focal_tip_start = c(0,   1,     1,     -1),
    id              = c(4L,  4L,    0L,    1L),
    parent_id       = c(-1L, -1L,   -1L,   -1L))
}

# A second shape: two missing lineages alive at once, so a segment carries three
# distinct tip starts.  Births at 1 (dead at 8) and 2 (dead at 6), an observed
# split at 4 of the lineage born at 1, the present at 9.  The alive multiset:
#   (0,1] {0,0}      (1,2] {0,1,1}    (2,4] {1,1,2,2}
#   (4,6] {1,2,2,4,4}  (6,8] {1,2,4,4}  (8,9] {2,4,4}
.hand_tree2 <- function() {
  data.frame(
    brts            = c(1,   2,   4,     6,     8,     9),
    n               = c(2,   3,   4,     5,     4,     3),
    t_ext           = c(8,   6,   T_TIP, T_EXT, T_EXT, T_TIP),
    pd              = c(2,   4,   10,    17,    21,    17),
    tip_start       = c(1,   2,   4,     2,     1,     9),
    focal_tip_start = c(0,   0,   1,     2,     1,     -1),
    id              = c(7L,  8L,  0L,    8L,    7L,    1L),
    parent_id       = c(-1L, -1L, -1L,   -1L,   -1L,   -1L))
}

# theta grids.  A and D put the kink inside a segment with b > 0, B and C with
# b < 0, E leaves every rate positive throughout (no clipping) and F adds an N
# term.  All six keep every event rate strictly positive, so log f is finite.
.thetas <- list(
  c(1.0,  0,     0, 0.50, 0.35,  0,     0, 0.20),
  c(1.0,  0,     0, -0.50, 0.35, 0,     0, -0.20),
  c(2.0,  0,     0, -0.30, 0.60, 0,     0, -0.45),
  c(1.5,  0,     0, 0.35, 0.25,  0,     0, 0.30),
  c(2.0,  0,     0, 0.05, 0.80,  0,     0, 0.02),
  c(1.2, -0.05,  0, 0.15, 0.40, -0.01,  0, 0.08)
)

.cpp_logf <- function(p, df, model = c(0L, 0L, 1L), link = 0L, rho = 1) {
  eval_logf(p, list(df), model = as.integer(model), link = as.integer(link),
            rho = rho)$logf
}


test_that("the ep_linear compensator is the exact integral of its own rate", {
  skip_on_cran()
  worst <- 0
  for (df in list(.hand_tree(), .hand_tree2())) {
    for (p in .thetas) {
      got <- .cpp_logf(p, df)
      ref <- .logf_ref(p, df)
      expect_true(is.finite(got))
      expect_true(is.finite(ref))
      expect_lt(abs(got - ref), 1e-10)
      worst <- max(worst, abs(got - ref))
    }
  }
  expect_lt(worst, 1e-10)
})


test_that("clipping really bites on these trees, in both directions", {
  skip_on_cran()
  # A theta only tests the kink if some lineage's rate changes sign inside a
  # segment.  Check that directly, so the test above cannot pass vacuously.
  crossings <- function(df, p) {
    n <- nrow(df); k <- 0
    for (i in seq_len(n)) {
      prev <- if (i == 1L) 0 else df$brts[i - 1L]
      if (df$brts[i] <= prev) next
      M  <- df$pd[i] / df$n[i]
      ts <- .alive_ts(df, i)
      for (sl in list(c(1, 2, 3, 4), c(5, 6, 7, 8))) {
        A <- p[sl[1]] + p[sl[2]] * df$n[i] + (p[sl[3]] - p[sl[4]]) * M - p[sl[4]] * ts
        b <- p[sl[4]]
        if (b != 0) {
          root <- -A / b
          k <- k + sum(root > prev & root < df$brts[i])
        }
      }
    }
    k
  }
  up   <- sum(vapply(list(.hand_tree(), .hand_tree2()),
                     function(d) crossings(d, .thetas[[1L]]), 0))
  down <- sum(vapply(list(.hand_tree(), .hand_tree2()),
                     function(d) crossings(d, .thetas[[2L]]), 0))
  none <- sum(vapply(list(.hand_tree(), .hand_tree2()),
                     function(d) crossings(d, .thetas[[5L]]), 0))
  expect_gt(up, 0)       # b > 0, the rate rises through zero
  expect_gt(down, 0)     # b < 0, it falls through zero
  expect_identical(none, 0)   # and the unclipped case is covered too
})


test_that("the kink split agrees with blind adaptive quadrature", {
  skip_on_cran()
  # Same compensator, integrated without being told where the kinks are.
  for (df in list(.hand_tree(), .hand_tree2())) {
    for (p in .thetas) {
      blind <- sum(vapply(seq_len(nrow(df)), function(i) .int_segment_blind(p, df, i), 0))
      split <- 0
      M <- df$pd / df$n
      for (i in seq_len(nrow(df))) {
        prev <- if (i == 1L) 0 else df$brts[i - 1L]
        if (df$brts[i] <= prev) next
        ts <- .alive_ts(df, i)
        split <- split +
          sum(.int_relu_num(p[1] + p[2] * df$n[i] + (p[3] - p[4]) * M[i] - p[4] * ts,
                            p[4], prev, df$brts[i])) +
          sum(.int_relu_num(p[5] + p[6] * df$n[i] + (p[7] - p[8]) * M[i] - p[8] * ts,
                            p[8], prev, df$brts[i]))
      }
      expect_lt(abs(split - blind), 1e-6 * max(1, abs(blind)))
    }
  }
})


test_that("the ep_linear compensator is exact on augmented trees too", {
  skip_on_cran()
  set.seed(5)
  phy  <- ape::rphylo(12L, 0.8, 0.2)
  brts <- emphasis:::.extract_brts(phy)
  aug <- augment_trees(as.numeric(brts), c(0.9, 0, 0, 0.10, 0.35, 0, 0, 0.03),
                       10L, 20000L, 300L, 1e6, 1L,
                       model = c(0L, 0L, 1L), link = 0L, rho = 1,
                       parent_tip_start = emphasis:::.pts(brts))
  skip_if(length(aug$trees) == 0L, "augmentation drew no tree")

  worst <- 0
  for (df in aug$trees) {
    for (p in .thetas[c(3L, 5L, 6L)]) {
      got <- .cpp_logf(p, df)
      ref <- .logf_ref(p, df)
      if (!is.finite(got) && !is.finite(ref)) next
      worst <- max(worst, abs(got - ref))
    }
  }
  expect_lt(worst, 1e-10)
})


# --------------------------------------------------------------------------- #
#  Gate: the D model reduces to the model it extends                           #
# --------------------------------------------------------------------------- #

test_that("beta_M = beta_D = gamma_M = gamma_D = 0 reproduces cr and dd", {
  skip_on_cran()
  set.seed(5)
  phy  <- ape::rphylo(12L, 0.8, 0.2)
  brts <- emphasis:::.extract_brts(phy)
  trees <- list()
  for (spec in list(list(c(0L, 0L, 0L), c(0.9, 0, 0, 0, 0.35, 0, 0, 0)),
                    list(c(0L, 0L, 1L), c(0.9, 0, 0, 0.1, 0.35, 0, 0, 0.03)))) {
    a <- augment_trees(as.numeric(brts), spec[[2L]], 8L, 20000L, 300L, 1e6, 1L,
                       model = spec[[1L]], link = 0L, rho = 1,
                       parent_tip_start = emphasis:::.pts(brts))
    trees <- c(trees, a$trees)
  }
  skip_if(length(trees) == 0L, "augmentation drew no tree")
  trees <- c(trees, list(.hand_tree(), .hand_tree2()))

  # The D coefficients are the only thing switched off, so the D branch has to
  # collapse onto dt * n * (lambda + mu) — the same number the model it extends
  # computes without ever entering the branch.
  #
  # Links 0 and 1 only.  The gaussian branch counts the two crown lineages out
  # of its alive set (H99, untouched here), so it does not reduce; see the
  # commit message.
  for (link in c(0L, 1L)) {
    grid <- if (link == 0L)
      list(list(c(0L, 0L, 0L), c(0.90,  0.00, 0, 0, 0.35, 0, 0, 0)),
           list(c(1L, 0L, 0L), c(0.90, -0.02, 0, 0, 0.35, 0, 0, 0)),
           list(c(1L, 0L, 0L), c(1.30, -0.06, 0, 0, 0.10, 0, 0, 0)))
    else
      list(list(c(0L, 0L, 0L), c(-0.10,  0.00, 0, 0, -1.20, 0, 0, 0)),
           list(c(1L, 0L, 0L), c(-0.10, -0.01, 0, 0, -1.20, 0, 0, 0)))
    for (g in grid) {
      base <- eval_logf(g[[2L]], trees, model = g[[1L]], link = link, rho = 1)$logf
      dmod <- eval_logf(g[[2L]], trees, model = c(g[[1L]][1L], 0L, 1L),
                        link = link, rho = 1)$logf
      expect_true(all(is.finite(base)))
      expect_lt(max(abs(dmod - base)), 1e-10)
    }
  }
})


test_that("cr and dd never enter the D branch", {
  skip_on_cran()
  # model_bin[2] == 0 is the gate; the D coefficients are ignored entirely, so
  # moving them may not move a cr or dd log f by a bit.
  df <- .hand_tree2()
  for (link in c(0L, 1L, 2L)) {
    for (mb in list(c(0L, 0L, 0L), c(1L, 0L, 0L))) {
      p0 <- c(0.9, -0.02, 0, 0.00, 0.35, 0, 0, 0.00)
      p1 <- c(0.9, -0.02, 0, 0.77, 0.35, 0, 0, -0.41)
      expect_identical(eval_logf(p0, list(df), model = mb, link = link, rho = 1),
                       eval_logf(p1, list(df), model = mb, link = link, rho = 1))
    }
  }
})


# --------------------------------------------------------------------------- #
#  Gate: the P the sampler reads                                               #
# --------------------------------------------------------------------------- #

test_that("nh_rate reads P off the node that governs the candidate's segment", {
  skip_on_cran()
  set.seed(5)
  phy  <- ape::rphylo(14L, 0.8, 0.2)
  brts <- emphasis:::.extract_brts(phy)
  p    <- c(0.9, 0, 0, 0.10, 0.35, 0, 0, 0.03)
  aug <- augment_trees(as.numeric(brts), p, 12L, 20000L, 300L, 1e6, 1L,
                       model = c(0L, 0L, 1L), link = 0L, rho = 1,
                       parent_tip_start = emphasis:::.pts(brts))
  skip_if(length(aug$trees) == 0L, "augmentation drew no tree")

  worst_ref <- 0; worst_lin <- 0; probes <- 0L
  for (df in aug$trees) {
    ab <- .alive_before(df)
    # node.n is the size of that alive set: the count on the segment that ends
    # at the node, which is the N the extrapolation uses as a slope.
    expect_equal(df$n, vapply(ab, length, 0L), tolerance = 0)
    for (i in seq_len(nrow(df))) {
      prev <- if (i == 1L) 0 else df$brts[i - 1L]
      if (df$brts[i] <= prev) next
      t <- seq(prev, df$brts[i], length.out = 9L)[-c(1L, 9L)]
      got <- eval_nh_rate(p, df, t, model = c(0L, 0L, 1L), link = 0L, rho = 1)$pd
      worst_ref <- max(worst_ref, max(abs(got - (length(ab[[i]]) * t - sum(ab[[i]])))))
      worst_lin <- max(worst_lin, max(abs(got - (df$pd[i] + df$n[i] * (t - df$brts[i])))))
      probes <- probes + length(t)
    }
  }
  expect_gt(probes, 500L)
  expect_lt(worst_lin, 1e-10)   # P(t) = node.pd + node.n * (t - node.brts)
  expect_lt(worst_ref, 1e-10)   # and that is the pendant PD of a real alive set
})


test_that("the pendant PD the sweep stores is the one a replay produces", {
  skip_on_cran()
  set.seed(19)
  phy  <- ape::rphylo(13L, 0.9, 0.3)
  brts <- emphasis:::.extract_brts(phy)
  aug <- augment_trees(as.numeric(brts), c(1.0, 0, 0, 0.12, 0.45, 0, 0, 0.05),
                       12L, 20000L, 300L, 1e6, 1L,
                       model = c(0L, 0L, 1L), link = 0L, rho = 1,
                       parent_tip_start = emphasis:::.pts(brts))
  skip_if(length(aug$trees) == 0L, "augmentation drew no tree")
  expect_true(any(vapply(aug$trees, nrow, 0L) > length(brts)))   # lineages added

  for (df in aug$trees) {
    ab  <- .alive_before(df)
    ref <- vapply(seq_len(nrow(df)),
                  function(i) length(ab[[i]]) * df$brts[i] - sum(ab[[i]]), 0)
    expect_equal(df$pd, ref, tolerance = 1e-10)
    expect_true(all(df$pd >= -1e-10))
    expect_true(all(df$pd <= df$n * df$brts + 1e-10))
  }
})

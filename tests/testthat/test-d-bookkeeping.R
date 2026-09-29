# The D compensator on a tree whose pendant ages can be read off the topology.
#
# On a fully observed mu = 0 tree every lineage's pendant start is the start of
# the edge it is on, and it restarts when the lineage splits.  The compensator
# has to integrate the rate over that alive set -- the same one the event terms
# and M are computed from -- under all three links.  A reference that never
# resets a start (each lineage keeps the time it was born) is what the check
# must reject: the two differ by an amount proportional to beta_D.

.edge_starts <- function(phy, t) {
  h  <- ape::node.depth.edgelength(phy)
  st <- h[phy$edge[, 1L]]
  en <- h[phy$edge[, 2L]]
  st[st < t & en >= t]
}

.rate <- function(link, intercept, arg) {
  if (link == 0L) pmax(0, intercept + arg)
  else if (link == 1L) exp(intercept + arg)
  else intercept * exp(-0.5 * (arg - 1)^2)
}

# log f for the model c(0, 0, 1) with beta_N = beta_M = gamma_* = 0 (so mu = 0
# and no extinction term), integrating the compensator numerically over the
# given alive starts per segment.
.ref_logf <- function(phy, df, pts, b0, bD, link, starts_of) {
  n <- nrow(df); ev <- 0; inte <- 0; prev <- 0
  for (i in seq_len(n)) {
    t  <- df$brts[i]
    st <- .edge_starts(phy, t)
    M  <- sum(t - st) / length(st)
    ts <- starts_of(i, st)
    f  <- function(u) vapply(u, function(v) sum(.rate(link, b0, bD * ((v - ts) - M))), 0)
    inte <- inte + stats::integrate(f, prev, t, rel.tol = 1e-10, subdivisions = 500L)$value
    if (i < n) ev <- ev + log(.rate(link, b0, bD * ((t - pts[i]) - M)))
    prev <- t
  }
  ev - inte
}

.observed_frame <- function(phy) {
  brts <- emphasis:::.extract_brts(phy)
  pts  <- emphasis:::.pts(brts)
  a <- augment_trees(as.numeric(brts), c(0.3, 0, 0, 0, 0, 0, 0, 0),
                     1L, 500L, 100L, 1e6, 1L,
                     model = c(0L, 0L, 0L), link = 0L, rho = 1,
                     parent_tip_start = pts)
  list(df = a$trees[[1L]], pts = pts)
}

test_that("the D compensator integrates over pendant ages that restart at a split", {
  skip_on_cran()
  worst <- 0; gap <- 0
  for (seed in c(4L, 11L, 23L)) {
    set.seed(seed)
    phy <- ape::rphylo(12L, 0.6, 0)
    fr  <- .observed_frame(phy); df <- fr$df; pts <- fr$pts
    expect_equal(nrow(df), length(emphasis:::.extract_brts(phy)))
    for (link in 0:2) {
      b0 <- if (link == 1L) log(0.4) else 0.4
      for (bD in c(-0.05, 0.05, 0.15)) {
        got <- eval_logf(c(b0, 0, 0, bD, if (link == 1L) -30 else 0, 0, 0, 0), list(df),
                         model = c(0L, 0L, 1L), link = link, rho = 1)$logf
        # under the exponential link mu = exp(-30) is not zero; drop its
        # compensator, which is (T * n-weighted) exp(-30) ~ 1e-12
        reset <- .ref_logf(phy, df, pts, b0, bD, link, function(i, st) st)
        birth <- .ref_logf(phy, df, pts, b0, bD, link, function(i, st) {
          n <- nrow(df); prev <- if (i == 1L) 0 else df$brts[i - 1L]
          j <- which(df$t_ext != 0 & seq_len(n) != n & df$brts <= prev & df$t_ext >= df$brts[i])
          c(0, 0, df$tip_start[j])
        })
        worst <- max(worst, abs(got - reset))
        gap   <- max(gap, abs(reset - birth))
      }
    }
  }
  expect_lt(worst, 1e-6)
  # and the check is not vacuous: the two bookkeepings disagree on these trees
  expect_gt(gap, 0.1)
})

test_that("the multiset the compensator reads reproduces node.pd on augmented trees", {
  skip_on_cran()
  set.seed(5)
  phy  <- ape::rphylo(12L, 0.8, 0.2)
  brts <- emphasis:::.extract_brts(phy)
  aug <- augment_trees(as.numeric(brts), c(0.9, 0, 0, 0.10, 0.35, 0, 0, 0.03),
                       10L, 20000L, 300L, 1e6, 1L,
                       model = c(0L, 0L, 1L), link = 0L, rho = 1,
                       parent_tip_start = emphasis:::.pts(brts))
  skip_if(length(aug$trees) == 0L, "augmentation drew no tree")
  replay <- function(df) {
    alive <- c(0, 0); out <- vector("list", nrow(df))
    for (i in seq_len(nrow(df))) {
      out[[i]] <- alive
      if (df$t_ext[i] == 0) alive <- alive[-which.min(abs(alive - df$tip_start[i]))]
      else if (i != nrow(df)) {
        want  <- max(df$focal_tip_start[i], 0)
        alive <- c(alive[-which.min(abs(alive - want))], df$brts[i], df$brts[i])
      }
    }
    out
  }
  for (df in aug$trees) {
    al <- replay(df)
    expect_equal(vapply(al, length, 1L), as.integer(df$n))
    expect_equal(vapply(seq_along(al), function(i) sum(df$brts[i] - al[[i]]), 0), df$pd, tolerance = 1e-8)
  }
})

# The D compensator on a tree whose pendant ages can be read off the topology.
#
# On a fully observed mu = 0 tree every lineage's pendant start is the start of
# the edge it is on, and it restarts when the lineage splits.  Between events
# every pendant age E_s(u) = u - ts_s and their mean M(u) grow at unit rate, so
# D_s(u) = E_s(u) - M(u) is constant on a segment; the forward simulator holds
# it there, and the compensator has to integrate the rate at that D over the
# alive set the event terms and M are computed from, under both links.
# Two references are what the check must reject: one that never resets a start
# (each lineage keeps the time it was born), and one that holds M at its value
# at the segment's end and lets D = (u - ts_s) - M grow within the segment (the
# convention the likelihood had before).  Both differ from the model by an
# amount proportional to beta_D.

.edge_starts <- function(phy, t) {
  h  <- ape::node.depth.edgelength(phy)
  st <- h[phy$edge[, 1L]]
  en <- h[phy$edge[, 2L]]
  st[st < t & en >= t]
}

.rate <- function(link, intercept, arg) {
  if (link == 0L) pmax(0, intercept + arg)
  else exp(intercept + arg)
}

# log f for the model c(0, 0, 1) with beta_N = beta_M = gamma_* = 0 (so mu = 0
# and no extinction term), integrating the compensator numerically over the
# given alive starts per segment.  D_s(u) = (u - ts_s) - M(u) with M(u) the
# mean pendant age of the edges alive at u, read at every u; with
# m_held = TRUE, M is held at its value at the segment's end instead.
.ref_logf <- function(phy, df, pts, b0, bD, link, starts_of, m_held = FALSE) {
  n <- nrow(df); ev <- 0; inte <- 0; prev <- 0
  for (i in seq_len(n)) {
    t  <- df$brts[i]
    st <- .edge_starts(phy, t)
    M  <- sum(t - st) / length(st)
    Mu <- function(v) if (m_held) M else sum(v - st) / length(st)
    ts <- starts_of(i, st)
    f  <- function(u) vapply(u, function(v) sum(.rate(link, b0, bD * ((v - ts) - Mu(v)))), 0)
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

test_that("the D compensator holds D constant on a segment, over pendant ages that restart at a split", {
  skip_on_cran()
  worst <- 0; gap <- 0; gap_held <- 0
  for (seed in c(4L, 11L, 23L)) {
    set.seed(seed)
    phy <- ape::rphylo(12L, 0.6, 0)
    fr  <- .observed_frame(phy); df <- fr$df; pts <- fr$pts
    expect_equal(nrow(df), length(emphasis:::.extract_brts(phy)))
    for (link in 0:1) {
      b0 <- if (link == 1L) log(0.4) else 0.4
      for (bD in c(-0.05, 0.05, 0.15)) {
        got <- eval_logf(c(b0, 0, 0, bD, if (link == 1L) -30 else 0, 0, 0, 0), list(df),
                         model = c(0L, 0L, 1L), link = link, rho = 1)$logf
        # under the exponential link mu = exp(-30) is not zero; drop its
        # compensator, which is (T * n-weighted) exp(-30) ~ 1e-12
        reset <- .ref_logf(phy, df, pts, b0, bD, link, function(i, st) st)
        held  <- .ref_logf(phy, df, pts, b0, bD, link, function(i, st) st, m_held = TRUE)
        birth <- .ref_logf(phy, df, pts, b0, bD, link, function(i, st) {
          n <- nrow(df); prev <- if (i == 1L) 0 else df$brts[i - 1L]
          j <- which(df$t_ext != 0 & seq_len(n) != n & df$brts <= prev & df$t_ext >= df$brts[i])
          c(0, 0, df$tip_start[j])
        })
        worst    <- max(worst, abs(got - reset))
        gap      <- max(gap, abs(reset - birth))
        gap_held <- max(gap_held, abs(reset - held))
      }
    }
  }
  expect_lt(worst, 1e-6)
  # and the check is not vacuous: neither the never-reset bookkeeping nor the
  # M-held-at-the-end convention agrees with the model on these trees
  expect_gt(gap, 0.1)
  expect_gt(gap_held, 0.1)
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

# --------------------------------------------------------------------------- #
#  A complete simulated tree: the likelihood is the simulator's               #
# --------------------------------------------------------------------------- #

# The events of a DDD L-table (crown age T) in forward time: births of every
# row but the two crown rows, deaths of every extinct row.
.l_events <- function(L, T) {
  b <- T - L[, 1L]; b[1:2] <- 0
  d <- ifelse(L[, 4L] == -1, Inf, T - L[, 4L])
  ev <- rbind(data.frame(t = b[-(1:2)], type = 1L, row = 3:nrow(L)),
              data.frame(t = d[is.finite(d)], type = rep(2L, sum(is.finite(d))), row = which(is.finite(d))))
  list(ev = ev[order(ev$t), ], prow = c(NA, NA, match(L[-(1:2), 2L], L[, 3L])), d = d)
}

# The complete tree as the estimator's frame, every lineage on record and
# nothing hidden: a birth node per lineage (t_ext its death time, or the tip
# marker when it is alive at T), an extinction node per death, the present.
# n, pd, tip_start and focal_tip_start are replayed here from the L-table.
.complete_frame <- function(L, T) {
  e <- .l_events(L, T); ev <- e$ev; prow <- e$prow
  ids <- c(-2L, -3L, seq_len(nrow(L) - 2L) - 1L)
  ts <- c(0, 0, rep(NA_real_, nrow(L) - 2L)); alive <- c(TRUE, TRUE, rep(FALSE, nrow(L) - 2L))
  out <- vector("list", nrow(ev) + 1L)
  for (k in seq_len(nrow(ev))) {
    t <- ev$t[k]; r <- ev$row[k]; a <- which(alive)
    n <- length(a); pd <- sum(t - ts[a])
    if (ev$type[k] == 1L) {
      p <- prow[r]
      out[[k]] <- data.frame(brts = t, n = n, t_ext = if (is.finite(e$d[r])) e$d[r] else 10e10,
                             pd = pd, tip_start = t, focal_tip_start = ts[p], clade = 0L,
                             id = ids[r], parent_id = ids[p])
      ts[p] <- t; ts[r] <- t; alive[r] <- TRUE
    } else {
      out[[k]] <- data.frame(brts = t, n = n, t_ext = 0, pd = pd, tip_start = ts[r],
                             focal_tip_start = ts[r], clade = 0L, id = ids[r], parent_id = -1L)
      alive[r] <- FALSE
    }
  }
  a <- which(alive)
  out[[nrow(ev) + 1L]] <- data.frame(brts = T, n = length(a), t_ext = 10e10, pd = sum(T - ts[a]),
                                     tip_start = T, focal_tip_start = -1, clade = 0L,
                                     id = -1L, parent_id = -1L)
  do.call(rbind, out)
}

# The complete-data log-likelihood of an N + D model, straight off the L-table
# and independent of the C++: on every segment between events each alive
# lineage's D is the constant mean(ts) - ts_s (the mean pendant start of the
# alive lineages minus its own), its rate is constant, the compensator is the
# segment's length times the summed rates, and each event is charged the rate
# of its own lineage just before it.  p8 = (b0, bN, bM, bD, g0, gN, gM, gD).
.ref_complete_logf <- function(L, T, p8, link) {
  e <- .l_events(L, T); ev <- e$ev; prow <- e$prow
  ts <- c(0, 0, rep(NA_real_, nrow(L) - 2L)); alive <- c(TRUE, TRUE, rep(FALSE, nrow(L) - 2L))
  rate <- function(c0, cN, cD, N, D) {
    eta <- c0 + cN * N + cD * D
    if (link == 1L) exp(eta) else pmax(0, eta)
  }
  ll <- 0; prev <- 0
  for (k in seq_len(nrow(ev) + 1L)) {
    t <- if (k <= nrow(ev)) ev$t[k] else T
    a <- which(alive); N <- length(a); D <- mean(ts[a]) - ts[a]
    lam <- rate(p8[1], p8[2], p8[4], N, D); mu <- rate(p8[5], p8[6], p8[8], N, D)
    ll <- ll - (t - prev) * sum(lam + mu)
    if (k > nrow(ev)) break
    r <- ev$row[k]
    if (ev$type[k] == 1L) {
      p <- prow[r]
      ll <- ll + log(lam[a == p])
      ts[p] <- t; ts[r] <- t; alive[r] <- TRUE
    } else {
      ll <- ll + log(mu[a == r])
      alive[r] <- FALSE
    }
    prev <- t
  }
  ll
}

test_that("on a complete simulated tree the d and nd likelihood is the constant-D one", {
  skip_on_cran()
  T <- 6
  cases <- list(
    list(model = "d",  link = "exponential", li = 1L, pars = c(log(0.45), 0.30, log(0.12), -0.25)),
    list(model = "nd", link = "exponential", li = 1L, pars = c(log(0.6), -0.01, -0.30, log(0.12), 0, 0.20)),
    list(model = "d",  link = "linear",      li = 0L, pars = c(0.45, 0.08, 0.12, -0.03)),
    list(model = "nd", link = "linear",      li = 0L, pars = c(0.55, -0.004, 0.08, 0.12, 0, -0.03)))
  worst <- 0; n_trees <- 0L; n_ext <- 0L
  for (cs in cases) {
    mb <- emphasis:::.resolve_model(cs$model)
    p8 <- emphasis:::.expand_pars(cs$pars, mb)
    set.seed(31)
    for (rep in 1:3) {
      s <- simulate_tree(pars = cs$pars, max_t = T, model = cs$model, link = cs$link,
                         max_lin = 2000L, max_tries = 20L, useDDD = FALSE)
      if (s$status != "done") next
      df  <- .complete_frame(s$L, T)
      got <- eval_logf(p8, list(df), model = mb[1:3], link = cs$li, rho = 1)$logf
      ref <- .ref_complete_logf(s$L, T, p8, cs$li)
      expect_true(is.finite(got))
      worst <- max(worst, abs(got - ref) / max(1, abs(ref)))
      n_trees <- n_trees + 1L; n_ext <- n_ext + sum(df$t_ext == 0)
    }
  }
  expect_gt(n_trees, 8L)
  expect_gt(n_ext, 0L)          # extinct lineages are on these trees
  expect_lt(worst, 1e-9)
})

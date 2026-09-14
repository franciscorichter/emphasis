# The parent an augmented lineage is drawn from (audit finding H45).
#
# THE JOINT DENSITY.  The thinning sampler draws, in forward time,
#
#   birth times  an inhomogeneous Poisson process of intensity nh(t);
#   a lifetime   for the lineage born at t_k;
#   a parent     for it.
#
# The first two are held by test-proposal-density.R.  This file holds the third.
#
# What f scores is a complete labelled history: every event names the lineage it
# happened to, one lineage of every split carries on under the splitting
# lineage's name and the other is born.  An augmented tree does not record the
# names, so several histories can land on one augmented tree, and the density
# the sampler must draw from is the one over histories:
#
#   a birth off an OBSERVED lineage has TWO histories.  The split leaves an
#   observed branch and a missing one, and either of them may be the one that
#   carries on under the observed lineage's name.  Both give the same augmented
#   tree, the same tip starts and therefore the same f.
#
#   a birth off an AUGMENTED lineage has ONE.  Its two branches differ by their
#   death times, so swapping which of them carries on is a different
#   augmentation, already a separate point of the space -- not a second name for
#   this one.
#
# So a birth at t has 2 * tips + Ne labelled attachments, tips being the
# observed lineages alive (the two crown lineages included) and Ne the augmented
# ones, and drawing one of them uniformly is the -log(2 * tips + Ne) that
# Model::sampling_prob already charges (audit finding H5, refuted; the term is
# not changed here).  E_q[f/q] is then the integral of f over augmentations with
# every labelled attachment counted once, which is what gate 1 integrates.
#
# The defect was in the sampler, not the density: the candidate list was built
# from the nodes of the tree, and the two crown lineages carry no node, so they
# were never candidates -- n - 2 of the n lineages alive, with parent -1 and a
# lineage dropped from the L-table when the list was empty.  Attachments f gives
# positive density had proposal probability zero, so E_q[f/q] could not be the
# integral of f for any model that reads the parent, which is every D model.
#
# Gate 1  E_q[f/q] is the brute-force integral of f on a tree WITH an observed
#         split, where the attachments differ in f and the crown lineages are
#         among them: a 3-tip tree (one split, K = 4 then 6) and a balanced
#         4-tip tree (two splits, K = 4, 6, 8).
# Gate 2  the tree that comes back is one application of the sweep driven by the
#         caller's parent_tip_start -- the state the sampler carried is the
#         state the finished tree reports.  A tree with one observed branching
#         cannot see this: before the first split both crown lineages report a
#         tip start of 0 and every pass breaks the tie the same way.
# Gate 3  every augmented lineage has a parent on record and reaches the tas.
# Gate 4  the candidates are the lineages alive, and the sampler and the density
#         count the same attachments.
# Gate 7  the same integral with no topology at all, where the legacy
#         convention dates every observed lineage from the crown.  The parent
#         draw changed that path too: a birth before the first observed node
#         used to carry no parent, and so fell back to D = 0.

T_TIP <- 10e10; T_EXT <- 0; T_UNS <- 5e10
.is_mis <- function(te) !(te == T_EXT | te == T_TIP | te == T_UNS)

CROWN_A <- -2L   # the crown lineage the sweep splits at the first observed event
CROWN_B <- -3L


# --------------------------------------------------------------------------- #
#  The 3-tip tree, by hand                                                     #
# --------------------------------------------------------------------------- #
#
# Crown at 0, one observed split at t1, the present at TT.  Write A for the
# crown lineage that splits at t1, B for the other, C for the lineage born at
# t1.  An augmentation with one missing lineage is (t, d, parent), and the
# attachments are
#
#   t < t1   A, A, B, B            (K = 4)
#   t > t1   A, A, C, C, B, B      (K = 6),
#
# where A and C carry tip start t1 after the split and B still carries 0, so f
# takes two values in each stratum.  P = sum over the alive lineages of
# (t - their tip start) is written out at every node below; a split resets the
# tip start of the splitting lineage as well as of its daughter.

.row3 <- function(brts, n, t_ext, pd, tip_start, focal)
  data.frame(brts = brts, n = n, t_ext = t_ext, pd = pd, tip_start = tip_start,
             focal_tip_start = focal, id = 0L, parent_id = -1L)

.aug0_3 <- function(TT, t1) rbind(
  .row3(t1, 2, T_TIP, 2 * t1,          t1, 0),
  .row3(TT, 3, T_TIP, 3 * TT - 2 * t1, TT, -1))

# `mode` is "A" or "B" when t < t1, and the parent's tip start (t1 or 0) when
# t > t1.
.aug1_3 <- function(t, d, TT, t1, mode) {
  if (t < t1) {
    A <- identical(mode, "A")
    if (d < t1) {
      rbind(.row3(t,  2, d,     2 * t,          t,  0),
            .row3(d,  3, T_EXT, 3 * d - 2 * t,  t,  t),
            .row3(t1, 2, T_TIP, 2 * t1 - t,     t1, if (A) t else 0),
            .row3(TT, 3, T_TIP,
                  if (A) 3 * TT - 2 * t1 else 3 * TT - (2 * t1 + t), TT, -1))
    } else {
      S <- if (A) 2 * t1 + t else 2 * t1 + 2 * t   # sum of tip starts after t1
      rbind(.row3(t,  2, d,     2 * t,            t,  0),
            .row3(t1, 3, T_TIP, 3 * t1 - 2 * t,   t1, if (A) t else 0),
            .row3(d,  4, T_EXT, 4 * d - S,        t,  t),
            .row3(TT, 3, T_TIP, 3 * TT - (S - t), TT, -1))
    }
  } else {
    tau <- as.numeric(mode)
    S <- 2 * t1 - tau + 2 * t
    rbind(.row3(t1, 2, T_TIP, 2 * t1,           t1, 0),
          .row3(t,  3, d,     3 * t - 2 * t1,   t,  tau),
          .row3(d,  4, T_EXT, 4 * d - S,        t,  t),
          .row3(TT, 3, T_TIP, 3 * TT - (S - t), TT, -1))
  }
}

.gl <- function(n) {                               # Golub-Welsch, mapped to [0,1]
  k <- 1:(n - 1); b <- k / sqrt(4 * k^2 - 1)
  J <- diag(0, n); J[cbind(k, k + 1)] <- b; J[cbind(k + 1, k)] <- b
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = (e$values[o] + 1) / 2, w = e$vectors[1, o]^2)
}

# int f over the one-missing stratum, summed over the labelled attachments.
.L1_3 <- function(TT, t1, p, mb, lk, n) {
  gq <- .gl(n)
  trees <- list(); wgt <- numeric(0)
  add <- function(tr, w) { trees[[length(trees) + 1L]] <<- tr; wgt <<- c(wgt, w) }
  for (tp in list(c(0, t1), c(t1, TT))) {          # t either side of the split
    ts <- tp[1] + (tp[2] - tp[1]) * gq$x
    for (ix in seq_len(n)) {
      t <- ts[ix]; wt <- gq$w[ix] * (tp[2] - tp[1])
      dps <- if (t < t1) list(c(t, t1), c(t1, TT)) else list(c(t, TT))
      for (dp in dps) {                            # and d either side of it
        ds <- dp[1] + (dp[2] - dp[1]) * gq$x
        for (iy in seq_len(n)) {
          d <- ds[iy]; w <- wt * gq$w[iy] * (dp[2] - dp[1])
          if (t < t1) {
            add(.aug1_3(t, d, TT, t1, "A"), 2 * w)
            add(.aug1_3(t, d, TT, t1, "B"), 2 * w)
          } else {
            add(.aug1_3(t, d, TT, t1, t1), 4 * w)  # A and C, two histories each
            add(.aug1_3(t, d, TT, t1, 0),  2 * w)  # B
          }
        }
      }
    }
  }
  sum(wgt * exp(eval_logf(p, trees, model = mb, link = lk, rho = 1)$logf))
}

.cfg3 <- list(TT = 3, t1 = 2.4, mb = c(0L, 0L, 1L), link = 0L,
              p = c(0.35, 0, 0, 0.35, 0.30, 0, 0, 0))


# --------------------------------------------------------------------------- #
#  The sweep, from the recursion                                               #
# --------------------------------------------------------------------------- #
#
# A tree with two observed branchings has too many cases to write out by hand,
# so the reference is the recursion itself, applied to an event list:
#
#   P(t) = N t - S,  S the sum of tip starts over the lineages alive
#   p splits at t :  S += (t - ts[p]) + t;  ts[p] <- ts[daughter] <- t;  N + 1
#   s dies at d   :  S -= ts[s];  N - 1
#
# node.pd is P at the node's own time over the lineages alive on the segment
# that ends there, node.n is that same count, and node.focal_tip_start is the
# tip start the splitting (or dying) lineage carried just before the event.
# The two crown lineages start at tip start 0 under the reserved ids.  The 3-tip
# references above are written out by hand and are checked against this below.
.sweep_ref <- function(events) {
  ts <- list(); ts[["-2"]] <- 0; ts[["-3"]] <- 0
  N <- 2; S <- 0; rows <- vector("list", length(events))
  for (k in seq_along(events)) {
    e <- events[[k]]; t <- e$t
    r <- list(brts = t, n = N, t_ext = e$t_ext, pd = N * t - S,
              id = e$id, parent_id = e$parent_id)
    if (e$kind == "spec") {
      r$focal_tip_start <- ts[[e$lin]]
      S <- S + (t - ts[[e$lin]]) + t
      ts[[e$lin]] <- t; ts[[e$new]] <- t; N <- N + 1
      r$tip_start <- t
    } else if (e$kind == "ext") {
      r$focal_tip_start <- ts[[e$lin]]; r$tip_start <- ts[[e$lin]]
      S <- S - ts[[e$lin]]; N <- N - 1; ts[[e$lin]] <- NULL
    } else {
      # the present: a terminal marker, not an event
      r$tip_start <- t; r$focal_tip_start <- -1
    }
    rows[[k]] <- as.data.frame(r)
  }
  do.call(rbind, rows)[, c("brts", "n", "t_ext", "pd", "tip_start",
                           "focal_tip_start", "id", "parent_id")]
}

.ev <- function(t, kind, lin = NA, new = NA, id, parent_id, t_ext)
  list(t = t, kind = kind, lin = lin, new = new, id = id,
       parent_id = as.integer(parent_id), t_ext = t_ext)

.obs3 <- function(TT, t1) list(
  .ev(t1, "spec", "-2", "0", 0L, -1L, T_TIP),
  .ev(TT, "last", id = 1L, parent_id = -1L, t_ext = T_TIP))

.aug3 <- function(t, d, TT, t1, p) {
  ev <- c(.obs3(TT, t1),
          list(.ev(t, "spec", p, "2", 2L, as.integer(p), d),
               .ev(d, "ext", "2", NA, 2L, as.integer(p), T_EXT)))
  .sweep_ref(ev[order(vapply(ev, function(e) e$t, 0))])
}


# --------------------------------------------------------------------------- #
#  The balanced 4-tip tree                                                     #
# --------------------------------------------------------------------------- #
#
# The crown at 0 splits into A and B; A splits again at t1 and B at t2, so both
# observed branchings are of a lineage that was last a tip at the crown and
# parent_tip_start is c(0, 0, -1).  A is crown id -2 -- the sweep splits it at
# the first observed event -- and B is -3; their daughters are node ids 0 and 1.
# The attachments a birth has are
#
#   t < t1        A, B              K = 4
#   t1 < t < t2   A, B, C0          K = 6
#   t2 < t        A, B, C0, C1      K = 8,
#
# every candidate an observed lineage and so worth two labelled histories.
#
# This is the smallest tree on which the sampler's forward pass and a replay of
# it can disagree about WHICH lineage splits at an observed event: at t2 the
# alive observed lineages carry three different tip starts, so a match key that
# has been moved by an augmented split can resolve the wrong one.
.cfg4 <- list(TT = 3, t1 = 0.5, t2 = 2.4)
.brts4 <- function(cf) c(cf$TT, cf$TT - cf$t1, cf$TT - cf$t2)
.pts4 <- c(0, 0, -1)

.obs4 <- function(TT, t1, t2) list(
  .ev(t1, "spec", "-2", "0", 0L, -1L, T_TIP),
  .ev(t2, "spec", "-3", "1", 1L, -1L, T_TIP),
  .ev(TT, "last", id = 2L, parent_id = -1L, t_ext = T_TIP))

.aug4 <- function(t, d, TT, t1, t2, p) {
  ev <- c(.obs4(TT, t1, t2),
          list(.ev(t, "spec", p, "3", 3L, as.integer(p), d),
               .ev(d, "ext", "3", NA, 3L, as.integer(p), T_EXT)))
  .sweep_ref(ev[order(vapply(ev, function(e) e$t, 0))])
}
.aug0_4 <- function(TT, t1, t2) .sweep_ref(.obs4(TT, t1, t2))

# The attachments of the one-missing stratum, enumerated against the topology as
# supplied: the lineages alive in each stratum, each carrying two histories.
.strata4 <- function(TT, t1, t2)
  list(list(tp = c(0, t1),  par = c("-2", "-3")),
       list(tp = c(t1, t2), par = c("-2", "-3", "0")),
       list(tp = c(t2, TT), par = c("-2", "-3", "0", "1")))

# int f over that stratum.  d is integrated piece by piece between the observed
# event times, where f is smooth.
.L1_strata <- function(strata, cuts, build, pars, mb, lk, n) {
  gq <- .gl(n)
  trees <- list(); wgt <- numeric(0)
  add <- function(tr, w) { trees[[length(trees) + 1L]] <<- tr; wgt <<- c(wgt, w) }
  for (st in strata) {
    tq <- st$tp[1] + (st$tp[2] - st$tp[1]) * gq$x
    for (ix in seq_len(n)) {
      t <- tq[ix]; wt <- gq$w[ix] * (st$tp[2] - st$tp[1])
      br <- unique(sort(c(t, cuts))); br <- br[br >= t]
      for (j in seq_len(length(br) - 1L)) {
        lo <- br[j]; hi <- br[j + 1L]
        ds <- lo + (hi - lo) * gq$x
        for (iy in seq_len(n)) {
          w <- wt * gq$w[iy] * (hi - lo)
          for (pp in st$par) add(build(t, ds[iy], pp), 2 * w)
        }
      }
    }
  }
  sum(wgt * exp(eval_logf(pars, trees, model = mb, link = lk, rho = 1)$logf))
}

.L1_4 <- function(TT, t1, t2, pars, mb, lk, n)
  .L1_strata(.strata4(TT, t1, t2), c(t1, t2, TT),
             function(t, d, p) .aug4(t, d, TT, t1, t2, p), pars, mb, lk, n)

# The three models gate 1 is run under on the 4-tip tree: D alone, N + D under
# the exponential link, and one with M in the predictor as well.
.cfgs4 <- list(
  list(name = "d linear", mb = c(0L, 0L, 1L), link = 0L,
       p = c(0.30, 0, 0, 0.45, 0.25, 0, 0, 0.12)),
  list(name = "nd exp",   mb = c(1L, 0L, 1L), link = 1L,
       p = c(log(0.45), -0.04, 0, 0.30, log(0.25), 0, 0, 0.12)),
  list(name = "M + D",    mb = c(0L, 1L, 1L), link = 0L,
       p = c(0.30, 0, 0.10, 0.40, 0.25, 0, 0, 0.10)))


# --------------------------------------------------------------------------- #
#  The same, with no topology at all                                           #
# --------------------------------------------------------------------------- #
#
# What a bare branching-time vector supports: every observed lineage dates from
# the crown, the two crown lineages are not in the sum, and a node counts
# itself, so
#
#   P(t) = count*t - S over the nodes consumed, plus the node's own arrival
#          (t - its tip start) or departure at its own time
#   an observed node   : tip start 0, no splitting lineage on record
#   an augmented node  : tip start its own birth; the tip start its parent last
#                        reported, which is 0 for a crown lineage
#   an extinction node : tip start and focal are the lineage's own birth
#
# node.n is the lineage count, which this sum does not carry.
.sweep_ref_legacy <- function(events) {
  count <- 0; S <- 0; N <- 2
  born <- list(); reset <- list(); reset[["-2"]] <- 0; reset[["-3"]] <- 0
  rows <- vector("list", length(events))
  for (k in seq_along(events)) {
    e <- events[[k]]; t <- e$t
    key <- as.character(e$id); pkey <- as.character(e$parent_id)
    if (e$kind == "ext") {
      ts <- born[[key]]
      r <- list(brts = t, n = N, t_ext = e$t_ext, pd = count * t - S - (t - ts),
                tip_start = e$t_spec, focal_tip_start = e$t_spec)
      count <- count - 1; S <- S - ts; N <- N - 1
      born[[key]] <- NULL; reset[[key]] <- NULL
    } else {
      has_p <- e$parent_id != -1L
      ts <- if (has_p) t else 0
      focal <- -1
      if (has_p) {
        focal <- if (is.null(reset[[pkey]])) 0 else reset[[pkey]]
        if (!is.null(reset[[pkey]])) reset[[pkey]] <- t
      }
      r <- list(brts = t, n = N, t_ext = e$t_ext, pd = count * t - S + (t - ts),
                tip_start = ts, focal_tip_start = focal)
      count <- count + 1; S <- S + ts; N <- N + 1
      born[[key]] <- ts; reset[[key]] <- t
    }
    r$id <- e$id; r$parent_id <- e$parent_id
    rows[[k]] <- as.data.frame(r)
  }
  do.call(rbind, rows)[, c("brts", "n", "t_ext", "pd", "tip_start",
                           "focal_tip_start", "id", "parent_id")]
}

.evl <- function(t, kind, id, parent_id, t_ext, t_spec = NA)
  list(t = t, kind = kind, id = id, parent_id = as.integer(parent_id),
       t_ext = t_ext, t_spec = t_spec)

.obs3l <- function(TT, t1) list(.evl(t1, "spec", 0L, -1L, T_TIP),
                                .evl(TT, "spec", 1L, -1L, T_TIP))

.aug3l <- function(t, d, TT, t1, p) {
  ev <- c(.obs3l(TT, t1),
          list(.evl(t, "spec", 2L, as.integer(p), d),
               .evl(d, "ext", 2L, as.integer(p), T_EXT, t_spec = t)))
  .sweep_ref_legacy(ev[order(vapply(ev, function(e) e$t, 0))])
}
.aug0_3l <- function(TT, t1) .sweep_ref_legacy(.obs3l(TT, t1))

# The candidates are the same lineages as with a topology -- the two crown
# lineages before the split, and the lineage born at it after -- so K is 4 and
# then 6.  What differs is that both crown lineages report a tip start of 0.
.L1_3l <- function(TT, t1, pars, mb, lk, n)
  .L1_strata(list(list(tp = c(0, t1),  par = c("-2", "-3")),
                  list(tp = c(t1, TT), par = c("-2", "-3", "0"))),
             c(t1, TT), function(t, d, p) .aug3l(t, d, TT, t1, p),
             pars, mb, lk, n)

.cfgs3l <- list(
  list(name = "d linear", mb = c(0L, 0L, 1L), link = 0L,
       p = c(0.35, 0, 0, 0.35, 0.30, 0, 0, 0.10)),
  list(name = "nd exp",   mb = c(1L, 0L, 1L), link = 1L,
       p = c(log(0.45), -0.05, 0, 0.25, log(0.25), 0, 0, 0.10)))


# The columns the sweep is responsible for, plus the ones that identify a node.
.SWEPT <- c("brts", "n", "t_ext", "pd", "tip_start", "focal_tip_start",
            "clade", "id", "parent_id")

.expect_tree_equal <- function(got, ref, cols = .SWEPT, info = NULL) {
  for (cl in cols)
    expect_equal(as.numeric(got[[cl]]), as.numeric(ref[[cl]]),
                 tolerance = 1e-10, info = paste(cl, info))
}


test_that("the 3-tip augmented tree is the one the sampler builds", {
  skip_on_cran()
  cf <- .cfg3
  a <- augment_trees(c(cf$TT, cf$TT - cf$t1), cf$p, 3000L, 300000L, 200L, 1e6, 1L,
                     model = cf$mb, link = cf$link, rho = 1,
                     parent_tip_start = c(0, -1))
  expect_equal(a$envelope_violations, 0)
  nm <- vapply(a$trees, function(d) sum(.is_mis(d$t_ext)), 0)
  expect_gt(sum(nm == 1L), 100L)
  expect_equal(as.numeric(a$trees[[which(nm == 0L)[1L]]]$pd),
               as.numeric(.aug0_3(cf$TT, cf$t1)$pd), tolerance = 1e-12)
  seen <- character(0)
  for (j in which(nm == 1L)) {
    df  <- a$trees[[j]]
    i   <- which(.is_mis(df$t_ext))
    t   <- df$brts[i]; d <- df$t_ext[i]; pid <- df$parent_id[i]
    # every attachment is one of the five the stratum has, and none is "-1"
    if (t < cf$t1) {
      expect_true(pid %in% c(CROWN_A, CROWN_B))
      mode <- if (pid == CROWN_A) "A" else "B"
    } else {
      expect_true(pid %in% c(CROWN_A, CROWN_B, 0L))
      mode <- if (pid == CROWN_B) 0 else cf$t1
    }
    seen <- union(seen, paste0(if (t < cf$t1) "lt" else "gt", pid))
    ref <- .aug1_3(t, d, cf$TT, cf$t1, mode)
    for (col in c("brts", "n", "t_ext", "pd", "tip_start", "focal_tip_start"))
      expect_equal(as.numeric(df[[col]]), as.numeric(ref[[col]]),
                   tolerance = 1e-12, info = col)
  }
  # all five attachments are actually drawn, the two crown lineages included
  expect_setequal(seen, c("lt-2", "lt-3", "gt-2", "gt-3", "gt0"))
})


# --------------------------------------------------------------------------- #
#  Gate 1: E_q[f/q] is the marginal likelihood on a tree with an observed split #
# --------------------------------------------------------------------------- #

test_that("E_q[f/q] is the brute-force marginal likelihood on a 3-tip tree", {
  skip_on_cran()
  cf <- .cfg3
  brts <- c(cf$TT, cf$TT - cf$t1)
  l0 <- exp(eval_logf(cf$p, list(.aug0_3(cf$TT, cf$t1)), model = cf$mb,
                      link = cf$link, rho = 1)$logf)
  l1 <- .L1_3(cf$TT, cf$t1, cf$p, cf$mb, cf$link, 40L)
  # the quadrature is converged: halving the order moves it by < 1e-4 relative
  expect_lt(abs(l1 - .L1_3(cf$TT, cf$t1, cf$p, cf$mb, cf$link, 20L)) / l1, 1e-4)

  tot1 <- ss1 <- totle <- ssle <- w2 <- 0; N <- 0L; nz <- 0L
  for (b in 1:3) {
    a <- augment_trees(brts, cf$p, 20000L, 1000000L, 200L, 1e6, 1L, model = cf$mb,
                       link = cf$link, rho = 1, parent_tip_start = c(0, -1))
    expect_equal(a$rejected_overruns, 0)
    expect_equal(a$rejected_lambda, 0)
    expect_equal(a$rejected_nonfinite, 0)
    expect_equal(a$envelope_violations, 0)
    nz <- nz + a$rejected_zero_weights
    w  <- exp(eval_logf(cf$p, a$trees, model = cf$mb, link = cf$link,
                        rho = 1)$logf - a$logg)
    nm <- vapply(a$trees, function(d) sum(.is_mis(d$t_ext)), 0)
    wle <- ifelse(nm <= 1, w, 0)
    N <- N + length(w)
    tot1 <- tot1 + sum(w[nm == 1]); ss1 <- ss1 + sum(w[nm == 1]^2)
    totle <- totle + sum(wle); ssle <- ssle + sum(wle^2)
    w2 <- w2 + sum(w[nm >= 2])
  }
  M <- N + nz
  e1  <- tot1 / M;  s1  <- sqrt(ss1 / M - e1^2) / sqrt(M)
  ele <- totle / M; sle <- sqrt(ssle / M - ele^2) / sqrt(M)
  # On the build before this one the sampler could not reach the crown
  # lineages, and the same integral came out 24 % above the m = 1 stratum
  # (0.00977 against 0.00789): about 30 standard errors at this sample size.
  expect_lt(abs(e1 - l1) / s1, 5,
            label = "|IS - L1| / se over the one-missing draws")
  expect_lt(abs(e1 - l1) / l1, 0.05,
            label = "|IS - L1| / L1 over the one-missing draws")
  expect_lt(abs(ele - l0 - l1) / sle, 5,
            label = "|IS - (L0 + L1)| / se over the draws with at most one")
  # the truncation, measured rather than assumed: the draws left out
  expect_gt(w2 / (totle + w2), 0)       # there are some
  expect_lt(w2 / (totle + w2), 0.30)    # and they are a stated fraction
})


# --------------------------------------------------------------------------- #
#  Gate 4: the candidates are the lineages alive                               #
# --------------------------------------------------------------------------- #

test_that("every lineage alive is a candidate, and q counts the same attachments", {
  skip_on_cran()
  set.seed(11)
  phy  <- ape::rphylo(8L, 0.8, 0.25)
  brts <- emphasis:::.extract_brts(phy)
  pts  <- emphasis:::.pts(brts)
  p    <- c(0.9, -0.02, 0, 0.10, 0.35, 0, 0, 0.03)
  # with the topology and without it, and on the 3-tip tree with its split
  cases <- list(list(b = as.numeric(brts), pts = pts, mb = c(1L, 0L, 1L)),
                list(b = as.numeric(brts), pts = numeric(0), mb = c(1L, 0L, 1L)),
                list(b = c(3, 0.6), pts = c(0, -1), mb = c(0L, 0L, 1L)))
  checked <- 0L
  for (cs in cases) {
    a <- augment_trees(cs$b, p, 60L, 20000L, 400L, 1e6, 1L, model = cs$mb,
                       link = 0L, rho = 1, parent_tip_start = cs$pts)
    skip_if(length(a$trees) == 0L, "augmentation drew no tree")
    for (df in a$trees) {
      rep <- eval_attachments(df)
      if (nrow(rep) == 0L) next
      # the alive count node.n carries, not node.n - 2
      expect_equal(rep$candidates, rep$n)
      # the parent drawn is one of them, and is never "no parent on record"
      expect_true(all(rep$parent_alive))
      expect_false(any(rep$parent_id == -1L))
      # and the labelled attachments are the 2 * tips + Ne of sampling_prob,
      # recomputed here from the node list
      tips <- df$n[1L]; Ne <- 0; K <- numeric(0)
      for (i in seq_len(nrow(df))) {
        tips <- tips + (df$t_ext[i] == T_TIP); Ne <- Ne - (df$t_ext[i] == T_EXT)
        if (.is_mis(df$t_ext[i])) { K <- c(K, 2 * tips + Ne); Ne <- Ne + 1 }
      }
      expect_equal(rep$attachments, K)
      checked <- checked + nrow(rep)
    }
  }
  expect_gt(checked, 200L)
})


# --------------------------------------------------------------------------- #
#  Gate 3: every augmented lineage reaches the tas                             #
# --------------------------------------------------------------------------- #

test_that("no augmented lineage is dropped from the reconstructed tas", {
  skip_on_cran()
  set.seed(4)
  for (phy in list(ape::rphylo(8L, 0.8, 0.25), ape::rphylo(20L, 0.7, 0.2))) {
    brts <- emphasis:::.extract_brts(phy)
    L_extant <- DDD::phylo2L(phy)
    a <- emphasis:::.augment_tree_internal(
      phy, pars = c(0.9, -0.02, 0.10, 0.35, 0, 0.03), model_bin = c(1L, 0L, 1L),
      sample_size = 120L, max_missing = 400, max_lambda = 1e6, maxN = 40000L,
      num_threads = 1L, link = 0L, rho = 1)
    skip_if(length(a$trees) == 0L, "augmentation drew no tree")
    n_aug <- 0L
    for (df in a$trees) {
      mis <- df$t_ext != 0 & df$t_ext < 1e11
      expect_false(any(df$parent_id[mis] == -1L))
      L <- emphasis:::.aug_to_Ltable(df, brts[1L], brts, L_extant)
      expect_false(is.null(L))
      # one row per augmented lineage: nothing is filtered out
      expect_equal(nrow(L), nrow(L_extant) + sum(mis))
      tas <- DDD::L2phylo(L, dropextinct = FALSE)
      expect_gte(min(tas$edge.length), 0)
      expect_equal(ape::Ntip(prune_to_extant(tas)), ape::Ntip(phy))
      n_aug <- n_aug + sum(mis)
    }
    expect_gt(n_aug, 40L)
  }
})


# --------------------------------------------------------------------------- #
#  The recursion reference agrees with the 3-tip trees written out by hand     #
# --------------------------------------------------------------------------- #

test_that("the sweep reference reproduces the hand-written 3-tip trees", {
  TT <- .cfg3$TT; t1 <- .cfg3$t1
  cols <- c("brts", "n", "t_ext", "pd", "tip_start", "focal_tip_start")
  .expect_tree_equal(.sweep_ref(.obs3(TT, t1)), .aug0_3(TT, t1), cols = cols)
  for (t in c(0.4, 1.7, 2.39)) for (d in c(t + 0.05, 2.5, 2.9)) {
    if (d <= t) next
    .expect_tree_equal(.aug3(t, d, TT, t1, "-2"),
                       .aug1_3(t, d, TT, t1, if (t < t1) "A" else t1),
                       cols = cols, info = paste("A", t, d))
    .expect_tree_equal(.aug3(t, d, TT, t1, "-3"),
                       .aug1_3(t, d, TT, t1, if (t < t1) "B" else 0),
                       cols = cols, info = paste("B", t, d))
  }
})


# --------------------------------------------------------------------------- #
#  The 4-tip tree: the draws are the trees the recursion builds                #
# --------------------------------------------------------------------------- #

test_that("the 4-tip augmented tree is the one the sampler builds", {
  skip_on_cran()
  cf <- .cfg4
  cols <- c("brts", "n", "t_ext", "pd", "tip_start", "focal_tip_start")
  # A constant-rate model, so that every attachment survives to be counted.
  # The tip starts and pd the sweep writes do not depend on the model at all,
  # but f does: under a strongly D-dependent set some attachments make a rate
  # negative and carry f = 0, and a tree of weight 0 is never returned.  That
  # is not a support gap -- an attachment f gives no density to needs none --
  # but it would hide which attachments the sampler can reach.
  a <- augment_trees(.brts4(cf), c(0.35, 0, 0, 0, 0.30, 0, 0, 0), 6000L,
                     600000L, 200L, 1e6, 1L, model = c(0L, 0L, 0L), link = 0L,
                     rho = 1, parent_tip_start = .pts4)
  expect_equal(a$rejected_zero_weights, 0)
  expect_equal(a$envelope_violations, 0)
  nm <- vapply(a$trees, function(d) sum(.is_mis(d$t_ext)), 0)
  expect_gt(sum(nm == 1L), 100L)
  .expect_tree_equal(a$trees[[which(nm == 0L)[1L]]],
                     .aug0_4(cf$TT, cf$t1, cf$t2), cols = cols)
  seen <- character(0)
  for (j in which(nm == 1L)) {
    df <- a$trees[[j]]
    i <- which(.is_mis(df$t_ext))
    t <- df$brts[i]; d <- df$t_ext[i]; pid <- df$parent_id[i]
    # the parent is one of the lineages alive at t, and never "no parent"
    stratum <- match(TRUE, t < c(cf$t1, cf$t2, Inf))
    alive <- list(c(CROWN_A, CROWN_B), c(CROWN_A, CROWN_B, 0L),
                  c(CROWN_A, CROWN_B, 0L, 1L))[[stratum]]
    expect_true(pid %in% alive)
    seen <- union(seen, paste0(stratum, "/", pid))
    .expect_tree_equal(df, .aug4(t, d, cf$TT, cf$t1, cf$t2, as.character(pid)),
                       cols = cols,
                       info = paste("t", t, "d", d, "pid", pid))
  }
  # all nine attachments are actually drawn
  expect_setequal(seen, c("1/-2", "1/-3",
                          "2/-2", "2/-3", "2/0",
                          "3/-2", "3/-3", "3/0", "3/1"))
})


# --------------------------------------------------------------------------- #
#  Gate 1 on the 4-tip tree: two observed splits, K = 4, 6 and 8               #
# --------------------------------------------------------------------------- #
#
# A 3-tip tree cannot separate a sampler that resolves the splitting lineage
# from the caller's topology from one that resolves it from a key the sweep has
# moved: before its one split both crown lineages report a tip start of 0.  On
# this tree the second split has three different tip starts alive, and with the
# key moved the replay resolved a different lineage on 61 % of the draws, for
# z = +18.0 (d linear), +8.3 (nd exp) and +17.6 (M + D) against this integral.

test_that("E_q[f/q] is the brute-force marginal likelihood on a 4-tip tree", {
  skip_on_cran()
  cf <- .cfg4
  brts <- .brts4(cf)
  for (cfg in .cfgs4) {
    l1 <- .L1_4(cf$TT, cf$t1, cf$t2, cfg$p, cfg$mb, cfg$link, 20L)
    # the quadrature is converged: halving the order moves it by < 1e-3 relative
    expect_lt(abs(l1 - .L1_4(cf$TT, cf$t1, cf$t2, cfg$p, cfg$mb, cfg$link, 10L)) / l1,
              1e-3, label = paste("quadrature drift,", cfg$name))

    tot1 <- ss1 <- totle <- w2 <- 0; N <- 0L; nz <- 0L
    for (b in 1:3) {
      a <- augment_trees(brts, cfg$p, 20000L, 1000000L, 200L, 1e6, 1L,
                         model = cfg$mb, link = cfg$link, rho = 1,
                         parent_tip_start = .pts4)
      expect_equal(a$rejected_overruns, 0)
      expect_equal(a$rejected_lambda, 0)
      expect_equal(a$rejected_nonfinite, 0)
      expect_equal(a$envelope_violations, 0)
      nz <- nz + a$rejected_zero_weights
      w <- exp(eval_logf(cfg$p, a$trees, model = cfg$mb, link = cfg$link,
                         rho = 1)$logf - a$logg)
      nm <- vapply(a$trees, function(d) sum(.is_mis(d$t_ext)), 0)
      N <- N + length(w)
      tot1 <- tot1 + sum(w[nm == 1]); ss1 <- ss1 + sum(w[nm == 1]^2)
      totle <- totle + sum(ifelse(nm <= 1, w, 0)); w2 <- w2 + sum(w[nm >= 2])
    }
    M <- N + nz
    e1 <- tot1 / M; s1 <- sqrt(ss1 / M - e1^2) / sqrt(M)
    expect_lt(abs(e1 - l1) / s1, 5, label = paste("|IS - L1| / se,", cfg$name))
    expect_lt(abs(e1 - l1) / l1, 0.05, label = paste("|IS - L1| / L1,", cfg$name))
    # the truncation, measured rather than assumed: the draws left out
    expect_gt(w2 / (totle + w2), 0)
    expect_lt(w2 / (totle + w2), 0.45)
  }
})


# --------------------------------------------------------------------------- #
#  Gate 2: the tree that comes back is the state the sampler carried           #
# --------------------------------------------------------------------------- #
#
# The sweep resolves an observed branching by the tip start the caller gave it,
# and writes that lineage's current tip start back over the same field.  The
# closing pass replays the whole tree, so it must be handed the caller's key
# again; driven by what the forward pass wrote it resolves a different lineage,
# on 61 % of the draws from the 4-tip tree here, 62 % from an 8-tip and 89 %
# from a 20-tip.  A 2- or 3-tip tree shows 0 % and cannot see it.

test_that("the augmented tree is one application of the sweep the caller drove", {
  skip_on_cran()
  cf <- .cfg4
  set.seed(11); phy8 <- ape::rphylo(8L, 0.8, 0.25)
  set.seed(4);  phy20 <- ape::rphylo(20L, 0.7, 0.2)
  b8 <- emphasis:::.extract_brts(phy8); b20 <- emphasis:::.extract_brts(phy20)
  cases <- list(list(nm = "4-tip", b = .brts4(cf), pts = .pts4),
                list(nm = "8-tip", b = as.numeric(b8), pts = emphasis:::.pts(b8)),
                list(nm = "20-tip", b = as.numeric(b20), pts = emphasis:::.pts(b20)))
  p <- c(0.9, -0.02, 0, 0.10, 0.35, 0, 0, 0.03)
  for (cs in cases) {
    a <- augment_trees(cs$b, p, 200L, 200000L, 400L, 1e6, 1L,
                       model = c(1L, 0L, 1L), link = 0L, rho = 1,
                       parent_tip_start = cs$pts)
    skip_if(length(a$trees) == 0L, "augmentation drew no tree")
    n_aug <- 0L
    for (df in a$trees) {
      n_aug <- n_aug + sum(.is_mis(df$t_ext))
      .expect_tree_equal(eval_pendant_sweep(df, cs$pts), df, info = cs$nm)
    }
    expect_gt(n_aug, 100L)
  }
})


# --------------------------------------------------------------------------- #
#  The attachment report does not depend on the pendant-age convention         #
# --------------------------------------------------------------------------- #
#
# eval_attachments replays the tree through whichever sweep the `clade` column
# records.  What it reports -- the alive set, the labelled attachments and
# whether the recorded parent is among them -- is carried by the lineage ids
# and the event order alone, which the two conventions share, so the report is
# the same under either.  Held here because the docstrings say so.

test_that("the attachment report is the same under both conventions", {
  skip_on_cran()
  cf <- .cfg4; cfg <- .cfgs4[[1]]
  a <- augment_trees(.brts4(cf), cfg$p, 100L, 100000L, 200L, 1e6, 1L,
                     model = cfg$mb, link = cfg$link, rho = 1,
                     parent_tip_start = .pts4)
  skip_if(length(a$trees) == 0L, "augmentation drew no tree")
  checked <- 0L
  for (df in a$trees) {
    expect_true("clade" %in% names(df))
    expect_equal(df$clade[nrow(df)], 1L)      # the topology sweep is on record
    legacy <- df; legacy$clade <- NULL        # read as a bare branching-time tree
    rep_t <- eval_attachments(df)
    expect_equal(rep_t, eval_attachments(legacy))
    checked <- checked + nrow(rep_t)
  }
  expect_gt(checked, 50L)
})


# --------------------------------------------------------------------------- #
#  Gate 7: the same integral with no topology                                  #
# --------------------------------------------------------------------------- #
#
# Drawing the parent from the alive set changed this path too: before it, a
# birth earlier than the first observed node had no parent on record, dated
# from the crown and fell back to D = 0.  It now names a crown lineage, dates
# from its own birth like every other augmented lineage, and reads the tip
# start the legacy convention gives an observed lineage.  No gate covered that,
# so here it is: E_q[f/q] against the integral of f over the same attachments.

test_that("E_q[f/q] is the brute-force marginal likelihood with no topology", {
  skip_on_cran()
  TT <- 3; t1 <- 1.1
  for (cfg in .cfgs3l) {
    l1 <- .L1_3l(TT, t1, cfg$p, cfg$mb, cfg$link, 20L)
    expect_lt(abs(l1 - .L1_3l(TT, t1, cfg$p, cfg$mb, cfg$link, 10L)) / l1, 1e-3,
              label = paste("quadrature drift,", cfg$name))
    l0 <- exp(eval_logf(cfg$p, list(.aug0_3l(TT, t1)), model = cfg$mb,
                        link = cfg$link, rho = 1)$logf)

    tot1 <- ss1 <- totle <- ssle <- w2 <- 0; N <- 0L; nz <- 0L
    for (b in 1:3) {
      a <- augment_trees(c(TT, TT - t1), cfg$p, 20000L, 1000000L, 200L, 1e6, 1L,
                         model = cfg$mb, link = cfg$link, rho = 1,
                         parent_tip_start = numeric(0))
      expect_equal(a$rejected_overruns, 0)
      expect_equal(a$rejected_lambda, 0)
      expect_equal(a$rejected_nonfinite, 0)
      expect_equal(a$envelope_violations, 0)
      nz <- nz + a$rejected_zero_weights
      w <- exp(eval_logf(cfg$p, a$trees, model = cfg$mb, link = cfg$link,
                         rho = 1)$logf - a$logg)
      nm <- vapply(a$trees, function(d) sum(.is_mis(d$t_ext)), 0)
      wle <- ifelse(nm <= 1, w, 0)
      N <- N + length(w)
      tot1 <- tot1 + sum(w[nm == 1]); ss1 <- ss1 + sum(w[nm == 1]^2)
      totle <- totle + sum(wle); ssle <- ssle + sum(wle^2); w2 <- w2 + sum(w[nm >= 2])
    }
    M <- N + nz
    e1 <- tot1 / M; s1 <- sqrt(ss1 / M - e1^2) / sqrt(M)
    ele <- totle / M; sle <- sqrt(ssle / M - ele^2) / sqrt(M)
    expect_lt(abs(e1 - l1) / s1, 5, label = paste("|IS - L1| / se,", cfg$name))
    expect_lt(abs(e1 - l1) / l1, 0.05, label = paste("|IS - L1| / L1,", cfg$name))
    expect_lt(abs(ele - l0 - l1) / sle, 5,
              label = paste("|IS - (L0 + L1)| / se,", cfg$name))
    expect_gt(w2 / (totle + w2), 0)
    expect_lt(w2 / (totle + w2), 0.45)
  }
})


# --------------------------------------------------------------------------- #
#  The drawn tree with no topology is the one the legacy recursion builds      #
# --------------------------------------------------------------------------- #

test_that("the no-topology augmented tree is the one the sampler builds", {
  skip_on_cran()
  TT <- 3; t1 <- 1.1
  cols <- c("brts", "n", "t_ext", "pd", "tip_start", "focal_tip_start")
  a <- augment_trees(c(TT, TT - t1), c(0.35, 0, 0, 0, 0.30, 0, 0, 0), 3000L,
                     300000L, 200L, 1e6, 1L, model = c(0L, 0L, 0L), link = 0L,
                     rho = 1, parent_tip_start = numeric(0))
  expect_equal(a$rejected_zero_weights, 0)
  nm <- vapply(a$trees, function(d) sum(.is_mis(d$t_ext)), 0)
  expect_gt(sum(nm == 1L), 100L)
  .expect_tree_equal(a$trees[[which(nm == 0L)[1L]]], .aug0_3l(TT, t1), cols = cols)
  seen <- character(0)
  for (j in which(nm == 1L)) {
    df <- a$trees[[j]]
    i <- which(.is_mis(df$t_ext))
    t <- df$brts[i]; d <- df$t_ext[i]; pid <- df$parent_id[i]
    expect_true(pid %in% (if (t < t1) c(CROWN_A, CROWN_B) else c(CROWN_A, CROWN_B, 0L)))
    seen <- union(seen, paste0(if (t < t1) "lt" else "gt", pid))
    .expect_tree_equal(df, .aug3l(t, d, TT, t1, as.character(pid)), cols = cols,
                       info = paste("t", t, "d", d, "pid", pid))
  }
  expect_setequal(seen, c("lt-2", "lt-3", "gt-2", "gt-3", "gt0"))
})


# --------------------------------------------------------------------------- #
#  cr and dd read none of the fields this convention touches                   #
# --------------------------------------------------------------------------- #
#
# Everything the parent draw moves -- the tip start of the splitting lineage,
# the pendant PD it feeds, and the clade flag that says which convention wrote
# them -- enters f and q only through M and D.  An N-only model has neither in
# its predictor, so its logf and logg cannot move with any of it.  That is the
# claim under which this change is safe for cr and dd, so it is held rather
# than asserted.

test_that("logf and logg for cr and dd do not read the pendant-age fields", {
  skip_on_cran()
  set.seed(5)
  phy <- ape::rphylo(8L, 0.8, 0.25)
  b <- emphasis:::.extract_brts(phy); pts <- emphasis:::.pts(b)
  a <- augment_trees(as.numeric(b), c(0.9, -0.02, 0, 0.10, 0.35, 0, 0, 0.03),
                     50L, 200000L, 400L, 1e6, 1L, model = c(1L, 0L, 1L),
                     link = 0L, rho = 1, parent_tip_start = pts)
  skip_if(length(a$trees) == 0L, "augmentation drew no tree")
  scrambled <- lapply(a$trees, function(df) {
    df$focal_tip_start <- ifelse(df$focal_tip_start < 0, -1, df$brts * 0.37)
    df$pd <- df$pd * 1.9 + 0.5
    df$clade <- NULL
    df
  })
  for (mb in list(c(0L, 0L, 0L), c(1L, 0L, 0L))) for (lk in c(0L, 1L)) {
    p <- if (lk == 0L) c(0.6, -0.01, 0, 0, 0.2, 0.002, 0, 0)
         else c(log(0.6), -0.01, 0, 0, log(0.2), 0.002, 0, 0)
    x <- eval_logf(p, a$trees, model = mb, link = lk, rho = 1)
    y <- eval_logf(p, scrambled, model = mb, link = lk, rho = 1)
    expect_equal(x$logf, y$logf, tolerance = 1e-12,
                 info = paste(paste(mb, collapse = ""), lk))
    expect_equal(x$logg, y$logg, tolerance = 1e-12,
                 info = paste(paste(mb, collapse = ""), lk))
  }
  # and a D model does read them, so the check above is not vacuous
  p <- c(0.6, 0, 0, 0.2, 0.2, 0, 0, 0.05)
  expect_false(isTRUE(all.equal(
    eval_logf(p, a$trees, model = c(0L, 0L, 1L), link = 0L, rho = 1)$logf,
    eval_logf(p, scrambled, model = c(0L, 0L, 1L), link = 0L, rho = 1)$logf)))
})

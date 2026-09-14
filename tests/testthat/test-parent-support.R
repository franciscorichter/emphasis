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
#         among them.
# Gate 3  every augmented lineage has a parent on record and reaches the tas.
# Gate 4  the candidates are the lineages alive, and the sampler and the density
#         count the same attachments.

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

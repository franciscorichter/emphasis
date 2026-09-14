# The pendant-age covariates M and D, once the estimator is given the observed
# topology.
#
#   P(t) = N*t - sum(tip_start)      M(t) = P/N      D(s,t) = (t - ts_s) - M(t)
#
# tip_start is the forward time at which a lineage last became a pendant tip,
# so it is reset every time that lineage speciates.  These are the definitions
# the simulator uses (inst/include/general_tree.hpp), and the estimator has to
# produce the same numbers on a tree where they can be read off the topology.

# --------------------------------------------------------------------------- #
#  Reference quantities, built from the edge list alone                        #
# --------------------------------------------------------------------------- #

# Lineages alive just before t: an edge that started strictly before t and has
# not ended before it.  That is the set node.n counts on the segment ending at
# an event, so P and N are read on the same side of the event as the node.
.alive_starts <- function(phy, t) {
  h  <- ape::node.depth.edgelength(phy)
  st <- h[phy$edge[, 1L]]
  en <- h[phy$edge[, 2L]]
  st[st < t & en >= t]
}
.true_P <- function(phy, t) sum(t - .alive_starts(phy, t))
.true_N <- function(phy, t) length(.alive_starts(phy, t))

# A tree whose reconstructed form IS the complete tree: every pendant age is
# known exactly, so M and D have no missing lineages in them.
.mu0_tree <- function(seed, n_tips = 12L) {
  set.seed(seed)
  ape::rphylo(n_tips, 0.6, 0)
}

# The observed tree as the estimator builds it, with no augmented lineage:
# lambda > 0 makes the log-likelihood finite, mu = 0 makes the thinning rate
# n * lambda * (1 - exp(-mu (T - t))) vanish, so nothing is inserted.
.observed_nodes <- function(phy, with_topology = TRUE) {
  brts <- emphasis:::.extract_brts(phy)
  pts  <- if (with_topology) emphasis:::.pts(brts) else numeric(0)
  a <- augment_trees(as.numeric(brts), c(0.3, 0, 0, 0, 0, 0, 0, 0),
                     1L, 500L, 100L, 1e6, 1L,
                     model = c(0L, 0L, 0L), link = 0L, rho = 1,
                     parent_tip_start = pts)
  df <- a$trees[[1L]]
  # No lineage was augmented, so the node list is the observed one.
  expect_equal(nrow(df), length(brts))
  df
}


# --------------------------------------------------------------------------- #
#  The correctness gate: M and D against the topology                          #
# --------------------------------------------------------------------------- #

test_that("on a fully observed mu = 0 tree, M and D match the topology", {
  skip_on_cran()
  for (seed in c(4L, 11L, 23L)) {
    phy  <- .mu0_tree(seed)
    brts <- emphasis:::.extract_brts(phy)
    pts  <- emphasis:::.pts(brts)
    df   <- .observed_nodes(phy)

    # Every node but the last is an observed branching event; the last marks
    # the present.
    ev <- seq_len(nrow(df) - 1L)
    t  <- df$brts[ev]

    N_true <- vapply(t, function(x) .true_N(phy, x), 0)
    M_true <- vapply(t, function(x) .true_P(phy, x), 0) / N_true

    expect_equal(df$n[ev], N_true, tolerance = 0)
    expect_equal(df$pd[ev] / df$n[ev], M_true, tolerance = 1e-10)

    # D of the lineage that splits: its own pendant age, centred on M.
    E_pkg  <- df$brts[ev] - df$focal_tip_start[ev]
    D_pkg  <- E_pkg - df$pd[ev] / df$n[ev]
    D_true <- (t - pts[ev]) - M_true
    expect_equal(D_pkg, D_true, tolerance = 1e-10)

    # Before the fix M was the crown age t itself (to within (j+1)/(j+2)) and
    # D was identically zero; neither survives on a tree with any depth.
    expect_gt(max(abs(D_pkg)), 1e-6)
    expect_lt(min(df$pd[ev] / df$n[ev] / t), 0.9)
  }
})


test_that("sum_s D(s, t) = 0 over the alive lineages at every event", {
  skip_on_cran()
  phy <- .mu0_tree(7L)
  df  <- .observed_nodes(phy)
  ev  <- seq_len(nrow(df) - 1L)
  # D is centred by construction only if M really is the mean pendant age over
  # the lineages alive at t, so this tests M against the whole alive set and
  # not just against the splitting lineage.
  sums <- vapply(ev, function(k) {
    t <- df$brts[k]
    sum((t - .alive_starts(phy, t)) - df$pd[k] / df$n[k])
  }, 0)
  expect_equal(sums, rep(0, length(ev)), tolerance = 1e-10)
})


test_that("the tip starts handed to C++ agree with an edge-based P(t)", {
  skip_on_cran()
  # The R-side reference is one per branching event; replaying the sweep
  # S += 2t - ts_p over those events has to reproduce the pendant PD read off
  # the edge list at arbitrary probe times, not only at the events.
  for (seed in 1:5) {
    phy  <- .mu0_tree(seed, 15L)
    brts <- emphasis:::.extract_brts(phy)
    fwd  <- (brts[1L] - brts)[-1L]                 # forward event times
    pts  <- emphasis:::.observed_parent_tip_start(phy)
    expect_length(pts, length(brts) - 1L)

    P_of <- function(t) {                          # sweep to just before t
      S <- 0; N <- 2
      for (k in seq_along(fwd)) {
        if (fwd[k] >= t) break
        S <- S + 2 * fwd[k] - pts[k]
        N <- N + 1
      }
      N * t - S
    }
    probes <- seq(1e-6, brts[1L] - 1e-6, length.out = 25)
    expect_equal(vapply(probes, P_of, 0),
                 vapply(probes, function(t) .true_P(phy, t), 0),
                 tolerance = 1e-10)
  }
})


test_that("an augmented tree keeps a pendant PD a real alive set could produce", {
  skip_on_cran()
  # P is a sum over N lineages of ages in [0, t], so 0 <= P <= N*t, with N the
  # count the node itself carries.  The bound fails as soon as the sweep
  # subtracts a tip_start no alive lineage has, which is what happens if the
  # parent the augmentation names is taken on trust.
  phy  <- .mu0_tree(11L, 14L)
  brts <- emphasis:::.extract_brts(phy)
  a <- augment_trees(as.numeric(brts), c(0.9, 0, 0, 0.10, 0.35, 0, 0, 0),
                     40L, 20000L, 300L, 1e6, 1L,
                     model = c(0L, 0L, 1L), link = 0L, rho = 1,
                     parent_tip_start = emphasis:::.pts(brts))
  skip_if(length(a$trees) == 0L, "augmentation drew no tree")
  expect_true(any(vapply(a$trees, nrow, 0L) > length(brts)))   # lineages added

  for (d in a$trees) {
    expect_true(all(is.finite(d$pd)))
    expect_true(all(d$pd >= -1e-10))
    expect_true(all(d$pd <= d$n * d$brts + 1e-10))
    # n is the alive count on the segment ending at the node: two crown
    # lineages, plus one per earlier birth, less one per earlier death.
    step <- ifelse(d$t_ext == 0, -1, 1)
    expect_equal(d$n, 2 + cumsum(c(0, utils::head(step, -1L))))
  }
})


# --------------------------------------------------------------------------- #
#  cr and dd are untouched                                                     #
# --------------------------------------------------------------------------- #

test_that("cr and dd log-likelihoods do not read the pendant-age columns", {
  skip_on_cran()
  phy  <- .mu0_tree(3L, 10L)
  brts <- emphasis:::.extract_brts(phy)
  aug <- augment_trees(as.numeric(brts), c(0.9, 0, 0, 0, 0.35, 0, 0, 0),
                       6L, 4000L, 500L, 1e6, 1L,
                       model = c(0L, 0L, 0L), link = 0L, rho = 1,
                       parent_tip_start = emphasis:::.pts(brts))
  skip_if(length(aug$trees) == 0L, "augmentation drew no tree")

  # Same trees with the pendant-age columns replaced by nonsense.  Neither M
  # nor D enters a cr or dd rate, so nothing downstream may move; under the
  # D model it must.
  wrecked <- lapply(aug$trees, function(d) {
    d$pd <- d$pd + 17
    d$tip_start <- d$tip_start * 0
    d$focal_tip_start <- rep(-1, nrow(d))
    d
  })

  for (link in c(0L, 1L)) {
    p_cr <- if (link == 0L) c(0.9, 0, 0, 0, 0.35, 0, 0, 0)
            else c(-0.1, 0, 0, 0, -1.0, 0, 0, 0)
    p_dd <- if (link == 0L) c(0.9, -0.02, 0, 0, 0.35, 0, 0, 0)
            else c(-0.1, -0.01, 0, 0, -1.0, 0, 0, 0)
    for (spec in list(list(c(0L, 0L, 0L), p_cr), list(c(1L, 0L, 0L), p_dd))) {
      a <- eval_logf(spec[[2L]], aug$trees, model = spec[[1L]], link = link, rho = 1)
      b <- eval_logf(spec[[2L]], wrecked,   model = spec[[1L]], link = link, rho = 1)
      expect_identical(a$logf, b$logf)
      expect_identical(a$logg, b$logg)
    }
  }

  # The same wrecking does move the D model, so the trees really do differ.
  p_d  <- c(0.9, 0, 0, 0.15, 0.35, 0, 0, 0)
  d_ok <- eval_logf(p_d, aug$trees, model = c(0L, 0L, 1L), link = 0L, rho = 1)$logf
  d_no <- eval_logf(p_d, wrecked,   model = c(0L, 0L, 1L), link = 0L, rho = 1)$logf
  expect_false(isTRUE(all.equal(d_ok, d_no)))
})


test_that("cr and dd MCEM fits return the estimates they did before the fix", {
  skip_on_cran()
  # Pinned from the build that preceded the topology change; the BDI sampler
  # is pure R, so a seeded fit is reproducible.  cr and dd have model_bin
  # c(0,0,0) and c(1,0,0): neither M nor D enters a rate, and both fits have
  # to come back unmoved.
  set.seed(42)
  phy <- ape::rphylo(12, 0.5, 0)

  set.seed(7)
  f_cr <- estimate_rates(phy, method = "mcem", model = "cr",
                         init_pars = c(0.5, 0.1),
                         control = list(lower_bound = c(0, 0),
                                        upper_bound = c(2, 1),
                                        num_trees = 30L, max_iter = 5L))
  expect_equal(unname(f_cr$pars), c(0.5847607702017, 0.0864009000326),
               tolerance = 1e-10)
  expect_equal(f_cr$loglik, -16.0992102587, tolerance = 1e-10)

  set.seed(9)
  f_dd <- estimate_rates(phy, method = "mcem", model = "dd",
                         init_pars = c(0.8, -0.02, 0.2, 0),
                         control = list(lower_bound = c(0.1, -0.5, 0, -0.01),
                                        upper_bound = c(3, 0.01, 1, 0.01),
                                        num_trees = 30L, max_iter = 5L))
  expect_equal(unname(f_dd$pars),
               c(0.595602943635, 0.01, 0.136509120974, 0.01),
               tolerance = 1e-10)
  expect_equal(f_dd$loglik, -15.6308203252, tolerance = 1e-10)
})


# --------------------------------------------------------------------------- #
#  No topology: branching times alone still work, and behave as they did       #
# --------------------------------------------------------------------------- #

test_that("a bare branching-time vector keeps the crown-age convention", {
  skip_on_cran()
  phy <- .mu0_tree(4L)
  df  <- .observed_nodes(phy, with_topology = FALSE)

  # No parent is on record, so every observed lineage is recorded as dating
  # from the crown and no event has a splitting lineage.
  expect_equal(df$tip_start, rep(0, nrow(df)))
  expect_equal(df$focal_tip_start, rep(-1, nrow(df)))
  # One lineage per node, all of age t: P(t_j) = (j + 1) * t_j.
  expect_equal(df$pd, seq_along(df$brts) * df$brts, tolerance = 1e-12)
  # so M = (j+1)/(j+2) * t, which is the crown age up to that factor.
  M <- df$pd / df$n
  expect_equal(M / df$brts, seq_along(df$brts) / (seq_along(df$brts) + 1),
               tolerance = 1e-12)
})


test_that("without topology D is zero at every observed event, with it it is not", {
  skip_on_cran()
  phy  <- .mu0_tree(4L)
  brts <- emphasis:::.extract_brts(phy)
  mk <- function(pts) {
    augment_trees(as.numeric(brts), c(0.3, 0, 0, 0, 0, 0, 0, 0),
                  1L, 500L, 100L, 1e6, 1L, model = c(0L, 0L, 0L),
                  link = 0L, rho = 1, parent_tip_start = pts)$trees[[1L]]
  }
  flat <- mk(numeric(0))
  real <- mk(emphasis:::.pts(brts))

  # D of the lineage that splits, at every observed event.
  D_at_events <- function(tr) {
    ev <- seq_len(nrow(tr) - 1L)
    M  <- tr$pd[ev] / tr$n[ev]
    E  <- ifelse(tr$focal_tip_start[ev] >= 0,
                 tr$brts[ev] - tr$focal_tip_start[ev], M)
    E - M
  }
  expect_equal(D_at_events(flat), rep(0, nrow(flat) - 1L), tolerance = 1e-12)
  expect_gt(max(abs(D_at_events(real))), 1e-6)

  # beta_D is not invisible on the flat tree even so.  The compensator
  # integrates each alive lineage's own D over the segment, and those are not
  # zero where the event's D is: only the event terms are constant in beta_D
  # here.  While the compensator was dt * n * (rate at the segment's end node),
  # the linear-link log f of a draw with no augmented lineage was exactly
  # constant in beta_D (H8) and the covariate was invisible on both trees.
  lf <- function(tr, beta_D)
    eval_logf(c(0.3, 0, 0, beta_D, 0, 0, 0, 0), list(tr),
              model = c(0L, 0L, 1L), link = 0L, rho = 1)$logf
  expect_false(isTRUE(all.equal(lf(flat, 0.0), lf(flat, 0.4))))
  expect_false(isTRUE(all.equal(lf(real, 0.0), lf(real, 0.4))))
})


test_that("estimate_rates accepts branching times with no topology", {
  skip_on_cran()
  phy  <- .mu0_tree(5L, 8L)
  brts <- as.numeric(sort(ape::branching.times(phy), decreasing = TRUE))
  expect_null(attr(brts, "parent_tip_start"))
  fit <- estimate_rates(brts, method = "mcem", model = "cr",
                        init_pars = c(0.5, 0.1),
                        control = list(lower_bound = c(0, 0),
                                       upper_bound = c(2, 1),
                                       num_trees = 10L, max_iter = 2L))
  expect_s3_class(fit, "emphasis_fit")
  expect_true(is.finite(fit$loglik))
})


# --------------------------------------------------------------------------- #
#  The R-side extraction                                                       #
# --------------------------------------------------------------------------- #

test_that(".observed_parent_tip_start is the parent node's forward time", {
  phy <- .mu0_tree(2L, 9L)
  h    <- ape::node.depth.edgelength(phy)
  root <- ape::Ntip(phy) + 1L
  int  <- setdiff(seq.int(root, max(phy$edge)), root)
  par  <- stats::setNames(phy$edge[, 1L], phy$edge[, 2L])
  want <- unname(h[par[as.character(int[order(h[int])])]])
  expect_equal(emphasis:::.observed_parent_tip_start(phy), want)
  # One per observed branching event; the crown split has no node of its own.
  expect_length(emphasis:::.observed_parent_tip_start(phy), ape::Nnode(phy) - 1L)
  # A child of the crown starts at the crown, forward time 0.
  expect_equal(min(emphasis:::.observed_parent_tip_start(phy)), 0)
})

test_that(".extract_brts carries the tip starts, and only for a tree", {
  phy <- .mu0_tree(2L, 9L)
  b   <- emphasis:::.extract_brts(phy)
  expect_length(emphasis:::.pts(b), length(b))
  expect_equal(utils::tail(emphasis:::.pts(b), 1L), -1)   # terminal marker
  expect_equal(emphasis:::.pts(as.numeric(b)), numeric(0))
  # Round-tripping a vector that already carries them keeps them.
  expect_equal(emphasis:::.pts(emphasis:::.extract_brts(b)), emphasis:::.pts(b))
  # A simulate_tree() result resolves to its extant tree.
  expect_equal(emphasis:::.pts(emphasis:::.extract_brts(list(tes = phy))),
               emphasis:::.pts(b))
})

test_that("a wrong-length parent_tip_start is refused", {
  phy  <- .mu0_tree(2L, 9L)
  brts <- emphasis:::.extract_brts(phy)
  expect_error(
    augment_trees(as.numeric(brts), c(0.3, 0, 0, 0, 0, 0, 0, 0),
                  1L, 100L, 100L, 1e6, 1L, model = c(0L, 0L, 0L), link = 0L,
                  rho = 1, parent_tip_start = c(0, 0)),
    "parent_tip_start"
  )
})

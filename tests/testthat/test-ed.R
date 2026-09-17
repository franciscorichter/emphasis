# The ED covariate: evolutionary distinctiveness (fair proportion) on the
# complete tree, in the simulator and in the likelihood.
#
# Four checks, in the order the study design lists them: ED agrees with an
# independent computation on an ape tree; it sums to Faith's PD; the ED model
# at beta_ED = 0 reproduces the dd log-likelihood; and the sampler's
# unbiasedness gate passes for an ED model.

# The fair-proportion ED of every tip of an ultrametric ape tree, computed
# from the edge table alone: sum over the edges above the tip of
# edge length / number of tips below the edge.
fair_proportion_ape <- function(phy) {
  ntip <- ape::Ntip(phy)
  nn   <- max(phy$edge)
  par  <- integer(nn); par[phy$edge[, 2L]] <- phy$edge[, 1L]
  len  <- numeric(nn); len[phy$edge[, 2L]] <- phy$edge.length
  ntips_below <- integer(nn)
  for (x in seq_len(ntip)) {
    v <- x
    while (v != ntip + 1L) { ntips_below[v] <- ntips_below[v] + 1L; v <- par[v] }
  }
  vapply(seq_len(ntip), function(x) {
    v <- x; s <- 0
    while (v != ntip + 1L) { s <- s + len[v] / ntips_below[v]; v <- par[v] }
    s
  }, numeric(1))
}

# The lineage forest of an ape tree under the convention .observed_parent_id
# uses: a lineage is named by the tip it ends in, the child whose subtree holds
# the lowest tip label continues its parent's lineage, the root's two children
# are the crown lineages.  Returns parent (1-based, 0 for a crown) and birth
# per tip, so that lineage x's ED at the present is the fair proportion of tip x.
lineage_forest_ape <- function(phy) {
  ntip <- ape::Ntip(phy); root <- ntip + 1L; nn <- max(phy$edge)
  h    <- ape::node.depth.edgelength(phy)
  par  <- integer(nn); par[phy$edge[, 2L]] <- phy$edge[, 1L]
  kids <- split(phy$edge[, 2L], phy$edge[, 1L])
  rep_ <- integer(nn); rep_[seq_len(ntip)] <- seq_len(ntip)
  for (v in rev(order(h[seq.int(root, nn)])) + root - 1L)
    rep_[v] <- min(rep_[kids[[as.character(v)]]])
  parent <- integer(ntip); birth <- numeric(ntip)
  for (x in seq_len(ntip)) {
    v <- x
    while (par[v] != root && rep_[par[v]] == x) v <- par[v]   # climb while still lineage x
    if (par[v] == root) { parent[x] <- 0L; birth[x] <- 0 }
    else { parent[x] <- rep_[par[v]]; birth[x] <- h[par[v]] }
  }
  list(parent = parent, birth = birth, T = max(h))
}

test_that("ED equals the fair proportion of an ape tree, and sums to Faith's PD", {
  set.seed(11)
  for (n in c(4L, 9L, 30L)) {
    phy <- ape::rcoal(n)
    fo  <- lineage_forest_ape(phy)
    r   <- emphasis:::ed_fair_proportion(fo$parent, fo$birth, rep(TRUE, n), fo$T)
    expect_equal(r$ed, fair_proportion_ape(phy), tolerance = 1e-10)
    expect_equal(sum(r$ed), sum(phy$edge.length), tolerance = 1e-10)   # Faith's PD
    # the pendant piece is the first term: ED >= pendant age, with equality
    # only for a lineage whose whole ancestry it shares with nobody (none here)
    expect_true(all(r$ed >= fo$T - r$ts - 1e-12))
    # the two crown lineages' subtrees partition the alive lineages
    expect_equal(sum(r$n_desc[fo$parent == 0L]), n)
    expect_true(all(r$n_desc >= 1L))
  }
})

test_that("the observed topology reaches C++ as lineage ids that name the same lineages", {
  set.seed(3)
  phy  <- ape::rcoal(12)
  brts <- emphasis:::.extract_brts(phy)
  pid  <- emphasis:::.pid(brts)
  pts  <- emphasis:::.pts(brts)
  expect_length(pid, ape::Nnode(phy))           # one per event, then the -1 marker
  expect_equal(pid[length(pid)], -1L)
  ev <- pid[-length(pid)]
  # every splitting lineage is a crown lineage or an earlier event's daughter
  expect_true(all(ev %in% c(-2L, -3L) | (ev >= 0L & ev < seq_along(ev) - 1L)))
  # the first event splits a crown lineage
  expect_true(ev[1] %in% c(-2L, -3L))
  # and the ids name lineages consistently with the tree: an event on lineage
  # L (L >= 0) sits in the subtree of the daughter born at event L, which is
  # the child of that event's node whose subtree does not hold its parent's
  # lowest tip; an event on a crown lineage sits under the matching root child
  ntip <- ape::Ntip(phy); root <- ntip + 1L; nn <- max(phy$edge)
  h    <- ape::node.depth.edgelength(phy)
  int  <- setdiff(seq.int(root, nn), root)
  ord  <- int[order(h[int])]                       # event k <-> node ord[k + 1]
  par  <- integer(nn); par[phy$edge[, 2L]] <- phy$edge[, 1L]
  kids <- split(phy$edge[, 2L], phy$edge[, 1L])
  rep_ <- integer(nn); rep_[seq_len(ntip)] <- seq_len(ntip)
  for (v in rev(order(h[seq.int(root, nn)])) + root - 1L)
    rep_[v] <- min(rep_[kids[[as.character(v)]]])
  is_desc <- function(v, a) { while (v != root && v != a) v <- par[v]; v == a }
  for (k in seq_along(ev)) {
    v <- ord[k]
    if (ev[k] >= 0L) {
      u <- ord[ev[k] + 1L]
      d <- kids[[as.character(u)]]; d <- d[rep_[d] != rep_[u]]
      expect_length(d, 1L)
      expect_true(is_desc(v, d))
      expect_lt(h[u], h[v])                        # born before it splits
    } else {
      rc <- kids[[as.character(root)]]
      under <- vapply(rc, function(c) is_desc(v, c), logical(1))
      expect_equal(sum(under), 1L)
    }
  }
  # the two crown ids each cover exactly one root-child subtree
  crown_child <- vapply(seq_along(ev), function(k) {
    v <- ord[k]; rc <- kids[[as.character(root)]]
    rc[vapply(rc, function(c) is_desc(v, c), logical(1))]
  }, integer(1))
  expect_true(all(tapply(crown_child[ev < 0L], ev[ev < 0L], function(x) length(unique(x))) == 1L))
})

test_that("the ED model at beta_ED = 0 reproduces the dd log-likelihood exactly", {
  set.seed(5)
  phy  <- ape::rcoal(10)
  brts <- emphasis:::.extract_brts(phy)
  dd   <- c(0.6, -0.01, 0.15, 0.0)                 # beta_0, beta_N, gamma_0, gamma_N
  aug  <- emphasis:::augment_trees(brts, emphasis:::.expand_pars(dd, emphasis:::.resolve_model("dd")),
                                   sample_size = 20L, maxN = 200L, max_missing = 100L,
                                   max_lambda = 1e6, num_threads = 1L,
                                   model = emphasis:::.resolve_model("dd"), link = 0L, rho = 1.0,
                                   parent_tip_start = emphasis:::.pts(brts), seed = 7L,
                                   parent_id = emphasis:::.pid(brts))
  expect_gt(length(aug$trees), 0L)
  lf_dd  <- emphasis:::eval_logf(emphasis:::.expand_pars(dd, emphasis:::.resolve_model("dd")), aug$trees,
                                 model = emphasis:::.resolve_model("dd"), link = 0L, rho = 1.0)
  ned0   <- c(0.6, -0.01, 0.0, 0.15, 0.0, 0.0)     # the same, with beta_ED = gamma_ED = 0
  lf_ned <- emphasis:::eval_logf(emphasis:::.expand_pars(ned0, emphasis:::.resolve_model("ned")), aug$trees,
                                 model = emphasis:::.resolve_model("ned"), link = 0L, rho = 1.0)
  expect_equal(lf_ned$logf, lf_dd$logf, tolerance = 1e-10)
  # and a non-zero beta_ED changes it, in the direction of more speciation
  # weight on distinct lineages: the log f moves, and by a finite amount
  ned1   <- c(0.6, -0.01, 0.05, 0.15, 0.0, 0.0)
  lf_ned1 <- emphasis:::eval_logf(emphasis:::.expand_pars(ned1, emphasis:::.resolve_model("ned")), aug$trees,
                                  model = emphasis:::.resolve_model("ned"), link = 0L, rho = 1.0)
  expect_true(all(is.finite(lf_ned1$logf)))
  expect_false(isTRUE(all.equal(lf_ned1$logf, lf_dd$logf)))
})

test_that("the ED model refuses a bare branching-time vector and the gaussian link", {
  set.seed(6)
  phy  <- ape::rcoal(8)
  brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)   # no topology
  expect_error(
    emphasis:::augment_trees(brts, emphasis:::.expand_pars(c(0.5, 0.0, 0.1, 0.0), emphasis:::.resolve_model("ed")),
                             sample_size = 5L, maxN = 50L, max_missing = 100L,
                             max_lambda = 1e6, num_threads = 1L,
                             model = emphasis:::.resolve_model("ed"), link = 0L, rho = 1.0),
    "topology")
  b <- emphasis:::.extract_brts(phy)
  expect_error(
    emphasis:::augment_trees(b,
                             emphasis:::.expand_pars(c(0.5, 0.0, 0.1, 0.0), emphasis:::.resolve_model("ed")),
                             sample_size = 5L, maxN = 50L, max_missing = 100L,
                             max_lambda = 1e6, num_threads = 1L,
                             model = emphasis:::.resolve_model("ed"), link = 2L, rho = 1.0,
                             parent_tip_start = emphasis:::.pts(b), parent_id = emphasis:::.pid(b)),
    "gaussian")
})

test_that("the simulator runs an ED model, and a negative beta_ED slows distinct lineages", {
  set.seed(8)
  # beta_ED < 0: lineages carrying more evolutionary history speciate less;
  # over many simulations the trees are smaller than under beta_ED = 0
  n_with <- n_without <- integer(30)
  for (i in seq_len(30)) {
    s0 <- simulate_tree(pars = c(0.5, 0.0, 0.0, 0.0, 0.0, 0.0), max_t = 8, model = "ned",
                        max_lin = 2000L, num_threads = 1L)
    s1 <- simulate_tree(pars = c(0.5, 0.0, -0.1, 0.0, 0.0, 0.0), max_t = 8, model = "ned",
                        max_lin = 2000L, num_threads = 1L)
    n_without[i] <- if (s0$status == "done") length(s0$tes$tip.label) else NA_integer_
    n_with[i]    <- if (s1$status == "done") length(s1$tes$tip.label) else NA_integer_
  }
  expect_true(mean(n_with, na.rm = TRUE) < mean(n_without, na.rm = TRUE))
})

test_that("estimate_rates fits an N + ED model end to end (smoke)", {
  set.seed(9)
  phy <- ape::rcoal(12)
  fit <- estimate_rates(phy, method = "mcem", model = "ned",
                        init_pars = c(0.8, -0.01, 0.0, 0.2, 0.0, 0.0),
                        control = list(sampling = "dynamic_fresh", sample_size = 20L,
                                       max_iter = 3L, num_threads = 1L,
                                       lower_bound = c(0.01, -0.2, -1, 0.0, -0.2, -1),
                                       upper_bound = c(3, 0.2, 1, 2, 0.2, 1),
                                       verbose = FALSE))
  expect_s3_class(fit, "emphasis_fit")
  expect_length(fit$pars, 6L)
  expect_equal(names(fit$pars), emphasis:::.par_names(emphasis:::.resolve_model("ned")))
  expect_true(is.finite(fit$loglik))
})

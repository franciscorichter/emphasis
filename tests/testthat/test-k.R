# The K covariate: the number of speciation events on the path from the crown
# to the lineage.  Constant between events, one more on both daughters at a
# split, hidden splits counted.  Per lineage like ED, so a K model takes the
# per-lineage likelihood and the per-lineage thinning proposal.

.node_depth_minus_crown <- function(phy) {
  par <- integer(max(phy$edge)); par[phy$edge[, 2]] <- phy$edge[, 1]
  root <- ape::Ntip(phy) + 1L
  vapply(seq_len(ape::Ntip(phy)), function(i) {
    d <- 0L; v <- i
    while (v != root) { v <- par[v]; d <- d + 1L }
    d - 1L
  }, 1L)
}

test_that("K resolves, expands and is named like the other covariates", {
  expect_equal(emphasis:::.resolve_model("k"),  c(0L, 0L, 0L, 0L, 1L))
  expect_equal(emphasis:::.resolve_model("nk"), c(1L, 0L, 0L, 0L, 1L))
  expect_equal(emphasis:::.resolve_model(~ N + K), c(1L, 0L, 0L, 0L, 1L))
  expect_equal(emphasis:::.resolve_model(~ K + ED), c(0L, 0L, 0L, 1L, 1L))
  mb <- emphasis:::.resolve_model("nk")
  compact <- c(0.5, -0.01, 0.2, 0.1, 0.02, -0.05)
  full <- emphasis:::.expand_pars(compact, mb)
  expect_length(full, 12L)
  expect_equal(full[c(1, 2, 5, 6)], c(0.5, -0.01, 0.1, 0.02))
  expect_equal(full[c(11, 12)], c(0.2, -0.05))
  expect_equal(full[c(3, 4, 7, 8, 9, 10)], rep(0, 6))
  expect_equal(emphasis:::.contract_pars(full, mb), compact)
  expect_equal(emphasis:::.par_names(mb), c("beta_0", "beta_N", "beta_K", "gamma_0", "gamma_N", "gamma_K"))
  expect_equal(emphasis:::.model_label(mb), "N + K")
  # a model without ED or K still expands to the 8-slot layout
  expect_length(emphasis:::.expand_pars(c(0.5, 0.02, 0.1, -0.01), c(1L, 0L, 0L, 0L, 0L)), 8L)
})

test_that("the simulator's K is the node depth of the tree it grew", {
  set.seed(11)
  s <- simulate_tree(pars = c(log(0.4), -0.004, 0.12, log(0.05), 0, 0), max_t = 10,
                     model = "nk", link = "exponential", rho = 1, max_tries = 5)
  skip_if(!identical(s$status, "done"), "no tree")
  tab <- lineage_table(s$tes)
  last <- tab[tab$seg == max(tab$seg), ]
  expect_equal(sort(as.integer(last$K)), sort(.node_depth_minus_crown(s$tes)))
  # K is the same under either link and needs no coefficient to exist
  tab2 <- lineage_table(ape::rphylo(20L, 0.6, 0.1))
  expect_true(all(tab2$K >= 0)); expect_true(all(tab2$K == round(tab2$K)))
})

test_that("the K likelihood is the Poisson form of the lineage table on a fully observed tree", {
  set.seed(4)
  phy <- ape::rphylo(25L, 0.7, 0)
  tab <- lineage_table(phy)
  fr  <- emphasis:::.observed_frame(phy)
  mb  <- emphasis:::.resolve_model("nk")
  dt  <- tab$t1 - tab$t0
  # exponential link
  p <- c(log(0.4), -0.004, 0.12, log(0.05), 0.001, -0.02)
  lf <- emphasis:::eval_logf(emphasis:::.expand_pars(p, mb), list(fr), model = mb, link = 1L, rho = 1)$logf
  eta <- p[1] + p[2] * tab$N + p[3] * tab$K
  etm <- p[4] + p[5] * tab$N + p[6] * tab$K
  expect_equal(lf, sum(eta[tab$event == 1L]) - sum(exp(eta) * dt) - sum(exp(etm) * dt), tolerance = 1e-10)
  # linear link, with a lineage clipped at zero somewhere
  p2 <- c(0.4, -0.002, -0.05, 0.05, 0, 0.01)
  lf2 <- emphasis:::eval_logf(emphasis:::.expand_pars(p2, mb), list(fr), model = mb, link = 0L, rho = 1)$logf
  eta2 <- pmax(0, p2[1] + p2[2] * tab$N + p2[3] * tab$K)
  etm2 <- pmax(0, p2[4] + p2[5] * tab$N + p2[6] * tab$K)
  expect_true(any(eta2 == 0))
  expect_equal(lf2, sum(log(eta2[tab$event == 1L])) - sum(eta2 * dt) - sum(etm2 * dt), tolerance = 1e-10)
  # K off collapses onto dd
  p0 <- c(log(0.4), -0.004, 0, log(0.05), 0, 0)
  a <- emphasis:::eval_logf(emphasis:::.expand_pars(p0, mb), list(fr), model = mb, link = 1L, rho = 1)$logf
  b <- emphasis:::eval_logf(emphasis:::.expand_pars(p0[c(1, 2, 4, 5)], emphasis:::.resolve_model("dd")),
                            list(fr), model = emphasis:::.resolve_model("dd"), link = 1L, rho = 1)$logf
  expect_equal(a, b, tolerance = 1e-12)
})

test_that("K models augment on the thinning proposal and are refused by the surrogate with a reason", {
  skip_on_cran()
  set.seed(5)
  phy <- ape::rphylo(30L, 0.6, 0.15)
  mb  <- emphasis:::.resolve_model("nk")
  p   <- c(log(0.5), -0.004, 0.1, log(0.1), 0, 0)
  a <- emphasis:::.augment_tree_internal(phy, pars = p, model_bin = mb, sample_size = 60L,
                                         link = 1L, rho = 1, maxN = 50000L)
  expect_gt(length(a$trees), 0L)
  lw <- a$logf - a$logg
  expect_true(all(is.finite(lw)))
  # the density the sampler charges is the one sampling_prob recomputes on the frame
  ev <- emphasis:::eval_logf(emphasis:::.expand_pars(p, mb), a$trees, model = mb, link = 1L, rho = 1)
  expect_equal(ev$logf, a$logf, tolerance = 1e-8)
  expect_false(emphasis:::.bdi_supported(mb, 1L, 1))
  expect_match(emphasis:::.bdi_unsupported_reason(mb, 1L, 1), "K-dependent")
  # the path offers K and runs
  cp <- covariate_path(a$trees, lw, n_trees = 20, covariates = c("N", "D", "ED", "K"), rate = "speciation")
  expect_true("K" %in% names(cp$speciation$path))
  expect_length(cp$speciation$coef, 5L)
  # a bare branching-time vector cannot carry K
  expect_error(emphasis:::augment_trees(as.numeric(emphasis:::.extract_brts(phy)), emphasis:::.expand_pars(p, mb),
                                        5L, 2000L, 100L, 1e6, 1L, model = mb, link = 1L, rho = 1),
               "topology")
})

# The lineage table is the design the likelihood is built on: rebuilding the
# exponential-link log-likelihood of the full model N + D + ED from it, with
# exact per-segment integrals, must reproduce eval_logf on augmented trees.
# That checks the table's alive set, pendant starts and ED values against the
# C++ likelihood, and the likelihood's ED path against the table's D starts.

.logf_from_table <- function(tab, df, p10) {
  # p10: beta_0 beta_N beta_M beta_D gamma_0 gamma_N gamma_M gamma_D beta_ED gamma_ED
  b <- p10[c(1, 2, 4, 9)]; g <- p10[c(5, 6, 8, 10)]
  # per row: eta(u) = a + s u on [t0, t1] with
  #   a = b0 + bN N + bD(-ts - M) + bED(ed0 - t0),  s = bD + bED
  int_exp <- function(a, s, t0, t1) if (abs(s) < 1e-12) exp(a) * (t1 - t0) else
                                      exp(a) * (exp(s * t1) - exp(s * t0)) / s
  a_l <- b[1] + b[2] * tab$N + b[3] * (-tab$ts - tab$M) + b[4] * (tab$ed0 - tab$t0)
  a_m <- g[1] + g[2] * tab$N + g[3] * (-tab$ts - tab$M) + g[4] * (tab$ed0 - tab$t0)
  inte <- sum(int_exp(a_l, b[3] + b[4], tab$t0, tab$t1)) + sum(int_exp(a_m, g[3] + g[4], tab$t0, tab$t1))
  ev <- tab[tab$event != 0L, ]
  eta_ev <- function(cf) cf[1] + cf[2] * ev$N + cf[3] * ((ev$t1 - ev$ts) - ev$M) + cf[4] * (ev$ed0 + ev$t1 - ev$t0)
  events <- sum(ifelse(ev$event == 1L, eta_ev(b), eta_ev(g)))
  events - inte
}

test_that("the lineage table reproduces the full model's log-likelihood under the exponential link", {
  skip_on_cran()
  set.seed(9)
  phy  <- ape::rphylo(14L, 0.7, 0.2)
  brts <- emphasis:::.extract_brts(phy)
  p10  <- c(log(0.6), -0.004, 0, 0.05, log(0.15), 0, 0, -0.03, -0.02, 0.01)
  aug <- emphasis:::augment_trees(as.numeric(brts), p10, 12L, 20000L, 300L, 1e6, 1L,
                                  model = c(1L, 0L, 1L, 1L), link = 1L, rho = 1,
                                  parent_tip_start = emphasis:::.pts(brts), seed = 5L,
                                  parent_id = emphasis:::.pid(brts))
  skip_if(length(aug$trees) == 0L, "augmentation drew no tree")
  got <- emphasis:::eval_logf(p10, aug$trees, model = c(1L, 0L, 1L, 1L), link = 1L, rho = 1)$logf
  ref <- vapply(aug$trees, function(df) .logf_from_table(lineage_table(df), df, p10), 0)
  expect_true(all(is.finite(got)))
  expect_equal(got, ref, tolerance = 1e-8)
  # the table's alive count is the node's n on every segment
  tab <- lineage_table(aug$trees[[1L]])
  n_seg <- tapply(tab$lineage, tab$seg, length)
  expect_equal(as.integer(n_seg), as.integer(aug$trees[[1L]]$n[as.integer(names(n_seg)) + 1L]))
})

test_that("the covariate path runs on a fully observed tree and returns an entry order", {
  skip_on_cran()
  skip_if_not_installed("dglars")
  set.seed(3)
  phy <- ape::rphylo(30L, 0.6, 0)
  cp <- covariate_path(phy, grid = 200L, rate = "speciation")
  expect_named(cp, "speciation")
  s <- cp$speciation
  expect_equal(s$n_events, sum(lineage_table(phy)$event == 1L))   # every split lands in a cell
  expect_true(all(names(s$entry_step) %in% c("N", "D", "ED")))
  expect_true(is.finite(s$coef[[1L]]))
  expect_equal(names(s$coef)[1L], "(Intercept)")
})

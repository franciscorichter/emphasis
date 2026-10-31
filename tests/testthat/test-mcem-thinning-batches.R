# The trees an E-step hands back (em_cpp(copy_trees = TRUE)) carry their topology, so that the
# batch M-steps of the "mc_error" rule can score a lineage-level model on them.  Before the fix,
# rcpp_mcem.cpp returned brts / n / t_ext only; every batch M-step of a D, ED or K model then
# failed on "the ED covariate needs the tree's topology", mc_se() returned NULL, mc_z was NA on
# every iteration, the draws never grew and no thinning fit of those models ever met the rule
# (E46's first tier: 311 of 395 fits on the 1000-iteration cap).
# num_threads = 1; a stored tree, so that the test does not depend on a simulation surviving.

sim_tree <- function() {
  for (s in 1:50) {
    set.seed(100 + s)
    r <- tryCatch(suppressWarnings(simulate_tree(pars = c(-0.7, -0.01, -0.2, -3, 0, 0), max_t = 10,
                  model = "nk", link = "exponential", max_lin = 400, max_tries = 1L, num_threads = 1L)),
                  error = function(e) NULL)
    if (!is.null(r) && identical(r$status, "done") && !is.null(r$tes) && ape::Ntip(r$tes) >= 15 && ape::Ntip(r$tes) <= 80) return(r$tes)
  }
  skip("no suitable nk tree in 50 seeds")
}

test_that("the trees an E-step returns carry the lineage ids the covariate tables need", {
  phy <- sim_tree(); brts <- emphasis:::.extract_brts(phy)
  mb <- emphasis:::.resolve_model("nk"); p_full <- emphasis:::.expand_pars(c(-0.7, -0.01, -0.2, -3, 0, 0), mb)
  res <- emphasis:::em_cpp(brts = brts, init_pars = p_full, sample_size = 20L, maxN = 2000L, max_missing = 1e4,
           max_lambda = 1e6, lower_bound = p_full - 1, upper_bound = p_full + 1, xtol_rel = 1e-3, num_threads = 1L,
           copy_trees = TRUE, rconditional = NULL,
           model = as.integer(mb), link = 1L, rho = 1, parent_tip_start = emphasis:::.pts(brts), parent_id = emphasis:::.pid(brts))
  expect_true(is.list(res$trees) && length(res$trees) >= 1)
  expect_true(all(c("id", "parent_id", "tip_start", "focal_tip_start", "clade", "pd") %in% names(res$trees[[1]])))
})

test_that("mc_error computes the Monte Carlo error of a lineage-level thinning fit", {
  phy <- sim_tree()
  f <- suppressWarnings(suppressMessages(estimate_rates(phy, method = "mcem", model = "nk", link = "exponential",
         init_pars = c(-0.7, -0.01, -0.2, -3, 0, 0),
         control = list(lower_bound = c(-4, -0.1, -1, -7, -0.1, -1), upper_bound = c(1, 0.1, 1, 1, 0.1, 1),
                        sample_size = 40, max_draws = 80, max_iter = 2, stop_rule = "mc_error",
                        sampling = "dynamic_fresh", num_threads = 1L, verbose = FALSE))))
  tr <- f$details$mcem
  expect_true(any(is.finite(tr$mc_z)))
})

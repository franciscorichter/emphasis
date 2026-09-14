# The importance-sampling E-step on a fixed tree (audit ids H94, H96, H104).
#
# This file held one block, skipped unconditionally behind a skip_on_cran(),
# that called mc_loglik() -- a function that exists in no version of the
# package.  What it meant to exercise is the path from a branching-time vector
# to an fhat: augment_trees draws the augmented trees and records logf and
# logg, and .is_fhat aggregates them.  The block is written against that path.
#
# The tree is a fixed branching-time vector rather than an ape::rphylo draw, so
# the observed data are the same on every run; the augmentation itself is
# clock-seeded and its assertions are structural or tolerances.

# Fixed 8-tip branching times, crown age 4, decreasing.
brts8 <- c(4, 3.4, 2.9, 2.3, 1.8, 1.2, 0.5)

# dd in the 8-parameter layout c(beta_0, beta_N, beta_M, beta_D,
#                                gamma_0, gamma_N, gamma_M, gamma_D)
pars_dd <- c(0.5, -0.01, 0, 0, 0.1, 0, 0, 0)
dd_bin  <- c(1L, 0L, 0L)

test_that("augment_trees + .is_fhat give a finite fhat on a small tree", {
  aug <- emphasis:::augment_trees(
    brts        = brts8,
    pars        = pars_dd,
    sample_size = 30L,
    maxN        = 2000L,
    max_missing = 500L,
    max_lambda  = 500,
    num_threads = 1L,
    model       = dd_bin,
    link        = 0L
  )

  expect_type(aug, "list")
  expect_length(aug$trees, 30L)
  expect_true(all(is.finite(aug$logf)))
  expect_true(all(is.finite(aug$logg)))

  fhat <- emphasis:::.is_fhat(aug$logf, aug$logg,
                              n_zero_weight = aug$rejected_zero_weights)
  expect_true(is.finite(fhat))
  expect_lt(fhat, 0)

  # ESS is bounded by the number of draws
  ess <- emphasis:::.ess_from_lw(aug$logf - aug$logg)
  expect_gt(ess, 0)
  expect_lte(ess, 30)
})

test_that("max_missing caps the augmentation and fills the overrun channel", {
  # With at most one extinct lineage per tree, about half the draws on this
  # tree overrun, so the overrun counter must move and every returned tree
  # must respect the cap.  That counter is what diagnose_mcem reports and what
  # the adaptive limit escalation reads.
  aug <- emphasis:::augment_trees(
    brts        = brts8,
    pars        = pars_dd,
    sample_size = 20L,
    maxN        = 20000L,
    max_missing = 1L,
    max_lambda  = 500,
    num_threads = 1L,
    model       = dd_bin,
    link        = 0L
  )

  expect_length(aug$trees, 20L)
  expect_gt(aug$rejected_overruns, 0L)
  for (tr in aug$trees) {
    # an extinct lineage contributes one death row, marked by t_ext == 0
    expect_lte(sum(tr$t_ext == 0), 1L)
  }
})

test_that("fhat is invariant to the order the draws come back in", {
  aug <- emphasis:::augment_trees(brts8, pars_dd, sample_size = 25L,
                                  maxN = 2000L, max_missing = 500L,
                                  max_lambda = 500, num_threads = 1L,
                                  model = dd_bin, link = 0L)
  ord <- rev(seq_along(aug$logf))
  expect_equal(emphasis:::.is_fhat(aug$logf, aug$logg),
               emphasis:::.is_fhat(aug$logf[ord], aug$logg[ord]),
               tolerance = 1e-12)
  # re-scoring the permuted trees returns the permuted densities
  ev <- emphasis:::eval_logf(pars_dd, aug$trees[ord], model = dd_bin, link = 0L)
  expect_identical(ev$logf, aug$logf[ord])
  expect_identical(ev$logg, aug$logg[ord])
})

# The M-step's first simplex is sized from the box (src/M_step.cpp), not from
# the distance to the nearer bound, which is what NLopt does when no initial
# step is set.  A start a hair inside a bound -- a nested fit's estimate
# clipped to the box by the warm start -- then gets a step of that hair and
# sbplx returns the start unchanged with XTOL_REACHED, at every EM iteration:
# on the ED simulation arm's main tier that froze four of ten fits of one cell
# at their starting beta_ED = 0.

e_step_at <- function(brts, pars10, mb, seed) {
  aug <- emphasis:::augment_trees(brts, pars10, sample_size = 60L, maxN = 3000L, max_missing = 5000L,
                                  max_lambda = 1e6, num_threads = 1L, model = mb, link = 0L, rho = 1.0,
                                  parent_tip_start = emphasis:::.pts(brts), seed = seed,
                                  parent_id = emphasis:::.pid(brts))
  lf <- emphasis:::eval_logf(pars10, aug$trees, model = mb, link = 0L, rho = 1.0)
  lw <- lf$logf - lf$logg
  w  <- exp(lw - max(lw)); w[!is.finite(lw)] <- 0
  list(trees = aug$trees, weights = w, rejected = 0L, rejected_overruns = 0L, rejected_lambda = 0L,
       rejected_zero_weights = 0L, time = 0, fhat = 0)
}

test_that("an M-step started a hair inside a bound moves off it", {
  set.seed(31)
  phy  <- ape::rcoal(20)
  brts <- emphasis:::.extract_brts(phy)
  mb   <- emphasis:::.resolve_model("dd")
  lo   <- c(0.05, -0.03, 0.0, -0.03)
  up   <- c(3.00,  0.03, 1.5,  0.03)
  ex   <- function(p) emphasis:::.expand_pars(p, mb)
  # beta_N one unit in the last place above its lower bound (the start that
  # froze on the ED arm sat 1.7e-18 above a bound of -0.0109), the rest interior
  start <- c(0.9, lo[2] * (1 - .Machine$double.eps), 0.2, 0.0)
  expect_gt(start[2], lo[2])
  es <- e_step_at(brts, ex(start), mb, seed = 5L)
  expect_gt(sum(es$weights > 0), 10L)
  m <- emphasis:::m_cpp(es, ex(start), "rpd5c", ex(lo), ex(up), 1e-3, 1L, model = mb, link = 0L, rho = 1.0)
  est <- emphasis:::.contract_pars(m$estimates, mb)
  expect_true(any(abs(est - start) > 1e-8))
  # and it lands where a start well inside the box lands, to the tolerance
  start2 <- start; start2[2] <- lo[2] + 0.2 * (up[2] - lo[2])
  m2 <- emphasis:::m_cpp(es, ex(start2), "rpd5c", ex(lo), ex(up), 1e-3, 1L, model = mb, link = 0L, rho = 1.0)
  est2 <- emphasis:::.contract_pars(m2$estimates, mb)
  expect_equal(est, est2, tolerance = 5e-2)
})

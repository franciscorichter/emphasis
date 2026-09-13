# Tests for the thinning E-step (src/E_step.cpp, E_step()).
# Pins wave-1 item 1.7: audit ids H1, H16 (E-step half), H62, H68.
#   - a completed augmentation is accepted iff its log-weight is finite
#     (no absolute cut at log_w = -745.13 from exp() underflow);
#   - log_w = -Inf counts as a zero-weight tree, +Inf/NaN separately;
#   - with num_threads > 1 the sampler stops at exactly sample_size trees;
#   - fhat divides by trees.size() + rejected_zero_weights.

# 60-tip constant-rate tree (TreeSim::sim.bd.taxa(60, 1, 0.2, 0.05), seed 60),
# branching times crown-first.  At lambda = 0.2, mu = 0.05 its log-weights sit
# near -150; rescaling time by s shifts every log-weight by -(n - 2) * log(s).
b60 <- c(20.2571, 18.9106, 18.8223, 14.7584, 13.9755, 11.9092, 11.3709,
         11.2664, 11.0612, 9.3217, 9.2315, 8.7633, 7.4891, 6.9616, 6.339,
         6.2372, 6.2025, 5.3442, 5.3046, 5.1328, 4.7774, 4.7511, 3.8254,
         3.8123, 3.3753, 3.2648, 2.7234, 2.7159, 2.6661, 2.3809, 2.2352,
         2.1715, 2.0224, 1.9892, 1.9652, 1.9154, 1.7791, 1.7301, 1.6915,
         1.529, 1.4641, 1.3952, 1.3293, 1.0584, 1.0575, 1.0053, 0.8894,
         0.8093, 0.7871, 0.6553, 0.6152, 0.4845, 0.4118, 0.3098, 0.2323,
         0.2081, 0.1329, 0.0911, 0.0353)

# 4-tip constant-rate tree: fast augmentation, used for the thread tests.
b4 <- c(4, 2.5, 1.2, 0.6)

cr    <- c(0L, 0L, 0L)
pars8 <- function(lambda, mu) c(lambda, 0, 0, 0, mu, 0, 0, 0)

aug <- function(brts, lambda, mu, N, maxN, num_threads = 1L) {
  augment_trees(brts, pars8(lambda, mu), as.integer(N), as.integer(maxN),
                max_missing = 10000L, max_lambda = 1e6,
                num_threads = as.integer(num_threads), model = cr, link = 0L,
                rho = 1.0)
}

# log-mean-exp of the IS weights over completed augmentations, as E_step() does
fhat_of <- function(r) {
  lw <- r$logf - r$logg
  m  <- max(lw)
  log(sum(exp(lw - m))) + m - log(length(lw) + r$rejected_zero_weights)
}

test_that("augment_trees returns num_trees and the non-finite counter", {
  r <- aug(b4, 0.5, 0.1, N = 20, maxN = 500)
  expect_true(all(c("num_trees", "rejected_nonfinite",
                    "rejected_zero_weights") %in% names(r)))
  expect_equal(r$num_trees, 20L)
  expect_equal(length(r$trees), 20L)
  expect_equal(length(r$logf), 20L)
  expect_equal(length(r$logg), 20L)
  expect_equal(r$rejected_nonfinite, 0L)
})

test_that("H1: a tree rescaled by 1e6 (log-weights near -950) is augmented", {
  N <- 100L; maxN <- 1000L; s <- 1e6
  n <- length(b60) + 1
  r <- aug(b60 * s, 0.2 / s, 0.05 / s, N, maxN)
  expect_equal(r$num_trees, N)
  expect_equal(length(r$logf), N)
  lw <- r$logf - r$logg
  expect_true(all(is.finite(lw)))
  # every accepted log-weight lies far below the former exp() underflow cut
  expect_true(all(lw < -745.13))
  expect_equal(r$rejected_zero_weights, 0L)
  expect_equal(r$rejected_nonfinite, 0L)
})

test_that("H1: fhat shifts by -(n - 2) * log(s) under time rescaling, within IS noise", {
  N <- 100L; maxN <- 1000L; s <- 1e6
  n <- length(b60) + 1
  reps <- 4L
  f1 <- replicate(reps, fhat_of(aug(b60, 0.2, 0.05, N, maxN)))
  fs <- replicate(reps, fhat_of(aug(b60 * s, 0.2 / s, 0.05 / s, N, maxN)))
  expect_true(all(is.finite(f1)))
  expect_true(all(is.finite(fs)))
  observed  <- mean(fs) - mean(f1)
  predicted <- -(n - 2) * log(s)
  # IS sd of fhat at N = 100 on this tree is ~0.3; the tolerance covers
  # the difference of two 4-replicate means with margin.
  expect_lt(abs(observed - predicted), 1.5)
})

test_that("zero-probability augmentations are counted as zero weights and the E-step errors", {
  # lambda = 0 under the linear link: log f = -Inf at every speciation node
  err <- tryCatch(aug(b4, 0, 0.1, N = 10, maxN = 40), error = function(e) e)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "40 zero weights")
  expect_match(conditionMessage(err), "0 non-finite weights")
  expect_match(conditionMessage(err), "Trees so far: 0")
})

test_that("em_cpp fhat is the log-mean-exp over trees.size() + rejected_zero_weights", {
  lb <- c(0.01, 0, 0, 0, 0.001, 0, 0, 0)
  ub <- c(5,    0, 0, 0, 5,     0, 0, 0)
  e <- em_cpp(brts = b4, init_pars = pars8(0.5, 0.1), sample_size = 50L,
              maxN = 5000L, max_missing = 1000L, max_lambda = 1e6,
              lower_bound = lb, upper_bound = ub, xtol_rel = 1e-3,
              num_threads = 1L, copy_trees = FALSE, model = cr, link = 0L,
              rho = 1, rconditional = NULL)
  expect_equal(e$trees, 50L)
  expect_equal(e$num_trees, 50L)
  expect_equal(length(e$logf), 50L)
  expect_equal(e$fhat, fhat_of(e), tolerance = 1e-10)
})

test_that("H62/H68: 50 replicate calls with num_threads = 8 all return exactly sample_size trees", {
  N <- 50L; maxN <- 5000L
  one <- function(nt) {
    r <- aug(b4, 0.5, 0.1, N, maxN, num_threads = nt)
    c(M = length(r$logf), num_trees = r$num_trees, fhat = fhat_of(r),
      rejected = r$rejected + r$rejected_overruns + r$rejected_lambda +
                 r$rejected_zero_weights + r$rejected_nonfinite)
  }
  m1 <- t(replicate(20, one(1L)))
  m8 <- t(replicate(50, one(8L)))
  expect_true(all(m1[, "M"] == N))
  expect_true(all(m8[, "M"] == N))
  expect_true(all(m8[, "num_trees"] == N))
  # fhat at 8 threads sits within the single-thread spread: the pre-fix bias
  # was log(M / N) >= 2.5 nats against a single-thread sd of ~0.05.
  sd1 <- sd(m1[, "fhat"])
  expect_lt(abs(mean(m8[, "fhat"]) - mean(m1[, "fhat"])), 4 * sd1 + 0.05)
  expect_lt(max(abs(m8[, "fhat"] - mean(m1[, "fhat"]))), 6 * sd1 + 0.1)
})

test_that("H62: em_cpp with num_threads = 8 returns sample_size trees and an unbiased fhat", {
  lb <- c(0.01, 0, 0, 0, 0.001, 0, 0, 0)
  ub <- c(5,    0, 0, 0, 5,     0, 0, 0)
  em1 <- function(nt) {
    r <- em_cpp(brts = b4, init_pars = pars8(0.5, 0.1), sample_size = 50L,
                maxN = 5000L, max_missing = 1000L, max_lambda = 1e6,
                lower_bound = lb, upper_bound = ub, xtol_rel = 1e-3,
                num_threads = nt, copy_trees = FALSE, model = cr, link = 0L,
                rho = 1, rconditional = NULL)
    c(trees = r$trees, num_trees = r$num_trees, fhat = r$fhat)
  }
  e1 <- t(replicate(10, em1(1L)))
  e8 <- t(replicate(20, em1(8L)))
  expect_true(all(e1[, "trees"] == 50L))
  expect_true(all(e8[, "trees"] == 50L))
  expect_true(all(e8[, "num_trees"] == 50L))
  expect_lt(abs(mean(e8[, "fhat"]) - mean(e1[, "fhat"])), 4 * sd(e1[, "fhat"]) + 0.05)
})

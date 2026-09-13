# Tests for log_sum (inst/include/model_helpers.hpp), the running sum of
# log(lambda) over speciation nodes inside model::loglik.
#
# Rule pinned here (audit H11(a)): a speciation node whose rate is <= 0 gives
# logf = -Inf whatever its position in the tree and whatever the running
# product/sum held at that point.  Before the fix the sign of a non-finite
# result was taken from the internal running sum, so a zero rate on the last
# scored node returned +Inf.
#
# Trees are hand-built in the C++ convention used by eval_logf: forward times,
# crown at 0 (not a node), observed speciation nodes at brts 1..k with
# n = 2..k+1, closing sentinel at tp.  Model dd (model_bin c(1,0,0)), linear
# link, so lambda(N) = max(0, beta0 + betaN * N).

mk_tree <- function(k, tp = k + 1) {
  data.frame(brts = c(seq_len(k), tp), n = c(2:(k + 1), k + 2), t_ext = 1e11,
             pd = 0, tip_start = 0, id = c(seq_len(k) - 1L, -1L), parent_id = -1L)
}
lam <- function(p8, N) pmax(0, p8[1] + p8[2] * N)
mu  <- function(p8, N) pmax(0, p8[5] + p8[6] * N)
hand_loglik <- function(p8, tr) {
  # sum log lambda over non-closing nodes - sum dt * n * (lambda + mu), rates at the end node
  dt <- diff(c(0, tr$brts))
  N  <- tr$n
  sum(log(lam(p8, N[-nrow(tr)]))) - sum(dt * N * (lam(p8, N) + mu(p8, N)))
}
ev <- function(p8, tr) {
  emphasis:::eval_logf(p8, list(tr), model = c(1L, 0L, 0L), link = 0L, rho = 1)$logf
}
p8 <- function(b0, bN, g0 = 0.1) c(b0, bN, 0, 0, g0, 0, 0, 0)

test_that("A0: finite tree agrees with the hand formula (tree builder is right)", {
  tr <- mk_tree(3); p <- p8(0.5, -0.05)
  expect_equal(ev(p, tr), hand_loglik(p, tr), tolerance = 1e-10)
})

test_that("A1: single node with lambda(2) = 0 on the last += gives -Inf", {
  tr <- mk_tree(1); p <- p8(0.2, -0.1)
  expect_identical(lam(p, 2), 0)
  expect_identical(ev(p, tr), -Inf)
})

test_that("A2: lambda(2) = 0 first, lambda(3) = 0.1 last gives -Inf", {
  tr <- mk_tree(2); p <- p8(-0.2, 0.1)
  expect_identical(lam(p, 2), 0)
  expect_gt(lam(p, 3), 0)
  expect_identical(ev(p, tr), -Inf)
})

test_that("A3: lambda(2) = 0.1 first, lambda(3) = 0 last gives -Inf", {
  tr <- mk_tree(2); p <- p8(0.3, -0.1)
  expect_gt(lam(p, 2), 0)
  expect_identical(lam(p, 3), 0)
  expect_identical(ev(p, tr), -Inf)
})

test_that("A4: zero on the last node after the running product crossed the upper threshold gives -Inf", {
  # lambda(N) = 104 - 4 N: lambda(2) = 96 ... lambda(25) = 4, lambda(26) = 0.
  # The product over 24 rates in [4, 96] exceeds 1e21, so sum_ > 0 when the zero arrives.
  tr <- mk_tree(25); p <- c(26 * 4, -4, 0, 0, 0.1, 0, 0, 0)
  expect_identical(lam(p, 26), 0)
  expect_gt(sum(log(lam(p, 2:25))), log(1e21))
  expect_identical(ev(p, tr), -Inf)
})

test_that("A5: zero on the last node after the running product crossed the lower threshold gives -Inf", {
  # lambda(N) = (31 - N) / 128, exact binary fractions so that beta0 + betaN * N is
  # exactly 0 at N = 31 even under FMA contraction.  The product over 29 rates in
  # [1/128, 29/128] falls below 1e-19, so sum_ < 0 when the zero arrives.
  tr <- mk_tree(30); p <- c(31 / 128, -1 / 128, 0, 0, 0.1, 0, 0, 0)
  expect_identical(lam(p, 31), 0)
  expect_lt(sum(log(lam(p, 2:30))), log(1e-19))
  expect_identical(ev(p, tr), -Inf)
})

test_that("A6: zero on a non-last speciation node of an augmented tree gives -Inf", {
  # Sequence of +=: obs n=2 (.3), obs n=3 (.2), missing n=4 (.1), missing n=5 (0),
  # two extinctions (skipped), obs n=4 (.1) last; closing sentinel excluded.
  p  <- p8(0.5, -0.1)
  tr <- data.frame(brts = 1:8, n = c(2, 3, 4, 5, 6, 5, 4, 5),
                   t_ext = c(1e11, 1e11, 5, 6, 0, 0, 1e11, 1e11),
                   pd = 0, tip_start = c(0, 0, 3, 4, 3, 4, 0, 0),
                   id = c(0L, 1L, 3L, 4L, 3L, 4L, 2L, -1L),
                   parent_id = c(-1L, -1L, 1L, 1L, 1L, 1L, -1L, -1L))
  expect_identical(lam(p, 5), 0)
  expect_gt(lam(p, 4), 0)
  expect_identical(ev(p, tr), -Inf)
})

test_that("a negative linear predictor clamps to rate 0 and also gives -Inf", {
  # lambda(N) = 0.1 - 0.1 N is negative for every N >= 2; speciation_rate clamps to 0.
  tr <- mk_tree(2); p <- p8(0.1, -0.1)
  expect_identical(ev(p, tr), -Inf)
})

test_that("a strictly negative rate (gaussian link, beta_0 < 0) gives -Inf", {
  # Every rate above is exactly 0: the linear link clamps with max(0, .), so it
  # cannot produce a strictly negative one.  The gaussian link does not clamp:
  # lambda = beta_0 * exp(-(eta_cov - 1)^2 / 2), so beta_0 = -0.4 gives
  # lambda = -0.2426 at every node.  Four scored nodes, so the running product
  # cannot end positive by sign cancellation.  Before the fix this returned
  # +Inf, and the "val <= 0" guard is what this pins: a guard weakened to
  # "val == 0" leaves every other expectation in this file green.
  p  <- c(-0.4, 0, 0, 0, 0.2, 0, 0, 0)
  lf <- emphasis:::eval_logf(p, list(mk_tree(4)), model = c(1L, 0L, 0L),
                             link = 2L, rho = 1)$logf
  expect_identical(lf, -Inf)
})

test_that("a zero rate gives -Inf for several trees evaluated in one call", {
  p <- p8(0.3, -0.1)
  trees <- list(mk_tree(1), mk_tree(2), mk_tree(3))
  lf <- emphasis:::eval_logf(p, trees, model = c(1L, 0L, 0L), link = 0L, rho = 1)$logf
  expect_true(is.finite(lf[1]))                 # lambda(2) = 0.1, the only node
  expect_identical(lf[2], -Inf)                 # lambda(3) = 0 last
  expect_identical(lf[3], -Inf)                 # lambda(3) = 0, then lambda(4) = 0
})

test_that("threshold crossings leave finite results unchanged", {
  # Same products as A4/A5 without the zero: eval_logf must match the hand formula.
  tr <- mk_tree(24); p <- c(26 * 4, -4, 0, 0, 0.1, 0, 0, 0)
  expect_equal(ev(p, tr), hand_loglik(p, tr), tolerance = 1e-10)
  tr <- mk_tree(29); p <- c(31 / 128, -1 / 128, 0, 0, 0.1, 0, 0, 0)
  expect_equal(ev(p, tr), hand_loglik(p, tr), tolerance = 1e-10)
})

# BDI E-step under DD: acceptance-corrected fhat, non-finite draws, draw
# counts (audit ids H10, H11(b), H59).  All sampling in .augment_tree_bdi is
# pure R (stats::rexp / runif), so set.seed reproduces every draw.
#
# Fixed trees (branching times, crown age first):
#   brts_dd20  DDD::dd_sim(c(0.8, 0.3, 25), age = 6), seed 11 -> 20 tips
#   brts_dd9   DDD::dd_sim(c(1.5, 0.4, 1.1/0.12), age = 6, ddmodel = 1),
#              seed 11 -> 9 tips
#   brts_cr20  ape::rphylo(20, 0.5, 0.1), seed 4 -> 20 tips

brts_dd20 <- c(6, 4.27428208320388, 4.1258975168475, 3.72207978339226,
               3.64993282835291, 3.17502914508657, 2.81163149722479,
               2.48700670286076, 2.38274529774947, 1.94421117175661,
               1.86664932095232, 1.81164887493753, 1.38503498765192,
               1.06976980095552, 0.460758774034712, 0.417936203054209,
               0.363178692045526, 0.253829586012209, 0.0517782331486725)

brts_dd9 <- c(6, 4.8484929804551, 4.40182149214092, 3.10816371045776,
              3.07391367934619, 2.835023069071, 1.8388281974863,
              0.504630059261689)

brts_cr20 <- c(5.48610854812521, 3.03448761731259, 3.02735597139097,
               2.95889244175666, 2.34868651175486, 2.14359507126436,
               0.769717103897072, 0.692247501721603, 0.648459687588115,
               0.624791735715604, 0.54690155959252, 0.489625935745245,
               0.309663547120309, 0.282539883338892, 0.250347796440948,
               0.244509063995814, 0.216886588327576, 0.0693130709982483,
               0.0143000508348123)

dd_bin <- c(1L, 0L, 0L)
cr_bin <- c(0L, 0L, 0L)

# DDD ddmodel 1: lambda(N) = l0 - (l0 - mu0) N / K, i.e. beta_N = -(l0 - mu0)/K
dd_pars <- function(l0, m0, K) c(l0, -(l0 - m0) / K, m0, 0)

test_that("BDI dd fhat matches DDD::dd_loglik up to a theta-independent constant (H10)", {
  skip_if_not_installed("DDD")
  grid <- expand.grid(K = c(20, 50, 1e4),
                      i = 1:3, stringsAsFactors = FALSE)
  lm <- rbind(c(0.8, 0.3), c(0.6, 0.1), c(1.0, 0.8))
  grid$l0 <- lm[grid$i, 1]
  grid$m0 <- lm[grid$i, 2]

  set.seed(101)
  rows <- lapply(seq_len(nrow(grid)), function(j) {
    l0 <- grid$l0[j]; m0 <- grid$m0[j]; K <- grid$K[j]
    ref <- DDD::dd_loglik(pars1 = c(l0, m0, K),
                          pars2 = c(300, 1, 0, 1, 0, 2),   # cond 0, btorph 1, soc 2
                          brts = brts_dd20, missnumspec = 0)
    e <- .augment_tree_bdi(brts_dd20, dd_pars(l0, m0, K), dd_bin,
                           sample_size = 500L, max_missing = 1e4L,
                           link = 0L, rho = 1)
    c(ref = ref, fhat = e$fhat, acc = e$acc, n_valid = e$n_valid,
      n_attempts = e$n_attempts, n_rejected = e$n_rejected,
      n_rejected_max_missing = e$n_rejected_max_missing)
  })
  res <- as.data.frame(do.call(rbind, rows))

  expect_true(all(is.finite(res$ref)))
  expect_true(all(is.finite(res$fhat)))
  expect_true(all(res$n_valid == 500L))
  expect_true(all(res$n_valid + res$n_rejected + res$n_rejected_max_missing ==
                    res$n_attempts))
  expect_true(all(res$n_rejected_max_missing == 0L))

  # Survivor rejections are real on this tree: acceptance well below 1.
  expect_lt(min(res$acc), 0.75)

  # Corrected estimator: fhat - dd_loglik is a constant (~0) across theta.
  gap <- res$fhat - res$ref
  expect_lt(diff(range(gap)), 0.5)
  expect_lt(max(abs(gap)), 0.3)

  # The estimator without log(acc) is off by -log(acc): > 0.3 nats at the
  # low-acceptance points of the grid.
  gap_uncorrected <- res$fhat - log(res$acc) - res$ref
  expect_gt(max(gap_uncorrected), 0.3)
})

test_that("BDI cr fhat has acc = 1 and equals DDD::bd_loglik(cond = 0, btorph = 1, soc = 2)", {
  skip_if_not_installed("DDD")
  set.seed(7)
  for (lm in list(c(0.8, 0.3), c(0.5, 0.1))) {
    e <- .augment_tree_bdi(brts_dd20, lm, cr_bin, sample_size = 20L,
                           max_missing = 1e4L, link = 0L, rho = 1)
    ref <- DDD::bd_loglik(pars1 = c(lm, 0, 0), pars2 = c(0, 0, 1, 0, 2),
                          brts = brts_dd20, missnumspec = 0)
    expect_equal(e$acc, 1)
    expect_equal(e$n_rejected, 0L)
    expect_equal(e$n_attempts, 20L)
    expect_equal(e$n_nonfinite, 0L)
    expect_lt(stats::sd(e$weights), 1e-8)
    expect_lt(abs(e$fhat - ref), 1e-8)
  }
})

test_that("non-finite log-weights leave the M-step set and stay in the fhat denominator (H11 B1)", {
  pars <- c(1.5, -0.12, 0.4, 0)   # lambda(N) = 0 for N >= 12.5 on a 9-tip tree
  set.seed(2)
  total_nonfinite <- 0L
  for (r in 1:3) {
    e <- .augment_tree_bdi(brts_dd9, pars, dd_bin, sample_size = 200L,
                           max_missing = 30L, link = 0L, rho = 1)
    total_nonfinite <- total_nonfinite + e$n_nonfinite

    expect_equal(e$n_valid, 200L)
    expect_true(all(is.finite(e$weights)))
    expect_true(all(is.finite(e$logf)))
    expect_length(e$trees, e$n_valid - e$n_nonfinite)
    expect_length(e$logf, length(e$trees))
    expect_length(e$logg, length(e$trees))
    expect_length(e$weights, length(e$trees))
    expect_true(is.finite(e$fhat))
    expect_true(is.finite(.ess_from_lw(e$weights)))

    # The self-normalised weights handed to m_cpp by .mcem_bdi contain no NaN.
    w_norm <- exp(e$weights - max(e$weights))
    w_norm <- w_norm / sum(w_norm) * length(w_norm)
    expect_false(anyNA(w_norm))
    expect_true(all(w_norm > 0))

    # fhat = log(sum_w / n_valid) + max_lw + log(acc): the dropped draws
    # count in the denominator, not in the sum.
    max_lw <- max(e$weights)
    expect_equal(e$fhat,
                 log(sum(exp(e$weights - max_lw)) / e$n_valid) + max_lw +
                   log(e$acc))
    expect_equal(e$acc, e$n_valid / (e$n_valid + e$n_rejected))
  }
  expect_gt(total_nonfinite, 0L)
})

test_that("all completed draws non-finite gives an empty M-step set and fhat = -Inf", {
  # lambda(N) = 0.6 - 0.5 N / 15 is zero at N = 18 and negative on the 20
  # observed lineages, so logf is non-finite on every augmented tree.
  set.seed(5)
  e <- .augment_tree_bdi(brts_dd20, dd_pars(0.6, 0.1, 15), dd_bin,
                         sample_size = 20L, max_missing = 1e4L,
                         link = 0L, rho = 1)
  expect_gt(e$n_valid, 0L)
  expect_equal(e$n_nonfinite, e$n_valid)
  expect_length(e$trees, 0L)
  expect_length(e$weights, 0L)
  expect_equal(e$fhat, -Inf)
})

test_that(".augment_tree_bdi returns draw counts and uses a 5 * sample_size budget under CR (H59)", {
  pars <- c(0.5, 0.4)   # high turnover: many missing lineages per draw
  cnt <- c("n_valid", "n_nonfinite", "n_attempts", "n_rejected",
           "n_rejected_max_missing", "acc")

  set.seed(3)
  expect_no_warning(
    e_full <- .augment_tree_bdi(brts_cr20, pars, cr_bin, sample_size = 50L,
                                max_missing = 1e4L, link = 0L, rho = 1))
  expect_true(all(cnt %in% names(e_full)))
  expect_equal(e_full$n_valid, 50L)
  expect_equal(e_full$n_attempts, 50L)
  expect_equal(e_full$n_rejected, 0L)
  expect_equal(e_full$n_rejected_max_missing, 0L)
  expect_equal(e_full$acc, 1)
  expect_length(e_full$trees, 50L)

  # max_missing = 10 rejects some draws; the budget above sample_size
  # still fills the request.
  e_10 <- .augment_tree_bdi(brts_cr20, pars, cr_bin, sample_size = 50L,
                            max_missing = 10L, link = 0L, rho = 1)
  expect_equal(e_10$n_valid, 50L)
  expect_gt(e_10$n_attempts, 50L)
  expect_lte(e_10$n_attempts, 250L)
  expect_gt(e_10$n_rejected_max_missing, 0L)
  expect_equal(e_10$n_valid + e_10$n_rejected + e_10$n_rejected_max_missing,
               e_10$n_attempts)

  # max_missing = 2 exhausts the budget: fewer trees than requested, a
  # warning, and the overflow count accounts for the shortfall.
  expect_warning(
    e_2 <- .augment_tree_bdi(brts_cr20, pars, cr_bin, sample_size = 50L,
                             max_missing = 2L, link = 0L, rho = 1),
    "over max_missing")
  expect_lt(e_2$n_valid, 50L)
  expect_equal(e_2$n_attempts, 250L)
  expect_equal(e_2$n_rejected, 0L)
  expect_equal(e_2$n_rejected_max_missing, 250L - e_2$n_valid)
  expect_equal(e_2$acc, 1)
  expect_length(e_2$trees, e_2$n_valid)
  expect_length(e_2$weights, e_2$n_valid)

  # max_missing overflows are not part of the fhat denominator: the CR
  # estimate (exact IS, constant weights) is the same at every max_missing
  # that returned at least one tree.
  if (e_2$n_valid > 0L) expect_lt(abs(e_2$fhat - e_full$fhat), 1e-8)
  expect_lt(abs(e_10$fhat - e_full$fhat), 1e-8)
})

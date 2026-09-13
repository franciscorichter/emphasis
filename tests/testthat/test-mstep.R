# Tests for the M-step objective (src/M_step.cpp): audit ids H4, H11(c).
#
# The objective is -Q(pars) with Q = sum_i w_i * loglik(pars, tree_i) over
# unnormalised weights; with a conditional it is -Q + sum(w) * conditional(pars),
# so the argmin does not depend on a positive rescaling of the weights.
# Terms with w_i == 0 are skipped and a non-finite loglik on a weighted tree
# makes the objective +Inf rather than NaN.

# ---- fixed inputs -----------------------------------------------------------

# 25-tip constant-rate tree (TreeSim seed 5, lambda = 1, mu = 0.5), crown age 4.088
brts_cr <- c(4.087931, 3.824527, 3.380355, 1.879518, 1.837783, 1.474373,
             1.182076, 1.172046, 1.067297, 0.847818, 0.516806, 0.493594,
             0.415446, 0.379772, 0.37824, 0.353089, 0.318518, 0.223403,
             0.222619, 0.181101, 0.172338, 0.161752, 0.076436, 0.058077)

# 8-tip diversity-dependent tree (DDD::dd_sim seed 11, pars (1.5, 0.4, K = 1.1/0.12), age 6)
brts_dd <- c(6, 4.848493, 4.401821, 3.108164, 3.073914, 2.835023, 1.838828, 0.50463)

# exact log P(both crown lineages survive to the present) under constant rates
log_psurv <- function(lam, mu, T) {
  if (abs(lam - mu) < 1e-8) mu <- lam - 1e-8
  r  <- lam - mu
  p0 <- mu * (1 - exp(-r * T)) / (lam - mu * exp(-r * T))
  2 * log(max(1 - p0, 1e-300))
}
cond_cr <- function(pars8) log_psurv(pars8[1], pars8[5], T = brts_cr[1])

# m_cpp on a fixed set of augmented trees with given weights
run_mstep <- function(trees, w, init, lb, ub, model_bin, cond = NULL, xtol_rel = 1e-6) {
  es <- list(trees = trees, weights = w, rejected = 0L, rejected_overruns = 0L,
             rejected_lambda = 0L, rejected_zero_weights = 0L, time = 0, fhat = 0)
  r <- m_cpp(e_step = es, init_pars = init, plugin = "rpd1",
             lower_bound = lb, upper_bound = ub, xtol_rel = xtol_rel,
             num_threads = 1L, model = model_bin, link = 0L, rho = 1,
             rconditional = cond)
  as.numeric(r$estimates)
}

# weights normalised to mean 1 from log-weights, as .mcem_bdi hands them to the M-step
mean_one <- function(lw) {
  w <- exp(lw - max(lw))
  w / sum(w) * length(w)
}

# ---- constant-rate E-step shared by the H4 tests ----------------------------

cr_bin <- c(0L, 0L, 0L)
lb_cr  <- c(1e-3, 0, 0, 0, 1e-3, 0, 0, 0)
ub_cr  <- c(5,    0, 0, 0, 5,    0, 0, 0)
theta0 <- c(1, 0, 0, 0, 0.5, 0, 0, 0)
N_cr   <- 200L

set.seed(1)
e_cr <- .augment_tree_bdi(brts_cr, theta0, model_bin = cr_bin, sample_size = N_cr,
                          max_missing = 1e4, link = 0L, rho = 1)
w_cr <- mean_one(e_cr$weights)

mstep_cr <- function(w, cond = NULL) {
  run_mstep(e_cr$trees, w, theta0, lb_cr, ub_cr, cr_bin, cond)[c(1, 5)]
}

test_that("M-step estimates are invariant to a positive rescaling of the weights", {
  expect_equal(sum(w_cr), N_cr)
  for (cond in list(NULL, cond_cr)) {
    ref <- mstep_cr(w_cr, cond)
    expect_true(all(is.finite(ref)))
    expect_equal(mstep_cr(w_cr / 256, cond), ref, tolerance = 1e-6)   # exact scaling
    expect_equal(mstep_cr(w_cr / N_cr, cond), ref, tolerance = 1e-6)  # sum(w) = 1
    expect_equal(mstep_cr(w_cr * 10, cond), ref, tolerance = 1e-6)    # sum(w) = 10 N
  }
})

test_that("the conditional changes the argmin by more than the dilution 1/sum(w) allowed", {
  est_uncond <- mstep_cr(w_cr)
  est_cond   <- mstep_cr(w_cr, cond_cr)
  # H4: with the penalty diluted by 1/N the two differed by about 1e-4
  expect_gt(max(abs(est_cond - est_uncond)), 0.01)
})

test_that("M-step argmin equals the optim maximiser of the self-normalised objective", {
  Qn <- function(p) {
    lf <- eval_logf(c(p[1], 0, 0, 0, p[2], 0, 0, 0), e_cr$trees,
                    model = cr_bin, link = 0L, rho = 1)$logf
    sum(w_cr * lf) / sum(w_cr)
  }
  ref <- function(pen) {
    optim(c(0.7, 0.4), function(p) {
      if (any(p < 1e-3) || any(p > 5)) return(1e10)
      -(Qn(p) - pen * log_psurv(p[1], p[2], T = brts_cr[1]))
    }, control = list(reltol = 1e-12, maxit = 5000))$par
  }
  expect_equal(mstep_cr(w_cr), ref(0), tolerance = 1e-3)             # argmax Q_norm
  expect_equal(mstep_cr(w_cr, cond_cr), ref(1), tolerance = 1e-3)    # argmax Q_norm - log P
})

# ---- diversity-dependent E-step with zero-density trees (H11) ---------------

test_that("an E-step containing a -Inf tree still moves the estimate", {
  dd_bin <- c(1L, 0L, 0L)
  pars   <- c(1.5, -0.12, 0.4, 0)            # lambda(N) = 1.5 - 0.12 N, zero at N = 12.5
  lb_dd  <- c(0.01, -1, 0.001, 0)
  ub_dd  <- c(5, 0, 2, 0)
  pars8  <- c(pars[1], pars[2], 0, 0, pars[3], pars[4], 0, 0)
  lb8    <- .expand_pars(lb_dd, dd_bin)
  ub8    <- .expand_pars(ub_dd, dd_bin)

  set.seed(2)
  e_dd <- .augment_tree_bdi(brts_dd, pars, model_bin = dd_bin, sample_size = 200L,
                            max_missing = 30L, link = 0L, rho = 1)
  # The sampler may or may not return zero-density draws; the M-step set here is
  # built from the finite-logf draws so that the test does not depend on it.
  keep <- is.finite(e_dd$logf)
  expect_gt(sum(keep), 20)
  trees_fin <- e_dd$trees[keep]
  w_fin <- mean_one(e_dd$weights[keep])

  # Hand-built tree in the augmented-tree convention (forward times, crown at 0,
  # n = lineages alive during the segment ending at brts): 7 observed
  # speciations (n 2..8), 6 missing births (n 8..13, lambda(13) = 0), their 6
  # extinctions (n 14..9), one more observed speciation at n = 8, closing
  # sentinel at 9. The zero-rate node is not the last speciation node.
  mb <- 7 + (1:6) / 10
  me <- 8 + (1:6) / 10
  tree_bad <- data.frame(brts = c(1:7, mb, me, 8.8, 9),
                         n = c(2:8, 8:13, 14:9, 8, 9),
                         t_ext = c(rep(1e11, 7), me, rep(0, 6), 1e11, 1e11),
                         pd = 0,
                         tip_start = c(rep(0, 7), mb, mb, 0, 0),
                         id = c(0:6, 8:13, 8:13, 7L, -1L),
                         parent_id = c(rep(-1L, 7), rep(0L, 12), -1L, -1L))
  logf_bad <- eval_logf(pars8, list(tree_bad), model = dd_bin, link = 0L, rho = 1)$logf
  expect_identical(logf_bad, -Inf)

  est_fin <- run_mstep(trees_fin, w_fin, pars8, lb8, ub8, dd_bin, xtol_rel = 1e-4)
  expect_true(all(is.finite(est_fin)))
  expect_gt(max(abs(est_fin - pars8)), 1e-3)

  # w == 0 on the -Inf tree: the term is skipped. The extra element changes
  # the order of the parallel reduction, so SBPLX can stop at a different
  # point of a flat direction; the invariant is the objective, not the
  # coordinates.
  est_zero <- run_mstep(c(trees_fin, list(tree_bad)), c(w_fin, 0), pars8,
                        lb8, ub8, dd_bin, xtol_rel = 1e-4)
  expect_true(all(is.finite(est_zero)))
  Q_fin <- function(th) sum(w_fin * eval_logf(th, trees_fin,
                                              model = dd_bin, link = 0L, rho = 1)$logf)
  expect_equal(Q_fin(est_zero), Q_fin(est_fin), tolerance = 1e-6)

  # H11.R part E: all sampler draws with weight 0 on the non-finite ones;
  # the weights differ from w_fin by the constant factor 200 / sum(keep)
  if (any(!keep)) {
    lw_fix <- ifelse(keep, e_dd$weights, -Inf)
    est_fix <- run_mstep(e_dd$trees, mean_one(lw_fix), pars8, lb8, ub8, dd_bin,
                         xtol_rel = 1e-4)
    expect_equal(est_fix, est_fin, tolerance = 1e-6)
  }

  # positive weight on the -Inf tree: the objective is +Inf at the init and the
  # optimizer ends where every weighted tree has finite density (lambda(13) > 0)
  est_pos <- run_mstep(c(trees_fin, list(tree_bad)), c(w_fin, 1), pars8,
                       lb8, ub8, dd_bin, xtol_rel = 1e-4)
  expect_true(all(is.finite(est_pos)))
  expect_gt(est_pos[1] + 13 * est_pos[2], 0)
  logf_at_est <- eval_logf(est_pos, c(trees_fin, list(tree_bad)),
                           model = dd_bin, link = 0L, rho = 1)$logf
  expect_true(all(is.finite(logf_at_est)))
})

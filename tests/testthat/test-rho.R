# Pins the sampling-fraction rho for H2, H29 and H79.
#
# H2 was that the BDI proposal drew only extinct and observed lineages, never an
# unsampled extant one, so at rho < 1 it sampled the rho = 1 conditioned process
# while eval_logf scored it under rho: logf shifted by n_tips*log(rho) and the
# estimator stayed on the rho = 1 surface.  The wave-2 answer was to refuse
# rho < 1 in .bdi_supported() and fall back to thinning.  The proposal now
# samples the rho-conditioned process itself -- p(t) is the probability of
# leaving a *sampled* descendant and a surviving missing lineage is emitted as
# an unsampled extant tip -- so the gate no longer reads rho, and what this file
# pins is the behaviour that replaced the fallback:
#   - .bdi_supported() is the same at rho = 1 and below it
#   - the BDI draws carry unsampled extant lineages, as the thinning draws do
#   - a cr fit at rho = 0.5 lands on the rho = 0.5 MLE, which is the original
#     H2 gate and is now met by the BDI sampler rather than by routing past it
#   - the fallback to thinning is announced unconditionally, not only under
#     control$verbose, for what is still out of scope
#   - emphasis_pipeline() inherits a top-level control$rho into every stage,
#     a nested control$<stage>$rho overrides it, and the rho used is recorded
#     in the stage log, in each stage's fit and in the pipeline result (H29)
#   - rho is validated in (0, 1] in R, where C++ (model.hpp:97) silently
#     substitutes 1 for anything outside that range (H79)
#
# Branching times are fixed (the 16-tip tree of dev/audit/checks/H2.R:
# ape::rphylo(16, 0.5, 0.3) with set.seed(11)).  num_threads = 1 throughout.
# set.seed() does not reach the augmentation RNG (a C++ mt19937 seeded from
# the clock and the thread id), so the fit-level expectations are margin-based:
# the rho = 1 and rho = 0.5 MLEs of this tree are 0.47 and 0.69, and the
# threshold sits between them.

brts16 <- c(4.225039, 3.654389, 3.210154, 3.001428, 2.078383, 1.516088,
            1.515195, 1.319225, 1.219620, 1.052087, 0.894478, 0.780463,
            0.622042, 0.431325, 0.062668)
n16 <- length(brts16) + 1L

cr_bin <- c(0L, 0L, 0L)
dd_bin <- c(1L, 0L, 0L)
nd_bin <- c(1L, 0L, 1L)
d_bin  <- c(0L, 0L, 1L)
ex <- function(p) emphasis:::.expand_pars(p, cr_bin)

# Nee et al. crown-conditioned log-density with Bernoulli sampling rho, up to a
# theta-independent constant (dev/audit/checks/H2.R).  Used only to locate the
# MLE the fit must land on, which is a property of the tree, not of the code.
nee_rho <- function(lam, mu, rho, brts) {
  r   <- lam - mu
  den <- rho * lam + (lam * (1 - rho) - mu) * exp(-r * brts)
  logp1 <- log(rho) + 2 * log(r) - r * brts - 2 * log(den)
  2 * logp1[1L] + sum(log(lam) + logp1[-1L])
}
nee_mle <- function(rho, brts) {
  stats::optim(c(0.4, 0.2),
               function(p) if (p[1L] <= p[2L] || p[2L] < 0) 1e10 else
                 -nee_rho(p[1L], p[2L], rho, brts),
               method = "L-BFGS-B", lower = c(1e-3, 0), upper = c(5, 5))$par
}

# Unsampled extant lineages carry the extinction sentinel 5e10 (model.hpp:383).
n_unsampled <- function(tr) sum(tr$t_ext == 5e10)


test_that(".bdi_supported does not depend on rho (H2)", {
  # The model/link decisions the validation corpus was fitted under, unchanged.
  expect_true(emphasis:::.bdi_supported(cr_bin, 0L))
  expect_true(emphasis:::.bdi_supported(cr_bin, 1L))
  expect_true(emphasis:::.bdi_supported(dd_bin, 0L))
  expect_true(emphasis:::.bdi_supported(dd_bin, 1L))
  expect_false(emphasis:::.bdi_supported(nd_bin, 0L))   # D covariate
  expect_false(emphasis:::.bdi_supported(d_bin,  0L))
  # cr on the gaussian link is in scope (the covariate part of eta is zero
  # there, so the rate is the constant beta_0*exp(-1/2)); dd on it is not.
  expect_true(emphasis:::.bdi_supported(cr_bin, 2L))
  expect_false(emphasis:::.bdi_supported(dd_bin, 2L))

  # The verdict is the same at every rho in range: the sampler conditions on
  # the rho-sampled survival probability rather than routing around it.
  for (rho in c(1, 0.999, 0.5, 0.1, 1e-6))
    for (mb in list(cr_bin, dd_bin, nd_bin, d_bin))
      for (lk in 0:2)
        expect_identical(emphasis:::.bdi_supported(mb, lk, rho = rho),
                         emphasis:::.bdi_supported(mb, lk))

  # A rho that is not a usable number is not treated as complete sampling.
  expect_false(emphasis:::.bdi_supported(cr_bin, 0L, rho = NA_real_))
  expect_false(emphasis:::.bdi_supported(cr_bin, 0L, rho = NULL))
  expect_false(emphasis:::.bdi_supported(cr_bin, 0L, rho = 0))
  expect_false(emphasis:::.bdi_supported(cr_bin, 0L, rho = 1.5))
})


test_that(".check_rho accepts (0, 1] and rejects everything else (H79)", {
  expect_identical(emphasis:::.check_rho(1), 1)
  expect_identical(emphasis:::.check_rho(0.5), 0.5)
  expect_identical(emphasis:::.check_rho(1L), 1)
  expect_identical(emphasis:::.check_rho(1e-8), 1e-8)

  # model.hpp:97 replaces each of these by 1 and reports nothing.
  for (bad in list(0, -1, 80, 1 + 1e-9, NA_real_, NaN, Inf, -Inf))
    expect_error(emphasis:::.check_rho(bad), "must be a single number in \\(0, 1\\]")
  expect_error(emphasis:::.check_rho(NULL), "got NULL")
  expect_error(emphasis:::.check_rho("0.5"), "must be a single number in \\(0, 1\\]")
  expect_error(emphasis:::.check_rho(c(0.5, 0.5)), "must be a single number in \\(0, 1\\]")
  # The message names the argument the caller used.
  expect_error(emphasis:::.check_rho(3, arg = "rho"), "^rho must be")
})


test_that("estimate_rates validates rho for every method (H79)", {
  ctrl <- function(rho, ...) c(list(rho = rho, lower_bound = c(1e-3, 0),
                                    upper_bound = c(3, 3), num_threads = 1L), list(...))
  for (m in c("mcem", "cem", "gam"))
    expect_error(estimate_rates(brts16, method = m, model = "cr",
                                init_pars = c(0.4, 0.2), control = ctrl(80)),
                 "must be a single number in \\(0, 1\\]")
  expect_error(estimate_rates(brts16, method = "mcem", model = "cr",
                              init_pars = c(0.4, 0.2), control = ctrl(-1)),
               "must be a single number in \\(0, 1\\]")
  expect_error(estimate_rates(brts16, method = "mcem", model = "cr",
                              init_pars = c(0.4, 0.2), control = ctrl(0)),
               "must be a single number in \\(0, 1\\]")
  # The check runs after the control merge, so it sees the user's value and
  # not the default; and it runs before the bound check, which is why an
  # otherwise-invalid call still reports rho.
  expect_error(estimate_rates(brts16, method = "mcem", model = "cr",
                              control = list(rho = 2)),
               "must be a single number in \\(0, 1\\]")
})


test_that("both proposals draw unsampled extant lineages at rho < 1 (H2)", {
  # The defect H2 recorded was that this count was 0 for the BDI draws and up
  # to 17 for the thinning draws on the same tree.  Both are now non-zero, and
  # the rho term of eval_logf carries the unsampled count on both.
  set.seed(101)
  bdi <- emphasis:::.augment_tree_bdi(brts16, pars = c(0.5, 0.3),
                                      model_bin = cr_bin, sample_size = 30L,
                                      link = 0L, rho = 0.5)
  expect_gt(length(bdi$trees), 0L)
  nb <- vapply(bdi$trees, n_unsampled, 1L)
  expect_gt(max(nb), 0L)

  thin <- emphasis:::.augment_tree_internal(brts16, pars = c(0.5, 0.3),
                                            model_bin = cr_bin, sample_size = 30L,
                                            maxN = 3000L, link = 0L, rho = 0.5,
                                            num_threads = 1L)
  nu <- vapply(thin$trees, n_unsampled, 1L)
  expect_gt(max(nu), 0L)

  p8 <- ex(c(0.5, 0.3))
  l1 <- emphasis:::eval_logf(p8, bdi$trees, model = cr_bin, link = 0L, rho = 1)$logf
  l5 <- emphasis:::eval_logf(p8, bdi$trees, model = cr_bin, link = 0L, rho = 0.5)$logf
  expect_equal(l5 - l1, n16 * log(0.5) + nb * log(0.5), tolerance = 1e-12)

  # On the thinning draws the same term carries the unsampled count as well.
  lt1 <- emphasis:::eval_logf(p8, thin$trees, model = cr_bin, link = 0L, rho = 1)$logf
  lt5 <- emphasis:::eval_logf(p8, thin$trees, model = cr_bin, link = 0L, rho = 0.5)$logf
  expect_equal(lt5 - lt1, n16 * log(0.5) + nu * log(0.5), tolerance = 1e-12)

  # At rho = 1 neither sampler writes the sentinel.
  set.seed(101)
  b1 <- emphasis:::.augment_tree_bdi(brts16, pars = c(0.5, 0.3),
                                     model_bin = cr_bin, sample_size = 30L,
                                     link = 0L, rho = 1)
  expect_equal(max(vapply(b1$trees, n_unsampled, 1L)), 0L)
})


test_that("estimate_rates(rho < 1) keeps the BDI sampler and says nothing (H2)", {
  ctrl <- list(rho = 0.5, sampling = "bdi", sample_size = 20L, max_iter = 2L,
               num_threads = 1L, lower_bound = c(1e-3, 0), upper_bound = c(3, 3))
  expect_no_message(
    fit <- estimate_rates(brts16, method = "mcem", model = "cr",
                          init_pars = c(0.4, 0.2), control = ctrl))
  # The trace is the BDI driver's (n_valid and rejected_max_missing); the
  # thinning driver writes maxN and four rejection channels instead.
  expect_true(all(c("n_valid", "rejected_max_missing") %in% names(fit$details$mcem)))
  expect_false("maxN" %in% names(fit$details$mcem))
  expect_equal(fit$rho, 0.5)

  # rho = 1 takes the same route.
  expect_no_message(
    fit1 <- estimate_rates(brts16, method = "mcem", model = "cr",
                           init_pars = c(0.4, 0.2),
                           control = utils::modifyList(ctrl, list(rho = 1))))
  expect_true(all(c("n_valid", "rejected_max_missing") %in% names(fit1$details$mcem)))
  expect_equal(fit1$rho, 1)
})


test_that("the fallback still fires, unconditionally, for what is out of scope", {
  # A D-dependent model at rho < 1 still goes to thinning, and the message
  # names the reason rather than only the fact.
  # model "nd" is c(1, 0, 1): beta_0, beta_N, beta_D, gamma_0, gamma_N, gamma_D.
  ctrl <- list(rho = 0.5, sampling = "bdi", sample_size = 20L, max_iter = 2L,
               maxN = 2000L, num_threads = 1L,
               lower_bound = c(1e-3, -1, -1, 0, 0, 0),
               upper_bound = c(3, 1, 1, 3, 0, 0))
  init <- c(0.4, -0.01, 0.01, 0.2, 0, 0)
  expect_message(
    fit <- estimate_rates(brts16, method = "mcem", model = "nd",
                          init_pars = init, control = ctrl),
    "D-dependent")
  # Not gated on verbose, and it names the sampler that ran.
  expect_message(
    estimate_rates(brts16, method = "mcem", model = "nd", init_pars = init,
                   control = utils::modifyList(ctrl, list(verbose = FALSE))),
    "dynamic_fresh")
  expect_true(all(c("maxN", "rejected_errors") %in% names(fit$details$mcem)))
  expect_equal(fit$rho, 0.5)
})


test_that("a cr fit at rho = 0.5 lands on the rho = 0.5 MLE (H2 gate)", {
  m1  <- nee_mle(1.0, brts16)
  m05 <- nee_mle(0.5, brts16)
  # The two surfaces this tree admits: lambda 0.4697 at rho = 1, 0.6872 at 0.5.
  expect_equal(m1[1L],  0.4697, tolerance = 1e-3)
  expect_equal(m05[1L], 0.6872, tolerance = 1e-3)

  fit <- suppressMessages(estimate_rates(
    brts16, method = "mcem", model = "cr", init_pars = c(0.4, 0.2),
    control = list(rho = 0.5, sampling = "bdi", sample_size = 200L,
                   max_iter = 25L, max_time = 600, num_threads = 1L,
                   lower_bound = c(1e-3, 0), upper_bound = c(3, 3))))
  lam <- unname(fit$pars[1L])
  # Measured before the fix (dev/audit/checks/H2.R): 0.4727, 0.4785, 0.4713 —
  # the rho = 1 MLE.  The fallback to thinning gave 0.6782 to 0.7071.  The BDI
  # sampler now meets the same gate itself, and exactly: its log-weights on this
  # tree are constant to 1e-12 at rho = 0.5, so the E-step contributes no Monte
  # Carlo error at all and what is left of the gap is the M-step's.
  expect_gt(lam, 0.60)
  expect_lt(abs(lam - m05[1L]), abs(lam - m1[1L]))
  expect_true("n_valid" %in% names(fit$details$mcem))   # the BDI driver ran
  expect_lt(stats::sd(fit$details$final_IS$lw), 1e-9)
})


test_that("emphasis_pipeline inherits top-level rho into its stages (H29)", {
  base <- list(lower_bound = c(1e-3, 0), upper_bound = c(3, 3), num_threads = 1L,
               mcem = list(sample_size = 20L, max_iter = 2L))

  res <- suppressMessages(emphasis_pipeline(
    brts16, model = "cr", link = "linear", stages = "mcem",
    control = utils::modifyList(base, list(rho = 0.5)), verbose = FALSE))
  expect_equal(res$fits$mcem$rho, 0.5)
  expect_equal(res$log$rho, 0.5)
  expect_equal(res$rho, 0.5)
  # rho reached the back end, not just the record.  The stage no longer changes
  # sampler on rho, so the evidence is the likelihood itself: the BDI E-step is
  # exact under constant rates, so fhat at the returned theta is the closed-form
  # rho-sampled likelihood there.  Had rho stopped at the record, fhat would be
  # the rho = 1 value, which differs by n_tips*log(rho) = 16*log(0.5) = -11.09.
  fit05 <- res$fits$mcem
  expect_lt(stats::sd(fit05$details$final_IS$lw), 1e-9)
  expect_lt(abs(fit05$loglik -
                nee_rho(fit05$pars[1L], fit05$pars[2L], 0.5, brts16)), 1e-8)
  expect_gt(abs(fit05$loglik -
                nee_rho(fit05$pars[1L], fit05$pars[2L], 1.0, brts16)), 1)

  # A nested per-stage value wins over the inherited one.
  ctrl_nested <- utils::modifyList(base, list(rho = 0.5))
  ctrl_nested$mcem$rho <- 1.0
  res2 <- emphasis_pipeline(brts16, model = "cr", link = "linear", stages = "mcem",
                            control = ctrl_nested, verbose = FALSE)
  expect_equal(res2$fits$mcem$rho, 1)
  expect_equal(res2$log$rho, 1)
  expect_true("n_valid" %in% names(res2$fits$mcem$details$mcem))

  # Default is complete sampling, recorded as such.
  res3 <- emphasis_pipeline(brts16, model = "cr", link = "linear", stages = "mcem",
                            control = base, verbose = FALSE)
  expect_equal(res3$rho, 1)
  expect_equal(res3$log$rho, 1)

  # And the pipeline validates rho before any stage runs.
  expect_error(emphasis_pipeline(brts16, model = "cr", stages = "mcem",
                                 control = utils::modifyList(base, list(rho = 2)),
                                 verbose = FALSE),
               "must be a single number in \\(0, 1\\]")
})


test_that("rho is reported only when it is below 1", {
  fit <- estimate_rates(brts16, method = "mcem", model = "cr",
                        init_pars = c(0.4, 0.2),
                        control = list(sample_size = 20L, max_iter = 2L,
                                       num_threads = 1L,
                                       lower_bound = c(1e-3, 0),
                                       upper_bound = c(3, 3)))
  expect_false(any(grepl("rho", capture.output(print(fit)))))
  fit$rho <- 0.5
  expect_match(paste(capture.output(print(fit)), collapse = "\n"),
               "rho:\\s+0\\.5 \\(incomplete sampling\\)")
})

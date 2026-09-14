# Pins the sampling-fraction rho after the wave-2 fix for H2, H29 and H79:
#   - .bdi_supported() takes rho and is FALSE below 1: the BDI proposal draws
#     only extinct and observed lineages, never an unsampled extant one, so at
#     rho < 1 it samples the rho = 1 conditioned process while eval_logf scores
#     it under rho, which shifts logf by n_tips*log(rho) and leaves the
#     estimator on the rho = 1 likelihood surface (H2)
#   - the fallback to thinning is announced unconditionally, not only under
#     control$verbose
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


test_that(".bdi_supported is unchanged at rho = 1 and FALSE below it (H2)", {
  # The model/link decisions the validation corpus was fitted under.
  expect_true(emphasis:::.bdi_supported(cr_bin, 0L))
  expect_true(emphasis:::.bdi_supported(cr_bin, 1L))
  expect_true(emphasis:::.bdi_supported(dd_bin, 0L))
  expect_true(emphasis:::.bdi_supported(dd_bin, 1L))
  expect_false(emphasis:::.bdi_supported(cr_bin, 2L))   # gaussian link
  expect_false(emphasis:::.bdi_supported(nd_bin, 0L))   # D covariate
  expect_false(emphasis:::.bdi_supported(d_bin,  0L))
  # The new argument defaults to complete sampling, so the two-argument calls
  # above and the explicit rho = 1 calls agree.
  for (mb in list(cr_bin, dd_bin, nd_bin, d_bin))
    for (lk in 0:2)
      expect_identical(emphasis:::.bdi_supported(mb, lk),
                       emphasis:::.bdi_supported(mb, lk, rho = 1))

  # Below 1 every model/link goes to thinning.
  for (rho in c(0.999, 0.5, 0.1, 1e-6))
    for (mb in list(cr_bin, dd_bin))
      for (lk in 0:1)
        expect_false(emphasis:::.bdi_supported(mb, lk, rho = rho))

  # A rho that is not a usable number is not treated as complete sampling.
  expect_false(emphasis:::.bdi_supported(cr_bin, 0L, rho = NA_real_))
  expect_false(emphasis:::.bdi_supported(cr_bin, 0L, rho = NULL))
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


test_that("the BDI proposal draws no unsampled extant lineage; thinning does (H2)", {
  bdi <- emphasis:::.augment_tree_bdi(brts16, pars = c(0.5, 0.3),
                                      model_bin = cr_bin, sample_size = 30L,
                                      link = 0L, rho = 0.5)
  expect_gt(length(bdi$trees), 0L)
  expect_equal(max(vapply(bdi$trees, n_unsampled, 1L)), 0L)

  thin <- emphasis:::.augment_tree_internal(brts16, pars = c(0.5, 0.3),
                                            model_bin = cr_bin, sample_size = 30L,
                                            maxN = 3000L, link = 0L, rho = 0.5,
                                            num_threads = 1L)
  nu <- vapply(thin$trees, n_unsampled, 1L)
  expect_gt(max(nu), 0L)

  # This is why scoring BDI draws at rho < 1 only shifts logf by a constant:
  # with no unsampled lineage the rho term of eval_logf is n_obs * log(rho).
  p8 <- ex(c(0.5, 0.3))
  l1 <- emphasis:::eval_logf(p8, bdi$trees, model = cr_bin, link = 0L, rho = 1)$logf
  l5 <- emphasis:::eval_logf(p8, bdi$trees, model = cr_bin, link = 0L, rho = 0.5)$logf
  expect_equal(l5 - l1, rep(n16 * log(0.5), length(l1)), tolerance = 1e-12)

  # On the thinning draws the same term carries the unsampled count as well.
  lt1 <- emphasis:::eval_logf(p8, thin$trees, model = cr_bin, link = 0L, rho = 1)$logf
  lt5 <- emphasis:::eval_logf(p8, thin$trees, model = cr_bin, link = 0L, rho = 0.5)$logf
  expect_equal(lt5 - lt1, n16 * log(0.5) + nu * log(0.5), tolerance = 1e-12)
})


test_that("estimate_rates(rho < 1) says the sampler changed and runs thinning (H2)", {
  ctrl <- list(rho = 0.5, sampling = "bdi", sample_size = 20L, max_iter = 2L,
               num_threads = 1L, lower_bound = c(1e-3, 0), upper_bound = c(3, 3))
  expect_message(
    fit <- estimate_rates(brts16, method = "mcem", model = "cr",
                          init_pars = c(0.4, 0.2), control = ctrl),
    "rho = 0.5")
  # The trace is the thinning driver's (maxN and the four rejection channels);
  # the BDI driver writes n_valid / rejected_max_missing instead.
  expect_true(all(c("maxN", "rejected_errors") %in% names(fit$details$mcem)))
  expect_false("n_valid" %in% names(fit$details$mcem))
  expect_equal(fit$rho, 0.5)

  # The message is not gated on verbose, and it names the sampler that ran.
  expect_message(
    estimate_rates(brts16, method = "mcem", model = "cr", init_pars = c(0.4, 0.2),
                   control = utils::modifyList(ctrl, list(verbose = FALSE))),
    "dynamic_fresh")

  # At rho = 1 the same call keeps BDI and says nothing.
  expect_no_message(
    fit1 <- estimate_rates(brts16, method = "mcem", model = "cr",
                           init_pars = c(0.4, 0.2),
                           control = utils::modifyList(ctrl, list(rho = 1))))
  expect_true(all(c("n_valid", "rejected_max_missing") %in% names(fit1$details$mcem)))
  expect_equal(fit1$rho, 1)
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
  # the rho = 1 MLE.  After it, thinning gives 0.6782 to 0.7071.  The threshold
  # is between the two clusters.
  expect_gt(lam, 0.60)
  expect_lt(abs(lam - m05[1L]), abs(lam - m1[1L]))
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
  # rho reached the back end, not just the record: the stage ran thinning,
  # which is a decision .run_mcem takes from rho alone on a cr/linear tree.
  expect_true("maxN" %in% names(res$fits$mcem$details$mcem))

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

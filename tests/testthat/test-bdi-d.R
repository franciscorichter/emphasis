# The surrogate sampler on a D model: the mean field reads D = 0 (D is
# centred) and the attachment tilt carries beta_D.  The proposal is not the
# model, so the importance weights must carry the difference: E_q[f/q] has to
# be the marginal likelihood, which on the 3-tip tree of test-parent-support
# is a two-dimensional quadrature (.L1_3 there, over one missing lineage).

.src_ps <- function() {
  env <- new.env()
  lines <- readLines(testthat::test_path("test-parent-support.R"))
  eval(parse(text = lines[seq_len(grep("^test_that", lines)[1L] - 1L)]), envir = env)
  env
}

test_that("the surrogate sampler admits D models on the linear and exponential links", {
  expect_true(emphasis:::.bdi_supported(c(1L, 0L, 1L, 0L), 0L, 1))
  expect_true(emphasis:::.bdi_supported(c(0L, 0L, 1L, 0L), 1L, 1))
  expect_true(emphasis:::.bdi_supported(c(1L, 0L, 1L, 1L), 1L, 1))
  expect_false(emphasis:::.bdi_supported(c(1L, 0L, 1L, 0L), 2L, 1))
  expect_false(emphasis:::.bdi_supported(c(1L, 1L, 0L, 0L), 0L, 1))
  expect_null(emphasis:::.bdi_unsupported_reason(c(1L, 0L, 1L, 0L), 0L, 1))
})

test_that("a D model runs on the thinning proposal by default and on the surrogate on request", {
  skip_on_cran()
  set.seed(8)
  phy <- ape::rphylo(20L, 0.5, 0.05)
  ctrl <- list(sample_size = 20L, max_iter = 1L, num_threads = 1L, verbose = FALSE,
               lower_bound = c(0.05, -0.05, -0.2, 0, -0.05, -0.2), upper_bound = c(2, 0.02, 0.2, 1, 0.02, 0.2),
               maxN = 10000L, max_missing = 1e4)
  f1 <- suppressWarnings(suppressMessages(estimate_rates(phy, method = "mcem", model = "nd",
                                                          init_pars = c(0.5, -0.01, 0, 0.05, 0, 0), control = ctrl)))
  expect_equal(f1$sampling, "dynamic_fresh")
  f2 <- suppressWarnings(suppressMessages(estimate_rates(phy, method = "mcem", model = "nd",
                                                          init_pars = c(0.5, -0.01, 0, 0.05, 0, 0),
                                                          control = c(ctrl, list(sampling = "bdi")))))
  expect_equal(f2$sampling, "bdi")
  # a model without D keeps the surrogate default
  f3 <- suppressWarnings(suppressMessages(estimate_rates(phy, method = "mcem", model = "dd",
                                                          init_pars = c(0.5, -0.01, 0.05, 0),
                                                          control = list(sample_size = 20L, max_iter = 1L, num_threads = 1L,
                                                                         lower_bound = c(0.05, -0.05, 0, -0.05), upper_bound = c(2, 0.02, 1, 0.02)))))
  expect_equal(f3$sampling, "bdi")
})

test_that("E_q[f/q] under the surrogate sampler is the marginal likelihood of a D model", {
  skip_on_cran()
  e <- .src_ps()
  cf <- e$.cfg3
  brts <- c(cf$TT, cf$TT - cf$t1)
  l0 <- exp(emphasis:::eval_logf(cf$p, list(e$.aug0_3(cf$TT, cf$t1)), model = cf$mb, link = cf$link, rho = 1)$logf)
  l1 <- e$.L1_3(cf$TT, cf$t1, cf$p, cf$mb, cf$link, 40L)
  # the tree with its topology: the observed split at t1 is of crown lineage
  # -2, whose tip start before it is 0 (one entry per event, plus a sentinel)
  tr <- brts; attr(tr, "parent_tip_start") <- c(0, -1); attr(tr, "parent_id") <- c(-2L, -1L)
  set.seed(7)
  a <- emphasis:::.augment_tree_bdi(tr, cf$p, model_bin = cf$mb, sample_size = 4000L, link = cf$link, rho = 1)
  skip_if(length(a$trees) < 100L, "too few draws")
  lw <- emphasis:::eval_logf(cf$p, a$trees, model = cf$mb, link = cf$link, rho = 1)$logf - a$logg
  nm <- vapply(a$trees, function(d) sum(e$.is_mis(d$t_ext)), 0)
  w  <- exp(lw)
  # draws with at most one missing lineage against L0 + L1, as the thinning gate does
  acc <- if (is.null(a$acc) || !is.finite(a$acc)) 1 else a$acc
  wle <- ifelse(nm <= 1, w, 0)
  est <- mean(wle) * acc; se <- stats::sd(wle) * acc / sqrt(length(wle))
  expect_lt(abs(est - (l0 + l1)) / se, 5)
  expect_lt(abs(est - (l0 + l1)) / (l0 + l1), 0.10)
})

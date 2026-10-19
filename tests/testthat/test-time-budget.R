# max_time binds inside the E-step: both samplers are given the fit's deadline
# and stop drawing when it passes, between draws and inside a draw, and the
# MCEM drivers then end the fit with stop_reason "time_budget" and a note.

elapsed_of <- function(expr) {
  t0 <- proc.time()[[3L]]
  force(expr)
  proc.time()[[3L]] - t0
}

# A surviving constant-rate tree of 120 to 250 tips: large enough that an
# E-step of a few hundred draws takes seconds.
big_tree <- local({
  phy <- NULL
  function() {
    if (is.null(phy)) {
      set.seed(3)
      repeat {
        tr <- simulate_tree(pars = log(c(0.6, 0.1)), max_t = 9, model = "cr", max_tries = 50)
        if (tr$status == "done" && length(tr$tes$tip.label) > 120 &&
            length(tr$tes$tip.label) < 250) break
      }
      phy <<- tr$tes
    }
    phy
  }
})

small_brts <- c(5, 3, 2, 1)

test_that("the thinning E-step returns its partial sample, flagged, when the clock has run out", {
  r <- emphasis:::em_cpp(brts = small_brts, init_pars = c(log(0.5), 0, 0, 0, log(0.1), 0, 0, 0),
                         sample_size = 50L, maxN = 500L, max_missing = 1e4, max_lambda = 1e6,
                         lower_bound = rep(-10, 8), upper_bound = rep(10, 8), xtol_rel = 1e-3,
                         num_threads = 1L, copy_trees = FALSE, model = c(0L, 0L, 0L), link = 1L,
                         time_left = 0)
  expect_true(r$timed_out)
  expect_identical(r$num_trees, 0L)
  # no M-step ran: the iterate is where the E-step sampled
  expect_equal(r$estimates, c(log(0.5), 0, 0, 0, log(0.1), 0, 0, 0))
  # no deadline (the default) and a generous one leave the E-step as it was
  r2 <- emphasis:::em_cpp(brts = small_brts, init_pars = c(log(0.5), 0, 0, 0, log(0.1), 0, 0, 0),
                          sample_size = 50L, maxN = 500L, max_missing = 1e4, max_lambda = 1e6,
                          lower_bound = rep(-10, 8), upper_bound = rep(10, 8), xtol_rel = 1e-3,
                          num_threads = 1L, copy_trees = FALSE, model = c(0L, 0L, 0L), link = 1L,
                          time_left = 3600, seed = 7L)
  r3 <- emphasis:::em_cpp(brts = small_brts, init_pars = c(log(0.5), 0, 0, 0, log(0.1), 0, 0, 0),
                          sample_size = 50L, maxN = 500L, max_missing = 1e4, max_lambda = 1e6,
                          lower_bound = rep(-10, 8), upper_bound = rep(10, 8), xtol_rel = 1e-3,
                          num_threads = 1L, copy_trees = FALSE, model = c(0L, 0L, 0L), link = 1L,
                          seed = 7L)
  expect_false(r2$timed_out)
  expect_identical(r2$num_trees, 50L)
  expect_identical(r2$logf, r3$logf)
})

test_that("the surrogate sampler returns its partial sample, flagged, when the clock has run out", {
  e <- emphasis:::.augment_tree_bdi(small_brts, pars = c(log(0.5), -0.01, 0, 0, log(0.1), 0, 0, 0),
                                    model_bin = c(1L, 0L, 0L), sample_size = 20L, link = 1L,
                                    deadline = proc.time()[[3L]] - 1)
  expect_true(e$timed_out)
  expect_identical(e$n_valid, 0L)
  expect_length(e$trees, 0L)
  e2 <- emphasis:::.augment_tree_bdi(small_brts, pars = c(log(0.5), -0.01, 0, 0, log(0.1), 0, 0, 0),
                                     model_bin = c(1L, 0L, 0L), sample_size = 20L, link = 1L)
  expect_false(e2$timed_out)
  expect_identical(e2$n_valid, 20L)
})

test_that("a single thinning draw that runs long is stopped inside the draw", {
  skip_on_cran()
  phy <- big_tree()
  b <- emphasis:::.extract_brts(phy)
  # turnover 0.975 at speciation 20 over the whole crown age: one draw would
  # insert of the order of 1e5 missing lineages before max_missing = 1e7
  t <- elapsed_of(r <- emphasis:::em_cpp(
    brts = b, init_pars = c(log(20), 0, 0, 0, log(19.5), 0, 0, 0),
    sample_size = 1L, maxN = 10L, max_missing = 1e7, max_lambda = 1e6,
    lower_bound = rep(-10, 8), upper_bound = rep(10, 8), xtol_rel = 1e-3,
    num_threads = 1L, copy_trees = FALSE, model = c(0L, 0L, 0L), link = 1L,
    parent_tip_start = emphasis:::.pts(b), parent_id = emphasis:::.pid(b),
    time_left = 0.5))
  expect_true(r$timed_out)
  expect_lt(t, 3)
})

test_that("a single surrogate draw that runs long is stopped inside the draw", {
  skip_on_cran()
  phy <- big_tree()
  t <- elapsed_of(e <- emphasis:::.augment_tree_bdi(
    phy, pars = c(log(20), 0, 0, 0, log(19.5), 0, 0, 0), model_bin = c(0L, 0L, 0L),
    sample_size = 1L, max_missing = 1e7, link = 1L, deadline = proc.time()[[3L]] + 0.5))
  expect_true(e$timed_out)
  expect_lt(t, 3)
})

test_that("a surrogate fit on a slow E-step stops on its clock, inside the E-step", {
  skip_on_cran()
  phy <- big_tree()
  # 400 draws of the dd surrogate on this tree take several seconds
  t <- elapsed_of(fit <- suppressWarnings(estimate_rates(
    phy, method = "mcem", model = "dd",
    init_pars = c(log(0.6), -0.001, log(0.3), 0.001),
    control = list(lower_bound = c(-3, -0.05, -6, -0.05), upper_bound = c(1, 0.05, 0, 0.05),
                   num_trees = 400L, max_iter = 50L, max_time = 0.5))))
  expect_identical(fit$sampling, "bdi")
  expect_identical(fit$stop_reason, "time_budget")
  expect_lt(t, 6)                         # the E-step alone runs past 4 s
  expect_true("clock" %in% names(fit$notes))
  expect_match(fit$notes[["clock"]], "inside the E-step of iteration 1")
  expect_output(print(fit), "inside the E-step")
})

test_that("a thinning fit on a slow E-step stops on its clock, inside the E-step", {
  skip_on_cran()
  phy <- big_tree()
  # 4000 thinning draws on this tree take well over 10 s
  t <- elapsed_of(fit <- suppressWarnings(estimate_rates(
    phy, method = "mcem", model = "cr", init_pars = c(log(0.6), log(0.5)),
    control = list(lower_bound = c(-3, -6), upper_bound = c(1, 0), sampling = "dynamic_fresh",
                   num_trees = 4000L, max_iter = 50L, max_time = 0.5))))
  expect_identical(fit$stop_reason, "time_budget")
  expect_identical(fit$iterations, 0L)
  expect_lt(t, 5)
  expect_match(fit$notes[["clock"]], "inside the E-step of iteration 1")
})

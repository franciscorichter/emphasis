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

# The M-step binds on the same deadline: its optimiser reads the clock before
# every evaluation of the objective and stops once it has passed.  A
# conditional that sleeps makes each evaluation slow while the E-step on this
# four-tip tree stays fast, so the M-step alone runs well past the budget.
mstep_setup <- function() {
  mb <- emphasis:::.resolve_model("cr")
  ex <- function(p) emphasis:::.expand_pars(p, mb)
  list(mb = mb, start = ex(c(log(1.5), log(0.5))),
       lo = ex(c(-3, -6)), up = ex(c(1, 0)))
}

test_that("m_cpp stops on its clock and keeps the best point only when it beats the start", {
  s <- mstep_setup()
  e <- emphasis:::em_cpp(brts = small_brts, init_pars = s$start, sample_size = 30L, maxN = 500L,
                         max_missing = 1e4, max_lambda = 1e6, lower_bound = s$lo, upper_bound = s$up,
                         xtol_rel = 1e-3, num_threads = 1L, copy_trees = TRUE, model = s$mb,
                         link = 1L, seed = 3L)
  counter <- new.env()
  counter$n <- 0L
  counting <- function(p) { counter$n <- counter$n + 1L; 0 }
  full <- emphasis:::m_cpp(e, s$start, "rpd1", s$lo, s$up, 1e-3, 1L, model = s$mb, link = 1L,
                           rconditional = counting)
  expect_false(full$timed_out)
  # at 0.05 s an evaluation the unbounded M-step takes over five times the budget
  expect_gt(counter$n * 0.05, 5 * 0.3)
  sleepy <- function(p) { Sys.sleep(0.05); 0 }
  t <- elapsed_of(cut <- emphasis:::m_cpp(e, s$start, "rpd1", s$lo, s$up, 1e-3, 1L, model = s$mb,
                                          link = 1L, rconditional = sleepy, time_left = 0.3))
  expect_true(cut$timed_out)
  expect_lt(t, 1)
  # a few evaluations from a start this far from the optimum find a better point
  expect_true(cut$improved)
  expect_false(isTRUE(all.equal(cut$estimates, s$start)))
  # no time left: the optimiser does not start, and the start is returned
  none <- emphasis:::m_cpp(e, s$start, "rpd1", s$lo, s$up, 1e-3, 1L, model = s$mb, link = 1L,
                           rconditional = sleepy, time_left = 0)
  expect_true(none$timed_out)
  expect_false(none$improved)
  expect_equal(none$estimates, s$start)
})

test_that("a fit whose M-step alone exceeds max_time stops within a small multiple of it", {
  skip_on_cran()
  s <- mstep_setup()
  # about a hundred evaluations at 0.1 s each (counted in the test above):
  # the unbounded M-step alone takes about ten times max_time
  sleepy <- function(p) { Sys.sleep(0.1); 0 }
  common <- list(brts = small_brts, pars = s$start, sample_size = 10L, max_missing = 1e4,
                 lower_bound = s$lo, upper_bound = s$up, max_iter = 50L, xtol = 1e-3,
                 tol = 1e-3, patience = 3L, num_threads = 1L, conditional = sleepy,
                 model = s$mb, link = 1L, max_time = 1)
  t_bdi <- elapsed_of(fit_bdi <- do.call(emphasis:::.mcem_bdi, common))
  t_thn <- elapsed_of(fit_thn <- do.call(emphasis:::.mcem_dynamic_fresh,
                                         c(common, list(maxN = 2000L))))
  for (fit in list(fit_bdi, fit_thn)) {
    expect_identical(fit$stop_reason, "time_budget")
    # the cut iteration keeps the better point it had found, and the fit ends there
    expect_identical(fit$iterations, 1L)
    expect_match(fit$time_note, "inside the M-step of iteration 1; that iteration keeps")
    expect_false(isTRUE(all.equal(fit$pars, s$start)))
  }
  # the clock, the final E-step's own clock (whose em_cpp call runs a bounded
  # M-step under the thinning sampler) and one evaluation past each
  expect_lt(t_bdi, 3.5)
  expect_lt(t_thn, 3.5)
})

test_that("an M-step cut before it improves on its start discards the iteration", {
  skip_on_cran()
  s <- mstep_setup()
  # the first evaluation, at the start, outlasts the whole budget
  counter <- new.env()
  counter$n <- 0L
  slow_first <- function(p) {
    counter$n <- counter$n + 1L
    if (counter$n == 1L) Sys.sleep(0.8)
    0
  }
  t <- elapsed_of(fit <- emphasis:::.mcem_bdi(
    brts = small_brts, pars = s$start, sample_size = 30L, max_missing = 1e4,
    lower_bound = s$lo, upper_bound = s$up, max_iter = 50L, xtol = 1e-3, tol = 1e-3,
    patience = 3L, num_threads = 1L, conditional = slow_first, model = s$mb, link = 1L,
    max_time = 0.4))
  expect_identical(fit$stop_reason, "time_budget")
  expect_identical(fit$iterations, 0L)
  expect_equal(fit$pars, s$start)
  expect_match(fit$time_note, "inside the M-step of iteration 1 before it found a point better")
  expect_lt(t, 2.5)
})

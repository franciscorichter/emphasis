# Tests for the seeding of the C++ samplers (audit id H47).
#
# The thinning augmenter, Model::extinction_time and the forward simulator each
# held an engine seeded from the wall clock XOR the thread id.  set.seed() did
# not reach any of them -- 0 of 4 repeats identical, against 4 of 4 for the
# pure-R BDI path -- and a thread_local engine is copied byte for byte into a
# forked child, so the children of a primed parent drew one shared stream.
#
# Now every Rcpp entry point takes a `seed`; the R layer defaults it to
# sample.int(.Machine$integer.max, 1L), and a seed of 0 draws one from R's
# generator inside C++, so set.seed() reaches the sampler from every caller.
# Within a call each thread draws from its own substream of that seed
# (inst/include/model.hpp, namespace emphasis::rng).
#
# What is pinned here:
#   - two runs under one seed are identical, two runs under different seeds
#     are not, for augment_trees, em_cpp and the forward simulator;
#   - the extinction times the augmentation draws are seeded too, not only its
#     birth times;
#   - retries of a forward simulation are distinct trees, not one tree redrawn;
#   - two forked children of a primed parent draw differently.

b6 <- c(5, 3.4, 2.6, 1.7, 1.1, 0.6)

cr    <- c(0L, 0L, 0L)
pars8 <- function(lambda, mu) c(lambda, 0, 0, 0, mu, 0, 0, 0)

aug <- function(N = 12L, maxN = 4000L, num_threads = 1L, rho = 1, ...) {
  augment_trees(b6, pars8(0.6, 0.2), as.integer(N), as.integer(maxN),
                max_missing = 10000L, max_lambda = 500,
                num_threads = as.integer(num_threads),
                model = cr, link = 0L, rho = rho, ...)
}


test_that("set.seed fixes the thinning augmenter", {
  set.seed(1); a <- aug()
  set.seed(1); b <- aug()
  set.seed(2); c <- aug()

  expect_identical(a$logf, b$logf)
  expect_identical(a$logg, b$logg)
  expect_identical(a$trees, b$trees)
  expect_false(identical(a$logf, c$logf))
})


test_that("the seed argument fixes the thinning augmenter without set.seed", {
  a <- aug(seed = 99L)
  b <- aug(seed = 99L)
  c <- aug(seed = 100L)

  expect_identical(a$logf, b$logf)
  expect_identical(a$trees, b$trees)
  expect_false(identical(a$logf, c$logf))
})


test_that("the extinction times an augmentation draws are seeded too", {
  # Model::extinction_time held a static thread_local engine of its own.  It
  # now draws from the stream the augmentation is on, so the lifetimes repeat
  # with the birth times rather than independently of them.
  ext <- function(r) unlist(lapply(r$trees, function(df) df$t_ext))
  a <- aug(seed = 7L)
  b <- aug(seed = 7L)
  c <- aug(seed = 8L)

  expect_gt(length(ext(a)), 0L)
  expect_identical(ext(a), ext(b))
  expect_false(identical(ext(a), ext(c)))
})


test_that("the unsampled-extant draw at rho < 1 is seeded", {
  # The rho < 1 branch of Model::extinction_time makes one extra uniform draw
  # per augmented birth; it is on the same stream.
  a <- aug(rho = 0.8, seed = 31L)
  b <- aug(rho = 0.8, seed = 31L)
  c <- aug(rho = 0.8, seed = 32L)

  expect_identical(a$trees, b$trees)
  expect_false(identical(a$logg, c$logg))
})


test_that("set.seed fixes em_cpp", {
  fit <- function() {
    em_cpp(brts = b6, init_pars = pars8(0.6, 0.2),
           sample_size = 12L, maxN = 4000L, max_missing = 10000L,
           max_lambda = 500,
           lower_bound = c(0.01, 0, 0, 0, 0.00, 0, 0, 0),
           upper_bound = c(3.00, 0, 0, 0, 2.00, 0, 0, 0),
           xtol_rel = 0.01, num_threads = 1L, copy_trees = FALSE,
           model = cr, link = 0L, rho = 1)
  }
  set.seed(4); a <- fit()
  set.seed(4); b <- fit()
  set.seed(5); c <- fit()

  expect_identical(a$logf, b$logf)
  expect_identical(a$estimates, b$estimates)
  expect_false(identical(a$logf, c$logf))
})


test_that("set.seed fixes the forward simulator", {
  sim <- function() simulate_tree(pars = c(0.6, 0.1), max_t = 3,
                                  max_tries = 50, model = "cr",
                                  useDDD = FALSE)
  set.seed(1); a <- sim()
  set.seed(1); b <- sim()
  set.seed(7); c <- sim()

  expect_identical(a$L, b$L)
  expect_identical(a$status, b$status)
  expect_false(identical(a$L, c$L))
})


test_that("simulate_div_tree_cpp honours an explicit seed", {
  one <- function(s) simulate_div_tree_cpp(pars8(0.6, 0.1), cr, 3, 1000000L,
                                           0L, 0L, seed = s)
  expect_identical(one(11L)$Ltable, one(11L)$Ltable)
  expect_false(identical(one(11L)$Ltable, one(12L)$Ltable))
})


test_that("a forward simulation's retries are distinct trees", {
  # The retry loop in simulate_tree() draws a seed per attempt.  One seed for
  # the whole loop would redraw the same extinct tree until max_tries ran out.
  set.seed(3)
  sizes <- vapply(seq_len(8L), function(i)
    NROW(simulate_tree(pars = c(0.3, 0.25), max_t = 4, max_tries = 20,
                       model = "cr", useDDD = FALSE)$L), 0L)
  expect_gt(length(unique(sizes)), 1L)
})


test_that("forked children of a primed parent draw different trees", {
  skip_on_os("windows")            # mclapply forks only on unix
  # The parent draws first, so both engines have a definite state; a
  # thread_local engine would then be inherited byte for byte and the children
  # would return one shared stream.
  invisible(aug(seed = 123L))
  out <- parallel::mclapply(1:4, function(i) aug()$logf,
                            mc.cores = 2L, mc.preschedule = TRUE)
  expect_false(isTRUE(all.equal(out[[1]], out[[2]])))   # different children
  expect_false(isTRUE(all.equal(out[[1]], out[[3]])))   # same child, next task
  expect_false(isTRUE(all.equal(out[[3]], out[[4]])))
})


test_that("forked children of a primed parent simulate different trees", {
  skip_on_os("windows")
  invisible(simulate_tree(pars = c(0.6, 0.1), max_t = 3, max_tries = 50,
                          model = "cr", useDDD = FALSE))
  sizes <- unlist(parallel::mclapply(1:4, function(i)
    NROW(simulate_tree(pars = c(0.6, 0.1), max_t = 3, max_tries = 50,
                       model = "cr", useDDD = FALSE)$L), mc.cores = 2L))
  expect_gt(length(unique(sizes)), 1L)
})


test_that("num_threads > 1 draws from seeded substreams but is not yet reproducible", {
  # Each TBB worker draws from its own substream of the seed, so what a worker
  # draws is a function of the seed alone.  Which attempts a worker completes
  # before the sample fills, and in which order the results take the mutex, is
  # not: E_step (src/E_step.cpp) keeps the first sample_size trees to arrive.
  # Two runs under one seed therefore share most of their draws without being
  # identical.  Pinning this needs the collection order fixed in E_step
  # (audit ids H62, H63), which this change does not touch.
  a <- aug(num_threads = 2L, seed = 11L)
  b <- aug(num_threads = 2L, seed = 11L)
  pooled <- unique(c(a$logf, b$logf))
  expect_length(a$logf, 12L)
  expect_lt(length(pooled), 2L * length(a$logf))
})

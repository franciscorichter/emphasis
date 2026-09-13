.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis)); suppressPackageStartupMessages(library(testthat))
run <- function(code) { r <- ListReporter$new(); with_reporter(r, code); d <- as.data.frame(r$get_results()); c(failed = sum(d$failed), error = sum(d$error), skipped = sum(d$skipped), nexp = sum(d$nb)) }
brts_em <- c(0.8, 0.6, 0.4, 0.2); pars8 <- c(0.5, -0.01, 0.01, 0, 0.1, 0, 0, 0); lb8 <- c(0, -0.1, -0.1, 0, 0, 0, 0, 0); ub8 <- c(2, 0.1, 0.1, 0, 0.5, 0, 0, 0)
cat("original test-em body, skip dropped:\n"); print(run(test_that("orig", {
  result <- mc_loglik(brts = brts_em, pars = pars8, sample_size = 10, maxN = 100, max_missing = 1000, max_lambda = 500,
                      lower_bound = lb8, upper_bound = ub8, xtol_rel = 1e-3, num_threads = 1, model = c(1L,1L,0L))
  expect_type(result, "list") })))
for (i in 1:3) { cat("fix sketch em_cpp rewrite rep", i, ":\n"); print(run(test_that("fix-em", {
  result <- emphasis:::em_cpp(brts = brts_em, init_pars = pars8, sample_size = 10, maxN = 100, max_missing = 1000, max_lambda = 500,
                              lower_bound = lb8, upper_bound = ub8, xtol_rel = 1e-3, num_threads = 1, copy_trees = FALSE, model = c(1L,1L,0L))
  expect_type(result, "list"); expect_true("fhat" %in% names(result)); expect_length(result$logf, 10) }))) }
for (i in 1:3) { cat("fix sketch E-step rewrite rep", i, ":\n"); print(run(test_that("fix-aug", {
  brts <- sort(as.numeric(ape::branching.times(ape::rphylo(8, 0.5, 0))), decreasing = TRUE)
  pars8 <- emphasis:::.expand_pars(c(0.5, 0.1), c(0L, 0L, 0L))
  aug <- emphasis:::augment_trees(brts, pars8, sample_size = 20, maxN = 500, max_missing = 500, max_lambda = 50000, num_threads = 1, model = c(0L,0L,0L))
  expect_length(aug$logf, 20); expect_true(is.finite(emphasis:::.is_fhat(aug$logf, aug$logg))) }))) }

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis)); suppressPackageStartupMessages(library(testthat))
td <- tempfile("h96"); dir.create(td)
w <- function(name, txt) writeLines(txt, file.path(td, name))
w("test-orig-em.R", '
test_that("orig em", {
  brts <- c(0.8, 0.6, 0.4, 0.2); pars8 <- c(0.5, -0.01, 0.01, 0, 0.1, 0, 0, 0)
  lb8 <- c(0, -0.1, -0.1, 0, 0, 0, 0, 0); ub8 <- c(2, 0.1, 0.1, 0, 0.5, 0, 0, 0)
  result <- mc_loglik(brts = brts, pars = pars8, sample_size = 10, maxN = 100, max_missing = 1000, max_lambda = 500,
                      lower_bound = lb8, upper_bound = ub8, xtol_rel = 1e-3, num_threads = 1, model = c(1L,1L,0L))
  expect_type(result, "list"); expect_true("fhat" %in% names(result)) })')
w("test-orig-aug.R", '
test_that("orig aug", {
  set.seed(42); tree <- ape::rphylo(8, 0.5, 0); brts <- ape::branching.times(tree)
  result <- mc_loglik(brts = brts, pars = c(0.1, 0.5, -0.01, 0.01), sample_size = 1, maxN = 50, max_missing = 500, max_lambda = 50000,
                      lower_bound = c(0, 0, -0.1, -0.1), upper_bound = c(0.5, 2, 0.1, 0.1), xtol_rel = 1e-5, num_threads = 1)
  expect_type(result, "list"); expect_true("fhat" %in% names(result)) })')
w("test-fix-em.R", '
for (i in 1:3) test_that(paste("fix em", i), {
  brts <- c(0.8, 0.6, 0.4, 0.2); pars8 <- c(0.5, -0.01, 0.01, 0, 0.1, 0, 0, 0)
  lb8 <- c(0, -0.1, -0.1, 0, 0, 0, 0, 0); ub8 <- c(2, 0.1, 0.1, 0, 0.5, 0, 0, 0)
  result <- emphasis:::em_cpp(brts = brts, init_pars = pars8, sample_size = 10, maxN = 100, max_missing = 1000, max_lambda = 500,
                              lower_bound = lb8, upper_bound = ub8, xtol_rel = 1e-3, num_threads = 1, copy_trees = FALSE, model = c(1L,1L,0L))
  expect_type(result, "list"); expect_true("fhat" %in% names(result)); expect_length(result$logf, 10) })')
w("test-fix-aug.R", '
for (i in 1:3) test_that(paste("fix aug", i), {
  brts <- sort(as.numeric(ape::branching.times(ape::rphylo(8, 0.5, 0))), decreasing = TRUE)
  pars8 <- emphasis:::.expand_pars(c(0.5, 0.1), c(0L, 0L, 0L))
  aug <- emphasis:::augment_trees(brts, pars8, sample_size = 20, maxN = 500, max_missing = 500, max_lambda = 50000, num_threads = 1, model = c(0L,0L,0L))
  expect_length(aug$logf, 20); expect_true(is.finite(emphasis:::.is_fhat(aug$logf, aug$logg))) })')
res <- as.data.frame(test_dir(td, reporter = "silent", stop_on_failure = FALSE, package = "emphasis", load_package = "none"))
print(res[, c("file", "test", "nb", "failed", "skipped", "error", "real")])

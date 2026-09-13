# H97 follow-up: (a) per-block error messages of the un-skipped test-em.R,
# (b) does the file pass once bounds are moved into control (the proposed fix)?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
suppressPackageStartupMessages(library(testthat))

src <- readLines("/Users/pancho/Code/emphasis/tests/testthat/test-em.R")
src <- src[!grepl("^\\s*skip\\(", src)]
tmp <- file.path(tempdir(), "test-em-unskipped.R")
writeLines(src, tmp)
cat("== error messages per block (un-skipped, as written) ==\n")
res <- testthat::test_file(tmp, reporter = "silent")
for (b in res) {
  errs <- Filter(function(x) inherits(x, "expectation_error") || inherits(x, "expectation_failure"), b$results)
  for (e in errs) cat(sprintf("[%s] %s: %s\n", b$test, class(e)[1], substr(gsub("\n", " ", conditionMessage(e)), 1, 160)))
}

cat("\n== control aliases: are tol/burnin accepted keys? ==\n")
print(names(emphasis:::estimate_rates_control("mcem")))
cat("has 'tol':", "tol" %in% names(emphasis:::estimate_rates_control("mcem")),
    " has 'burnin':", "burnin" %in% names(emphasis:::estimate_rates_control("mcem")), "\n")

cat("\n== proposed fix: bounds into control, mc_loglik block dropped ==\n")
fixed <- c(
'test_that("estimate_rates CR smoke test", {',
'  set.seed(42)',
'  tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")',
'  fit <- estimate_rates(tr, method = "mcem", model = "cr",',
'    control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),',
'                   sample_size = 20, max_iter = 3, max_time = 60, num_threads = 1))',
'  expect_named(fit$pars, c("beta_0", "gamma_0"))',
'  expect_length(fit$pars, 2)',
'  expect_true(is.numeric(fit$loglik))',
'})',
'test_that("estimate_rates DD smoke test", {',
'  set.seed(42)',
'  tr <- simulate_tree(pars = c(0.5, -0.005, 0.1, 0), max_t = 8, model = "dd")',
'  fit <- estimate_rates(tr, method = "mcem", model = "dd",',
'    control = list(lower_bound = c(0.1, -0.1, 0, -0.01), upper_bound = c(2, 0.01, 0.5, 0.01),',
'                   sample_size = 20, max_iter = 3, max_time = 60, num_threads = 1))',
'  expect_named(fit$pars, c("beta_0", "beta_N", "gamma_0", "gamma_N"))',
'  expect_length(fit$pars, 4)',
'})',
'test_that("estimate_rates EP smoke test", {',
'  set.seed(1)',
'  tr <- simulate_tree(pars = c(0.5, 0.05, 0.1, 0.01), max_t = 5, model = "ep")',
'  fit <- estimate_rates(tr, method = "mcem", model = "ep",',
'    control = list(lower_bound = c(0.1, -0.5, 0, -0.5), upper_bound = c(2, 0.5, 0.5, 0.5),',
'                   sample_size = 20, max_iter = 3, max_time = 60, num_threads = 1))',
'  expect_named(fit$pars, c("beta_0", "beta_D", "gamma_0", "gamma_D"))',
'  expect_true(is.numeric(fit$loglik))',
'})',
'test_that("estimate_rates EP + exponential link smoke test", {',
'  set.seed(1)',
'  tr <- simulate_tree(pars = c(0.5, 0.05, 0.1, 0.01), max_t = 5, model = "ep")',
'  fit <- estimate_rates(tr, method = "mcem", model = "ep", link = "exponential",',
'    control = list(lower_bound = c(-5, -5, -5, -5), upper_bound = c(2, 2, 2, 2),',
'                   sample_size = 20, max_iter = 3, max_time = 60, num_threads = 1))',
'  expect_named(fit$pars, c("beta_0", "beta_D", "gamma_0", "gamma_D"))',
'  expect_true(is.numeric(fit$loglik))',
'})')
tmp2 <- file.path(tempdir(), "test-em-fixed.R")
writeLines(fixed, tmp2)
res2 <- testthat::test_file(tmp2, reporter = "silent")
df <- as.data.frame(res2)
print(df[, c("test", "nb", "failed", "error")])
for (b in res2) {
  errs <- Filter(function(x) inherits(x, "expectation_error") || inherits(x, "expectation_failure"), b$results)
  for (e in errs) cat(sprintf("[%s] %s: %s\n", b$test, class(e)[1], substr(gsub("\n", " ", conditionMessage(e)), 1, 200)))
}

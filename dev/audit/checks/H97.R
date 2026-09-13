# H97: tests/testthat/test-em.R calls estimate_rates(..., lower_bound=, upper_bound=)
# which are not formals. Decide whether the calls would fail if un-skipped and
# whether moving bounds into control fixes them.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
suppressPackageStartupMessages(library(testthat))

cat("== formals(estimate_rates) ==\n")
print(names(formals(estimate_rates)))
cat("has '...':", "..." %in% names(formals(estimate_rates)), "\n")
cat("has lower_bound formal:", "lower_bound" %in% names(formals(estimate_rates)), "\n")

cat("\n== exact test-em.R:33-36 call (CR smoke) ==\n")
set.seed(42)
tr <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
r <- tryCatch(
  estimate_rates(tr, method = "mcem", model = "cr",
                 lower_bound = c(0, 0), upper_bound = c(2, 1),
                 control = list(sample_size = 20, tol = 0.5, burnin = 2)),
  error = function(e) e)
cat("class:", paste(class(r), collapse = "/"), "\n")
if (inherits(r, "error")) cat("message:", conditionMessage(r), "\n")

cat("\n== same call with bounds moved into control (small budget) ==\n")
r2 <- tryCatch(
  estimate_rates(tr, method = "mcem", model = "cr",
                 control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                                sample_size = 20, max_iter = 3, max_time = 60,
                                num_threads = 1)),
  error = function(e) e)
cat("class:", paste(class(r2), collapse = "/"), "\n")
if (inherits(r2, "error")) cat("message:", conditionMessage(r2), "\n") else {
  cat("names(fit$pars):", paste(names(r2$pars), collapse = ","), "\n")
  cat("loglik numeric:", is.numeric(r2$loglik), "\n")
}

cat("\n== other stale pieces in the same blocks ==\n")
cat("mc_loglik exists:", exists("mc_loglik", envir = asNamespace("emphasis")), "\n")
cat(".resolve_model('ep'):", tryCatch(paste(emphasis:::.resolve_model("ep"), collapse=","),
                                     error = function(e) paste("ERROR:", conditionMessage(e))), "\n")
cat("simulate_tree dd 4 pars:", tryCatch({
  s <- simulate_tree(pars = c(0.5, -0.005, 0.1, 0), max_t = 3, model = "dd"); "ok"
}, error = function(e) paste("ERROR:", conditionMessage(e))), "\n")

cat("\n== run test-em.R with every skip() removed ==\n")
src <- readLines("/Users/pancho/Code/emphasis/tests/testthat/test-em.R")
src <- src[!grepl("^\\s*skip\\(", src)]
tmp <- file.path(tempdir(), "test-em-unskipped.R")
writeLines(src, tmp)
res <- as.data.frame(testthat::test_file(tmp, reporter = "silent"))
print(res[, c("test", "nb", "failed", "error", "skipped")])

# H98: tests/testthat/test-inference.R:74,88,101 call
#   simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr")
# which binds c(0.5, 0.1) to the first formal `tree`, leaving `pars` missing.
# Self-contained check; no compilation.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(testthat) })

cat("formals(simulate_tree)[1:2]:", paste(names(formals(simulate_tree))[1:2], collapse = ", "), "\n\n")

# 1. The positional call as written in the tests / examples
r1 <- tryCatch(simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr"),
               error = function(e) conditionMessage(e))
cat("1. positional call  -> ", if (is.character(r1)) paste("ERROR:", r1) else "returned a list", "\n")
stopifnot(is.character(r1), grepl("pars", r1))

# 2. The named form works
r2 <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr")
cat("2. named  pars= call -> status:", r2$status, "| class(tes):", class(r2$tes), "\n")
stopifnot(is.list(r2), r2$status %in% c("done", "extinct", "overflow"))

# 3. The cited test body (test-inference.R:70-82) with the skip() lines removed.
#    It must fail at the simulate_tree line, before estimate_rates is ever reached.
res <- tryCatch(
  testthat::with_reporter("silent",
    testthat::test_that("H98 cited test body without skip", {
      set.seed(42)
      sim <- simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr")
      fit <- estimate_rates(sim, method = "mcem", model = "cr",
        control = list(lower_bound = c(0, 0), upper_bound = c(2, 1),
                       max_iter = 5, sample_size = 50))
      expect_s3_class(fit, "emphasis_fit")
    })),
  error = function(e) FALSE)
cat("3. test body without skip passes?", isTRUE(res), "\n")
stopifnot(!isTRUE(res))

# 4. The shipped test file as-is: the three tests are skipped, so they report SKIP, never FAIL.
tf <- "/Users/pancho/Code/emphasis/tests/testthat/test-inference.R"
out <- testthat::test_file(tf, reporter = "silent")
df  <- as.data.frame(out)
cat("\n4. shipped test-inference.R: failed =", sum(df$failed), " skipped =", sum(df$skipped),
    " passed =", sum(df$passed), "\n")
print(df[grepl("end-to-end", df$test), c("test", "failed", "skipped", "passed")])

# 5. Same positional form in roxygen examples (all inside \dontrun, so R CMD check never runs them)
ex <- system("grep -n 'simulate_tree(c(' /Users/pancho/Code/emphasis/R/*.R /Users/pancho/Code/emphasis/man/*.Rd", intern = TRUE)
cat("\n5. positional-form occurrences in R/ examples and man/:", length(ex), "\n")
cat(paste0("   ", ex), sep = "\n")

# 6. Corrected form for each documented example also runs
r6a <- simulate_tree(pars = c(0.5, -0.005, 0.1, 0), max_t = 8, model = "dd")
r6b <- simulate_tree(pars = c(0.6, -0.05, 0.15, -0.03), max_t = 10, model = "d")
cat("\n6. named dd example status:", r6a$status, " | named d example status:", r6b$status, "\n")
cat("\nH98: CONFIRMED (latent test/doc defect; dead code because of skip()/dontrun)\n")

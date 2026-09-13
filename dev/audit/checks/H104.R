.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis); library(testthat)
cat("NOT_CRAN env:", Sys.getenv("NOT_CRAN", "<unset>"), "\n")
show <- function(file) {
  r <- testthat::test_file(file, reporter = "silent", package = "emphasis")
  df <- as.data.frame(r)
  for (i in seq_len(nrow(df))) if (df$skipped[i])
    cat(sprintf("SKIPPED '%s': %s\n", df$test[i], conditionMessage(r[[i]]$results[[1]])))
}
show("/Users/pancho/Code/emphasis/tests/testthat/test-inference.R")
show("/Users/pancho/Code/emphasis/tests/testthat/test-augment.R")
cat("--- now with NOT_CRAN=true ---\n"); Sys.setenv(NOT_CRAN = "true")
show("/Users/pancho/Code/emphasis/tests/testthat/test-inference.R")
show("/Users/pancho/Code/emphasis/tests/testthat/test-augment.R")

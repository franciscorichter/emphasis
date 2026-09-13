lib <- commandArgs(trailingOnly = TRUE)[1]
Sys.setenv(TESTTHAT_MAX_FAILS = "1000")
.libPaths(c(lib, .libPaths()))
library(emphasis)
e <- new.env(parent = asNamespace("emphasis"))
r <- testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-bdi-dd.R",
                    env = e, reporter = "silent")
df <- as.data.frame(r)
cat("build:", lib, "\n")
print(df[, c("test", "nb", "failed", "error", "warning", "skipped")])

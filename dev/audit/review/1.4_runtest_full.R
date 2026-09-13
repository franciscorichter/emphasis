args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
library(emphasis)
options(testthat.summary.max_reports = Inf)
r <- testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-bdi-dd.R",
                    env = new.env(parent = asNamespace("emphasis")), reporter = "silent")
df <- as.data.frame(r)
print(df[, c("test","nb","failed","skipped","error","warning")])

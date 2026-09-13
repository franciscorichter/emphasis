lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths()))
options(testthat.progress.max_fails = 1000)
suppressMessages(library(emphasis))
cat("build:", lib, "\n")
res <- testthat::test_file("/Users/pancho/Code/emphasis/dev/audit/review/1.3_test-bdi-cr.R",
                    env = new.env(parent = asNamespace("emphasis")),
                    reporter = "silent")
d <- as.data.frame(res)
print(d[, c("test","nb","failed","error","skipped","warning")])
cat("TOTAL failed:", sum(d$failed), " errors:", sum(d$error), " passed:", sum(d$passed), "\n")

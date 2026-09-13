args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
library(emphasis)
cat("emphasis loaded from:", dirname(system.file(package = "emphasis")), "\n")
res <- testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-mcem-thinning.R",
                           env = new.env(parent = asNamespace("emphasis")), reporter = "summary")
df <- as.data.frame(res)
print(df[, c("test", "nb", "failed", "skipped", "error", "real")])

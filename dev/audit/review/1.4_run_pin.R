lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, .libPaths()))
library(emphasis)
cat("build:", lib, "\n")
cat("desc:", as.character(packageVersion("emphasis")), "\n")
e <- new.env(parent = asNamespace("emphasis"))
testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-bdi-dd.R",
                    env = e, reporter = "summary")

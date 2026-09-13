lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths()))
library(emphasis)
cat("lib:", lib, "\n")
testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-estep.R",
                    env = new.env(parent = asNamespace("emphasis")),
                    reporter = "summary")

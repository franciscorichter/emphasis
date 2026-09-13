lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths()))
library(emphasis)
cat("lib:", lib, " version:", as.character(packageVersion("emphasis")), "\n")
testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-mcem-bdi.R",
                    env = new.env(parent = asNamespace("emphasis")),
                    reporter = "summary")

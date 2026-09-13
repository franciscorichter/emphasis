args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
library(emphasis)
cat("build:", args[1], "\n")
testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-bdi-dd.R",
                    env = new.env(parent = asNamespace("emphasis")), reporter = "summary")

lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths()))
suppressMessages(library(emphasis))
cat("build:", lib, " pkg ver:", as.character(packageVersion("emphasis")), "\n")
testthat::test_file("/Users/pancho/Code/emphasis/dev/audit/review/1.3_test-bdi-cr.R",
                    env = new.env(parent = asNamespace("emphasis")),
                    reporter = "summary")

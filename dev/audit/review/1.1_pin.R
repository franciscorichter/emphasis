args <- commandArgs(TRUE); lib <- args[1]
.libPaths(c(lib, .libPaths())); suppressMessages(library(emphasis))
cat("lib:", lib, " pkg:", dirname(system.file(package="emphasis")), "\n")
testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-logsum.R",
  env = new.env(parent = asNamespace("emphasis")), reporter = "summary")

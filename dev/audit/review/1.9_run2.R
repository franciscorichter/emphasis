lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths()))
options(testthat.progress_max_fails = 200)
library(emphasis)
cat("### lib:", lib, "\n")
testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-estimate-rates.R",
  env = new.env(parent = asNamespace("emphasis")), reporter = "summary")

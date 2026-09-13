lib <- commandArgs(trailingOnly=TRUE)[1]
.libPaths(c(lib, .libPaths()))
library(emphasis)
cat("### lib:", lib, "\n")
rep <- testthat::SummaryReporter$new(max_reports = 100L)
testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-estimate-rates.R",
  env = new.env(parent = asNamespace("emphasis")), reporter = rep)

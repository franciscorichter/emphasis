args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages(library(emphasis))
cat("emphasis loaded from:", find.package("emphasis"), "\n")
options(testthat.summary.max_reports = Inf, testthat.progress.max_fails = Inf)
testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-estep.R",
                    env = new.env(parent = asNamespace("emphasis")),
                    reporter = "summary")

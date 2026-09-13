args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(emphasis))
cat("emphasis loaded from:", find.package("emphasis"), "\n")
cat("has .n0:", exists(".n0", envir = asNamespace("emphasis")), "\n")
cat("formals tol:", deparse(formals(emphasis:::.mcem_bdi)$tol), "\n")
set.seed(1)
res <- testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-mcem-bdi.R",
  env = new.env(parent = asNamespace("emphasis")), reporter = "summary")
df <- as.data.frame(res)
print(df[, c("test", "nb", "failed", "skipped", "error", "warning")])

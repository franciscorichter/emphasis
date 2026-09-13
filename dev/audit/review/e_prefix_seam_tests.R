# Discrimination: the seam-relevant wave-1 test files against a given build.
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
library(emphasis)
files <- c("test-mcem-bdi.R", "test-mcem-thinning.R", "test-estimate-rates.R",
           "test-mstep.R", "test-logsum.R", "test-estep.R")
for (f in files) {
  p <- file.path("/Users/pancho/Code/emphasis/tests/testthat", f)
  cat("\n===== ", f, " =====\n", sep = "")
  r <- tryCatch(
    testthat::test_file(p, env = new.env(parent = asNamespace("emphasis")),
                        reporter = "summary", stop_on_failure = FALSE),
    error = function(e) { cat("FILE ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) {
    d <- as.data.frame(r)
    cat(sprintf("[%s] pass=%d fail=%d error=%d skip=%d\n", f,
                sum(d$passed), sum(d$failed), sum(d$error), sum(d$skipped)))
  }
}

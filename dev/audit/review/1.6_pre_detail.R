args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
library(emphasis)
res <- testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-mcem-thinning.R",
                           env = new.env(parent = asNamespace("emphasis")), reporter = testthat::ListReporter$new())
for (t in res) {
  cat("\n## ", t$test, "\n")
  for (r in t$results) {
    cls <- class(r)[1]
    loc <- if (!is.null(r$srcref)) as.integer(r$srcref[1]) else NA
    cat(sprintf("  line %s  %-22s %s\n", loc, cls, substr(gsub("\n", " ", conditionMessage(r)), 1, 90)))
  }
}

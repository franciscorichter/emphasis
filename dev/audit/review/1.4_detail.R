lib <- commandArgs(trailingOnly = TRUE)[1]
Sys.setenv(TESTTHAT_MAX_FAILS = "1000")
.libPaths(c(lib, .libPaths())); library(emphasis)
e <- new.env(parent = asNamespace("emphasis"))
r <- testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-bdi-dd.R",
                    env = e, reporter = "silent")
for (tr in r) {
  cat("\n### ", tr$test, "\n")
  for (res in tr$results) {
    cls <- class(res)[1]
    loc <- if (!is.null(res$srcref)) paste0("L", res$srcref[1]) else "?"
    msg <- if (cls == "expectation_success") "" else gsub("\n", " | ", substr(conditionMessage(res), 1, 130))
    cat(sprintf("  %-26s %-6s %s\n", cls, loc, msg))
  }
}

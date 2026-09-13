args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
suppressWarnings(library(emphasis))
res <- testthat::test_file("/Users/pancho/Code/emphasis/dev/audit/review/1.3_test-bdi-cr.R",
                           env = new.env(parent = asNamespace("emphasis")),
                           reporter = "silent")
df <- as.data.frame(res)
cat(sprintf("BUILD %s\n", args[1]))
print(df[, c("test", "nb", "failed", "skipped", "error", "warning")], row.names = FALSE)
cat(sprintf("TOTAL nb=%d failed=%d error=%d skipped=%d\n", sum(df$nb), sum(df$failed), sum(df$error), sum(df$skipped)))

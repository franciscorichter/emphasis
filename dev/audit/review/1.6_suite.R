args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
library(emphasis)
files <- list.files("/Users/pancho/Code/emphasis/tests/testthat", pattern = "^test-.*\\.R$", full.names = TRUE)
for (f in files) {
  t0 <- proc.time()[3]
  res <- tryCatch(testthat::test_file(f, env = new.env(parent = asNamespace("emphasis")), reporter = "silent"), error = function(e) e)
  if (inherits(res, "error")) { cat(sprintf("%-40s ERROR %s\n", basename(f), conditionMessage(res))); next }
  df <- as.data.frame(res)
  cat(sprintf("%-40s tests=%2d exp=%3d failed=%2d error=%2d skipped=%2d  %.1fs\n", basename(f), nrow(df), sum(df$nb), sum(df$failed), sum(df$error), sum(df$skipped), proc.time()[3] - t0))
  bad <- df[df$failed > 0 | df$error, ]
  if (nrow(bad)) for (i in seq_len(nrow(bad))) cat("    FAIL:", bad$test[i], "\n")
}

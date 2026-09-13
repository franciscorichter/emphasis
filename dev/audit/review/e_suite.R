# (e) whole suite against the POST-FIX build.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis)
options(warn = 1)
t0 <- Sys.time()
res <- testthat::test_dir("/Users/pancho/Code/emphasis/tests/testthat",
                          env = new.env(parent = asNamespace("emphasis")),
                          reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat("\n#### TOTALS ####\n")
cat("files  :", nrow(df), "\n")
cat("PASS   :", sum(df$passed), "\n")
cat("FAIL   :", sum(df$failed), "\n")
cat("ERROR  :", sum(df$error), "\n")
cat("SKIP   :", sum(df$skipped), "\n")
cat("elapsed:", round(difftime(Sys.time(), t0, units = "secs")), "s\n")
print(df[df$failed > 0 | df$error, c("file", "passed", "failed", "skipped", "error")])

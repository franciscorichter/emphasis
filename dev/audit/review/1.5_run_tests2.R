lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, .libPaths()))
options(testthat.progress.max_fails = 1000)
library(emphasis)
res <- testthat::test_file("/Users/pancho/Code/emphasis/tests/testthat/test-mcem-bdi.R",
                    env = new.env(parent = asNamespace("emphasis")),
                    reporter = "silent")
df <- as.data.frame(res)
for (i in seq_len(nrow(df))) cat(sprintf("%-95s pass=%d fail=%d err=%s\n",
   substr(df$test[i],1,95), df$passed[i], df$failed[i], df$error[i]))

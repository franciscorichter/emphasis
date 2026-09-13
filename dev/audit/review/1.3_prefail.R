.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
options(testthat.progress.max_fails=1000)
suppressMessages(library(emphasis))
res <- testthat::test_file("/Users/pancho/Code/emphasis/dev/audit/review/1.3_test-bdi-cr.R",
        env=new.env(parent=asNamespace("emphasis")), reporter="silent")
for (f in res) {
  for (r in f$results) {
    cl <- class(r)[1]
    if (grepl("failure|error", cl)) cat(sprintf("[%s] %s :: %s\n", cl, f$test,
        gsub("\n"," | ", substr(conditionMessage(r),1,160))))
  }
}

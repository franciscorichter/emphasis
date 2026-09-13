.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({ library(emphasis); library(testthat) })
src <- "/Users/pancho/Code/emphasis/tests/testthat"
tmp <- file.path(tempdir(), "H94err"); dir.create(tmp, showWarnings = FALSE)
for (f in list.files(src, pattern = "^test-.*\\.R$", full.names = TRUE)) {
  l <- readLines(f); l <- l[!grepl("^\\s*skip\\(\"", l)]
  writeLines(l, file.path(tmp, basename(f)))
}
res <- testthat::test_dir(tmp, package = "emphasis", load_package = "installed",
                          reporter = "silent", stop_on_failure = FALSE)
for (r in res) for (e in r$results) if (!inherits(e, "expectation_success"))
  cat(sprintf("%s | %s | %s | %s\n", r$file, r$test, paste(class(e)[1:2], collapse="/"),
              substr(gsub("\n", " ", conditionMessage(e)), 1, 220)))
cat("exports augment_trees/eval_logf/mc_loglik:",
    c("augment_trees","eval_logf","mc_loglik") %in% getNamespaceExports("emphasis"), "\n")
ns <- asNamespace("emphasis")
aug <- ns$augment_trees(c(4,2.5,1.2,0.6), c(0.5,0,0,0,0.1,0,0,0), sample_size = 3L, maxN = 200L,
        max_missing = 1000L, max_lambda = 100, num_threads = 1L, model = c(0L,0L,0L), link = 0L)
ev <- ns$eval_logf(pars = c(0.5,0,0,0,0.1,0,0,0), trees = aug$trees, model = c(0L,0L,0L), link = 0L)
cat("eval_logf names:", paste(names(ev), collapse=","), " logf match:", isTRUE(all.equal(ev$logf, aug$logf)), "\n")

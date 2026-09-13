## H94 replication — independent variation: strip skip() lines, run each file
## with test_file (not test_dir), and print the *error messages* of every
## non-passing re-enabled block; also verify each claimed drift category by
## direct calls against the installed API.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({ library(emphasis); library(testthat) })
src <- "/Users/pancho/Code/emphasis/tests/testthat"
tmp <- file.path(tempdir(), "H94rep"); dir.create(tmp, showWarnings = FALSE)
for (f in list.files(src, pattern = "^test-.*\\.R$", full.names = TRUE)) {
  l <- readLines(f); l <- l[!grepl("^\\s*skip\\(\"", l)]
  writeLines(l, file.path(tmp, basename(f)))
}
for (f in list.files(tmp, pattern = "^test-.*\\.R$", full.names = TRUE)) {
  res <- as.data.frame(testthat::test_file(f, package = "emphasis",
                        load_package = "installed", reporter = "silent"))
  for (i in seq_len(nrow(res))) {
    if (res$failed[i] > 0 || res$error[i]) {
      for (e in res$result[[i]]) if (!inherits(e, "expectation_success"))
        cat(sprintf("%s | %s | %s | %s\n", basename(f), res$test[i], class(e)[1],
                    substr(gsub("\n", " ", conditionMessage(e)), 1, 200)))
    }
  }
}
cat("\n--- direct API checks ---\n")
cat("mc_loglik exists:", exists("mc_loglik", envir = asNamespace("emphasis")), "\n")
cat("estimate_rates formals:", paste(names(formals(emphasis::estimate_rates)), collapse = ","), "\n")
cat("simulate_tree formals:", paste(names(formals(emphasis::simulate_tree)), collapse = ","), "\n")
tr <- tryCatch(simulate_tree(c(0.5, 0.1), max_t = 5, model = "cr"), error = function(e) conditionMessage(e)); print(tr)
tr2 <- simulate_tree(pars = c(0.5, 0.1), max_t = 5, model = "cr"); cat("named pars ok, class:", class(tr2), "\n")
aug <- augment_trees(c(4,2.5,1.2,0.6), c(0.5,0,0,0,0.1,0,0,0), sample_size = 3L, maxN = 200L,
                     max_missing = 1000L, max_lambda = 100, num_threads = 1L, model = c(0L,0L,0L), link = 0L)
ev <- eval_logf(pars = c(0.5,0,0,0,0.1,0,0,0), trees = aug$trees, model = c(0L,0L,0L), link = 0L)
cat("eval_logf names:", paste(names(ev), collapse = ","), "\n")
cat("eval_logf$logf == aug$logf:", isTRUE(all.equal(ev$logf, aug$logf)), "\n")

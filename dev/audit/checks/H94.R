## H94 — are all C++-reaching tests unconditionally skip()-ed, and would the
## skipped bodies still pass against the current build if re-enabled?
##
## Part A: run the shipped test suite against the installed build; count
##         pass / skip / fail per file, and count how many calls reach any of
##         the five .Call entry points (RcppExports) during the non-skipped tests.
## Part B: strip the unconditional skip() lines into a temp copy of the test
##         directory and run it again -> which of the skipped bodies pass under
##         the current signatures (signature-drift check).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
Sys.setenv(NOT_CRAN = "true")  # so skip_on_cran() does not mask the result
suppressPackageStartupMessages({ library(emphasis); library(testthat) })

src_dir <- "/Users/pancho/Code/emphasis/tests/testthat"
cpp_entry <- c("simulate_div_tree_cpp", "eval_logf", "augment_trees", "em_cpp", "m_cpp")
cpp_entry <- intersect(cpp_entry, ls(asNamespace("emphasis")))
cat("C++ entry points in namespace:", paste(cpp_entry, collapse = ", "), "\n")

## --- static count ---------------------------------------------------------
files <- list.files(src_dir, pattern = "^test-.*\\.R$", full.names = TRUE)
static <- do.call(rbind, lapply(files, function(f) {
  l <- readLines(f)
  data.frame(file = basename(f),
             n_test_that = sum(grepl("^test_that\\(", l)),
             n_uncond_skip = sum(grepl("^\\s*skip\\(\"", l)),
             stringsAsFactors = FALSE)
}))
print(static); cat("TOTAL test_that:", sum(static$n_test_that),
                   " unconditional skip():", sum(static$n_uncond_skip), "\n\n")

## --- Part A: run as shipped, count C++ reaches ----------------------------
.cpp_hits <- new.env(); .cpp_hits$n <- 0L; .cpp_hits$by <- character()
for (fn in cpp_entry) {
  suppressMessages(trace(fn, where = asNamespace("emphasis"), print = FALSE,
        tracer = bquote({ .cpp_hits$n <- .cpp_hits$n + 1L
                          .cpp_hits$by <- c(.cpp_hits$by, .(fn)) })))
}
run_suite <- function(dir) {
  res <- testthat::test_dir(dir, package = "emphasis", load_package = "installed",
                            reporter = "silent", stop_on_failure = FALSE)
  df <- as.data.frame(res)
  agg <- aggregate(cbind(nb, failed, skipped, error) ~ file, data = df, FUN = sum)
  agg$n_tests <- as.vector(table(df$file)[agg$file])
  list(df = df, agg = agg)
}
A <- run_suite(src_dir)
cat("=== Part A: suite as shipped (installed build) ===\n")
print(A$agg)
cat(sprintf("tests: %d  skipped: %d  failed: %d  errored: %d\n",
            nrow(A$df), sum(A$df$skipped), sum(A$df$failed > 0), sum(A$df$error)))
cat(sprintf("C++ entry-point calls during shipped suite: %d  (%s)\n\n",
            .cpp_hits$n, if (.cpp_hits$n) paste(names(table(.cpp_hits$by)),
                                                 table(.cpp_hits$by), collapse = ", ") else "none"))
skipped_names <- A$df$test[A$df$skipped]

## --- Part B: unskip and re-run --------------------------------------------
tmp <- file.path(tempdir(), "H94_unskipped"); dir.create(tmp, showWarnings = FALSE)
for (f in files) {
  l <- readLines(f)
  l <- l[!grepl("^\\s*skip\\(\"", l)]
  writeLines(l, file.path(tmp, basename(f)))
}
.cpp_hits$n <- 0L; .cpp_hits$by <- character()
t0 <- Sys.time()
B <- run_suite(tmp)
cat("=== Part B: same suite with unconditional skip() lines removed ===\n")
print(B$agg)
cat(sprintf("elapsed: %.0f s; C++ entry-point calls: %d\n",
            as.numeric(Sys.time() - t0, units = "secs"), .cpp_hits$n))
prev <- B$df[B$df$test %in% skipped_names, c("file", "test", "nb", "failed", "error", "skipped")]
cat("\nPreviously-skipped blocks, status when re-enabled:\n")
print(prev, row.names = FALSE)
cat(sprintf("\nre-enabled: %d  pass: %d  fail: %d  error: %d\n",
            nrow(prev), sum(prev$failed == 0 & !prev$error & !prev$skipped),
            sum(prev$failed > 0), sum(prev$error)))

## failure messages for the re-enabled blocks that did not pass
bad <- prev$test[prev$failed > 0 | prev$error]
if (length(bad)) {
  cat("\n--- failure details ---\n")
  for (i in which(B$df$test %in% bad)) for (e in B$df$result[[i]])
    if (!inherits(e, "expectation_success"))
      cat(sprintf("[%s] %s :: %s\n", B$df$test[i], class(e)[1],
                  substr(gsub("\n", " | ", conditionMessage(e)), 1, 300)))
}

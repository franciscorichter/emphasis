## H96_replicate -- independent replication, varying what H96.R did not:
##  a different tree (20 tips), sample_size = 1 (test-augment's own value),
##  brts exactly as branching.times() returns them (unsorted, named),
##  the exponential link, and the fix sketch's augment_trees + .is_fhat block.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
suppressPackageStartupMessages(library(testthat))

cat("== A. is mc_loglik anywhere on any libPath / any installed emphasis? ==\n")
for (lp in .libPaths()) {
  d <- file.path(lp, "emphasis")
  if (dir.exists(d)) {
    v <- tryCatch(read.dcf(file.path(d, "DESCRIPTION"), "Version")[1], error = function(e) NA)
    ns <- tryCatch({ e <- new.env(); lazyLoad(file.path(d, "R", "emphasis"), envir = e); e },
                   error = function(e) NULL)
    cat(sprintf("  %s : emphasis %s : mc_loglik bound = %s\n", lp, v,
                if (is.null(ns)) "?" else exists("mc_loglik", envir = ns, inherits = FALSE)))
  } else cat(sprintf("  %s : no emphasis\n", lp))
}
cat("  attached namespace path:", getNamespaceInfo("emphasis", "path"), "\n")
cat("  packageVersion:", as.character(packageVersion("emphasis")), "\n")
cat("  exists('mc_loglik') from globalenv:", exists("mc_loglik"), "\n")

cat("\n== B. verbatim test-augment body, new tree, sample_size = 1, brts as branching.times() returns ==\n")
set.seed(7)
tree <- ape::rphylo(20, 0.6, 0.1)
brts <- ape::branching.times(tree)              # named, node order, NOT sorted
cat("  brts sorted decreasing as given?", !is.unsorted(rev(brts)), "\n")
r <- tryCatch(mc_loglik(brts = brts, pars = c(0.1, 0.5, -0.01, 0.01), sample_size = 1,
                        maxN = 50, max_missing = 500, max_lambda = 50000,
                        lower_bound = c(0, 0, -0.1, -0.1), upper_bound = c(0.5, 2, 0.1, 0.1),
                        xtol_rel = 1e-5, num_threads = 1),
              error = function(e) conditionMessage(e))
cat("  mc_loglik ->", if (is.character(r)) paste("ERROR:", r) else "ran", "\n")

cat("\n== C. does em_cpp accept the same inputs the test passes (unsorted brts, sample_size 1)? ==\n")
mb <- c(1L, 0L, 0L)
p8 <- emphasis:::.expand_pars(c(0.5, -0.01, 0.1, 0.01), mb)
lb <- emphasis:::.expand_pars(c(0, -0.1, 0, -0.1), mb)
ub <- emphasis:::.expand_pars(c(2, 0.1, 0.5, 0.1), mb)
for (srt in c(FALSE, TRUE)) for (ss in c(1L, 20L)) {
  b <- if (srt) sort(as.numeric(brts), decreasing = TRUE) else as.numeric(brts)
  out <- tryCatch(emphasis:::em_cpp(brts = b, init_pars = p8, sample_size = ss, maxN = 500,
                                    max_missing = 500, max_lambda = 50000, lower_bound = lb,
                                    upper_bound = ub, xtol_rel = 1e-5, num_threads = 1,
                                    copy_trees = FALSE, model = mb, link = 0L, rho = 1),
                  error = function(e) conditionMessage(e))
  if (is.character(out)) cat(sprintf("  sorted=%s ss=%d -> ERROR: %s\n", srt, ss, out))
  else cat(sprintf("  sorted=%s ss=%d -> list, fhat in names=%s, fhat=%.3f, n_logf=%d, rejected=%d\n",
                   srt, ss, "fhat" %in% names(out), out$fhat, length(out$logf), out$rejected))
}

cat("\n== D. exponential link (link = 1L), other thing the verifier did not vary ==\n")
b <- sort(as.numeric(brts), decreasing = TRUE)
out <- tryCatch(emphasis:::em_cpp(brts = b, init_pars = emphasis:::.expand_pars(c(log(0.5), log(0.1)), c(0L,0L,0L)),
                                  sample_size = 10, maxN = 500, max_missing = 500, max_lambda = 50000,
                                  lower_bound = emphasis:::.expand_pars(c(-5,-5), c(0L,0L,0L)),
                                  upper_bound = emphasis:::.expand_pars(c(2,2), c(0L,0L,0L)),
                                  xtol_rel = 1e-4, num_threads = 1, copy_trees = FALSE,
                                  model = c(0L,0L,0L), link = 1L, rho = 1),
                error = function(e) conditionMessage(e))
if (is.character(out)) cat("  ERROR:", out, "\n") else
  cat(sprintf("  link=1: fhat in names=%s, fhat=%.3f\n", "fhat" %in% names(out), out$fhat))

cat("\n== E. fix sketch's E-step block, run under testthat 3 times (clock-seeded C++ RNG) ==\n")
ok <- replicate(3, {
  with_reporter("silent", test_that("H96 fix sketch E-step", {
    brts <- sort(as.numeric(ape::branching.times(ape::rphylo(8, 0.5, 0))), decreasing = TRUE)
    pars8 <- emphasis:::.expand_pars(c(0.5, 0.1), c(0L, 0L, 0L))
    aug <- emphasis:::augment_trees(brts, pars8, sample_size = 20, maxN = 500, max_missing = 500,
                                    max_lambda = 50000, num_threads = 1, model = c(0L, 0L, 0L))
    expect_length(aug$logf, 20)
    fh <- emphasis:::.is_fhat(aug$logf, aug$logg)
    expect_true(is.finite(if (is.list(fh)) fh$fhat else fh[1]))
  }))
})
cat("  testthat pass x3:", paste(ok, collapse = ","), "\n")
fh <- emphasis:::.is_fhat(rnorm(5), rnorm(5))
cat("  .is_fhat returns:", class(fh), "of length", length(fh),
    if (is.list(fh)) paste("names:", paste(names(fh), collapse = ",")) else "", "\n")
cat("  is.finite(.is_fhat(...)) as written in fix sketch:",
    tryCatch(paste(is.finite(fh), collapse=","), error = function(e) paste("ERROR:", conditionMessage(e))), "\n")

cat("\n== F. fix sketch's em_cpp rewrite of test-em, 3 replicates ==\n")
brts_em <- c(0.8, 0.6, 0.4, 0.2)
pars8 <- c(0.5, -0.01, 0.01, 0, 0.1, 0, 0, 0); lb8 <- c(0, -0.1, -0.1, 0, 0, 0, 0, 0); ub8 <- c(2, 0.1, 0.1, 0, 0.5, 0, 0, 0)
ok <- replicate(3, with_reporter("silent", test_that("H96 fix sketch em", {
  t0 <- Sys.time()
  result <- emphasis:::em_cpp(brts = brts_em, init_pars = pars8, sample_size = 10, maxN = 100,
                              max_missing = 1000, max_lambda = 500, lower_bound = lb8, upper_bound = ub8,
                              xtol_rel = 1e-3, num_threads = 1, copy_trees = FALSE, model = c(1L, 1L, 0L))
  expect_type(result, "list"); expect_true("fhat" %in% names(result)); expect_length(result$logf, 10)
  cat("   elapsed:", format(as.numeric(Sys.time() - t0, units = "secs"), digits = 3), "s\n")
})))
cat("  testthat pass x3:", paste(ok, collapse = ","), "\n")

## H96 -- `mc_loglik()` is undefined; tests/testthat/test-augment.R:10 and
##        test-em.R:12 call it inside unconditionally skip()-ed blocks.
## Decides: (1) is the symbol bound anywhere in the installed package?
##          (2) what happens when the two test bodies run without skip()?
##          (3) does the em_cpp() rewrite of each body return a list with fhat?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages(library(emphasis))
library(testthat)

cat("== 1. symbol lookup ==\n")
ns <- asNamespace("emphasis")
cat("exists('mc_loglik', envir = ns):     ", exists("mc_loglik", envir = ns, inherits = FALSE), "\n")
cat("'mc_loglik' %in% ls(ns, all=TRUE):   ", "mc_loglik" %in% ls(ns, all.names = TRUE), "\n")
cat("in NAMESPACE exports:                 ", "mc_loglik" %in% getNamespaceExports("emphasis"), "\n")
ga <- getAnywhere("mc_loglik")
cat("getAnywhere() hits (all loaded pkgs): ", length(ga$objs), "\n")
cat("grep -rn mc_loglik R/ src/ NAMESPACE man/ ->",
    length(system2("grep", c("-rln", "mc_loglik",
                             "/Users/pancho/Code/emphasis/R",
                             "/Users/pancho/Code/emphasis/src",
                             "/Users/pancho/Code/emphasis/NAMESPACE",
                             "/Users/pancho/Code/emphasis/man"),
                   stdout = TRUE, stderr = FALSE)), "files\n")
cat("C++ symbols exported via Rcpp:        ",
    paste(grep("_cpp$|^augment_trees$|^eval_logf$", ls(ns, all.names = TRUE), value = TRUE),
          collapse = ", "), "\n")

cat("\n== 2. the two test bodies, skip() removed, run verbatim ==\n")
## test-augment.R:6-22 (skip_on_cran + skip dropped)
r_aug <- tryCatch({
  set.seed(42)
  tree <- ape::rphylo(8, 0.5, 0)
  brts <- ape::branching.times(tree)
  mc_loglik(brts = brts, pars = c(0.1, 0.5, -0.01, 0.01), sample_size = 1,
            maxN = 50, max_missing = 500, max_lambda = 50000,
            lower_bound = c(0, 0, -0.1, -0.1), upper_bound = c(0.5, 2, 0.1, 0.1),
            xtol_rel = 1e-5, num_threads = 1)
}, error = function(e) conditionMessage(e))
cat("test-augment.R:10 ->", if (is.character(r_aug)) paste("ERROR:", r_aug) else "ran", "\n")

## test-em.R:6-24 (skip dropped)
brts_em <- c(0.8, 0.6, 0.4, 0.2)
pars8 <- c(0.5, -0.01, 0.01, 0, 0.1, 0, 0, 0)
lb8 <- c(0, -0.1, -0.1, 0, 0, 0, 0, 0)
ub8 <- c(2, 0.1, 0.1, 0, 0.5, 0, 0, 0)
r_em <- tryCatch({
  mc_loglik(brts = brts_em, pars = pars8, sample_size = 10, maxN = 100,
            max_missing = 1000, max_lambda = 500, lower_bound = lb8,
            upper_bound = ub8, xtol_rel = 1e-3, num_threads = 1,
            model = c(1L, 1L, 0L))
}, error = function(e) conditionMessage(e))
cat("test-em.R:12      ->", if (is.character(r_em)) paste("ERROR:", r_em) else "ran", "\n")

## Same thing through testthat, so the failure class is the one CI would see.
tr <- with_reporter("silent", {
  test_that("H96 unskipped test-em body", {
    expect_type(mc_loglik(brts = brts_em, pars = pars8, sample_size = 10,
                          maxN = 100, max_missing = 1000, max_lambda = 500,
                          lower_bound = lb8, upper_bound = ub8, xtol_rel = 1e-3,
                          num_threads = 1, model = c(1L, 1L, 0L)), "list")
  })
})
cat("testthat outcome:", if (isTRUE(tr)) "PASS" else "FAIL (error)", "\n")

cat("\n== 3. em_cpp() rewrite of each body ==\n")
## test-em.R body: pars8 is already in the 8-slot layout; rename mc_loglik -> em_cpp,
## add copy_trees = FALSE (a required formal of em_cpp).
fit_em <- emphasis:::em_cpp(brts = brts_em, init_pars = pars8, sample_size = 10,
                            maxN = 100, max_missing = 1000, max_lambda = 500,
                            lower_bound = lb8, upper_bound = ub8, xtol_rel = 1e-3,
                            num_threads = 1, copy_trees = FALSE,
                            model = c(1L, 1L, 0L), link = 0L, rho = 1)
cat("em_cpp(test-em inputs): type =", typeof(fit_em), "; names =",
    paste(names(fit_em), collapse = ","), "\n")
cat("  'fhat' %in% names:", "fhat" %in% names(fit_em),
    "; fhat =", format(fit_em$fhat, digits = 6),
    "; length(logf) =", length(fit_em$logf), "\n")

## test-augment.R body: its 4-vector c(0.1, 0.5, -0.01, 0.01) is a pre-8-slot
## layout (beta_0, gamma_0, beta_N, gamma_N?) with no model argument; the current
## API needs a model_bin and the 8-slot layout. Take model = dd (N only) and
## compact order c(beta_0, beta_N, gamma_0, gamma_N) = c(0.5, -0.01, 0.1, 0.01).
set.seed(42)
tree <- ape::rphylo(8, 0.5, 0)
brts_aug <- as.numeric(sort(ape::branching.times(tree), decreasing = TRUE))
mb <- c(1L, 0L, 0L)
p8  <- emphasis:::.expand_pars(c(0.5, -0.01, 0.1, 0.01), mb)
lb  <- emphasis:::.expand_pars(c(0, -0.1, 0, -0.1), mb)
ub  <- emphasis:::.expand_pars(c(2, 0.1, 0.5, 0.1), mb)
fit_aug <- emphasis:::em_cpp(brts = brts_aug, init_pars = p8, sample_size = 20,
                             maxN = 500, max_missing = 500, max_lambda = 50000,
                             lower_bound = lb, upper_bound = ub, xtol_rel = 1e-5,
                             num_threads = 1, copy_trees = FALSE,
                             model = mb, link = 0L, rho = 1)
cat("em_cpp(test-augment inputs, dd): 'fhat' %in% names:", "fhat" %in% names(fit_aug),
    "; fhat =", format(fit_aug$fhat, digits = 6),
    "; estimates =", paste(format(fit_aug$estimates[c(1, 2, 5, 6)], digits = 4), collapse = " "), "\n")

## E-step-only alternative (augment_trees + .is_fhat), same inputs, cr model:
aug <- emphasis:::augment_trees(brts = brts_aug, pars = emphasis:::.expand_pars(c(0.5, 0.1), c(0L,0L,0L)),
                                sample_size = 20, maxN = 500, max_missing = 500,
                                max_lambda = 50000, num_threads = 1,
                                model = c(0L, 0L, 0L), link = 0L, rho = 1)
fh <- emphasis:::.is_fhat(aug$logf, aug$logg)
cat("augment_trees + .is_fhat (cr): fhat =", format(if (is.list(fh)) fh$fhat else fh[1], digits = 6),
    "; trees =", length(aug$logf), "\n")

cat("\n== 4. what the rest of the suite does today ==\n")
res <- as.data.frame(test_dir("/Users/pancho/Code/emphasis/tests/testthat",
                              package = "emphasis", load_package = "none",
                              reporter = "silent", stop_on_failure = FALSE))
cat(sprintf("test_that blocks: %d ; skipped: %d ; failed: %d ; errored: %d\n",
            nrow(res), sum(res$skipped), sum(res$failed), sum(res$error)))
cat("skipped blocks per file:\n")
print(table(res$file[res$skipped]))

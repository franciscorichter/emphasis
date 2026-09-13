.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
library(emphasis); ns <- asNamespace("emphasis"); attach(ns, name = "ns", warn.conflicts = FALSE)
src <- readLines("/Users/pancho/Code/emphasis/tests/testthat/test-bdi-dd.R")
eval(parse(text = src[1:34]))  # trees + dd_bin/cr_bin/dd_pars
dd_pars <- function(l0, m0, K) c(l0, -(l0 - m0) / K, m0, 0)
cat("PRE-FIX behavioural pins\n")
# H11 B1
set.seed(2); e <- .augment_tree_bdi(brts_dd9, c(1.5,-0.12,0.4,0), dd_bin, sample_size=200L, max_missing=30L, link=0L, rho=1)
cat("B1: n trees", length(e$trees), " nonfinite weights", sum(!is.finite(e$weights)), " fhat", e$fhat, "\n")
# all-nonfinite
set.seed(5); e <- .augment_tree_bdi(brts_dd20, dd_pars(0.6,0.1,15), dd_bin, sample_size=20L, max_missing=1e4L, link=0L, rho=1)
cat("all-nonfinite: n trees", length(e$trees), " nonfinite", sum(!is.finite(e$weights)), " fhat", e$fhat, "\n")
# H59
set.seed(3); e2 <- .augment_tree_bdi(brts_cr20, c(0.5,0.4), cr_bin, sample_size=50L, max_missing=2L, link=0L, rho=1)
cat("H59 max_missing=2: n trees", length(e2$trees), " names:", paste(names(e2), collapse=","), "\n")

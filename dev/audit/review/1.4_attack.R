.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib-wave1", .libPaths()))
library(emphasis); ns <- asNamespace("emphasis"); attach(ns, name = "ns", warn.conflicts = FALSE)
src <- readLines("/Users/pancho/Code/emphasis/tests/testthat/test-bdi-dd.R")
eval(parse(text = src[1:34]))
dd_pars <- function(l0, m0, K) c(l0, -(l0 - m0) / K, m0, 0)
show <- function(tag, e) cat(sprintf("%-42s fhat=%9.4f acc=%6.3f nv=%4d nnf=%3d att=%4d rej=%4d mm=%4d ntree=%4d ESS=%7.1f\n",
  tag, e$fhat, e$acc, e$n_valid, e$n_nonfinite, e$n_attempts, e$n_rejected, e$n_rejected_max_missing, length(e$trees), .ess_from_lw(e$weights)))
W <- function(expr) withCallingHandlers(expr, warning = function(w) { cat("   [warning] ", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
thin <- function(brts, pars, mb, N, link = 0L, mm = 1e4) { a <- .augment_tree_internal(brts, pars, mb, sample_size = N, max_missing = mm, link = link, num_threads = 1L); c(fhat = a$fhat, n = length(a$trees), zero = a$rejected_zero_weights, ov = a$rejected_overruns) }

cat("\n== 1. DD link 0 vs thinning vs DDD (20-tip) ==\n")
set.seed(1)
for (th in list(c(0.8,0.3,20), c(1.0,0.8,20), c(0.6,0.1,50), c(0.8,0.3,1e4))) {
  p <- dd_pars(th[1],th[2],th[3]); ref <- DDD::dd_loglik(c(th), c(300,1,0,1,0,2), brts_dd20, 0)
  e <- W(.augment_tree_bdi(brts_dd20, p, dd_bin, 500L, 1e4L, 0L, 1)); t <- thin(brts_dd20, p, dd_bin, 500L)
  show(sprintf("dd0 (%.1f,%.1f,K=%g) ref=%.3f thin=%.3f", th[1],th[2],th[3], ref, t["fhat"]), e)
  cat(sprintf("   identity check: %.3e\n", e$fhat - (log(sum(exp(e$weights - max(e$weights)))/(e$n_valid + e$n_rejected)) + max(e$weights))))
}
cat("\n== 2. DD link 1 (exponential) vs thinning ==\n")
set.seed(2)
for (p in list(c(log(0.8), -0.02, log(0.3), 0), c(log(1.0), -0.03, log(0.6), 0))) {
  e <- W(.augment_tree_bdi(brts_dd20, p, dd_bin, 500L, 1e4L, 1L, 1)); t <- thin(brts_dd20, p, dd_bin, 500L, link = 1L)
  show(sprintf("dd1 (%s) thin=%.3f", paste(round(p,2), collapse=","), t["fhat"]), e)
}
cat("\n== 3. CR boundaries (1.3 interaction): mu=lambda, mu>lambda, mu=0, tiny rates ==\n")
set.seed(3)
for (lm in list(c(0.5,0.5), c(0.3,0.5), c(0.3,0), c(1e-4,5e-5), c(5, 4.9))) {
  ref <- tryCatch(DDD::bd_loglik(c(lm,0,0), c(0,0,1,0,2), brts_dd20, 0), error = function(e) NA)
  e <- tryCatch(W(.augment_tree_bdi(brts_dd20, lm, cr_bin, 30L, 1e4L, 0L, 1)), error = function(err) { cat("  ERROR", conditionMessage(err), "\n"); NULL })
  if (!is.null(e)) { show(sprintf("cr (%g,%g) ref=%.4f sd_lw=%.1e", lm[1], lm[2], ref, sd(e$weights)), e) }
}
cat("\n== 4. DD boundaries: mu=lambda, mu>lambda, K<n ==\n")
set.seed(4)
for (th in list(c(0.5,0.5,30), c(0.3,0.5,30), c(0.8,0.3,10), c(0.8,0.3,19))) {
  p <- dd_pars(th[1],th[2],th[3])
  e <- tryCatch(W(.augment_tree_bdi(brts_dd20, p, dd_bin, 100L, 1e4L, 0L, 1)), error = function(err) { cat("  ERROR", conditionMessage(err), "\n"); NULL })
  if (!is.null(e)) show(sprintf("dd (%.1f,%.1f,K=%g)", th[1],th[2],th[3]), e)
}
cat("\n== 5. Tiny / large trees, sample_size 1, max_missing 0 ==\n")
set.seed(5)
e <- W(.augment_tree_bdi(c(5), c(0.5,0.2), cr_bin, 10L, 1e4L, 0L, 1)); show("cr 2-tip N=10", e)
e <- W(.augment_tree_bdi(c(5), dd_pars(0.5,0.2,30), dd_bin, 10L, 1e4L, 0L, 1)); show("dd 2-tip N=10", e)
e <- W(.augment_tree_bdi(c(5,2), dd_pars(0.5,0.2,30), dd_bin, 10L, 1e4L, 0L, 1)); show("dd 3-tip N=10", e)
e <- W(.augment_tree_bdi(brts_dd20, c(0.8,0.3), cr_bin, 1L, 1e4L, 0L, 1)); show("cr N=1", e)
e <- W(.augment_tree_bdi(brts_dd20, dd_pars(0.8,0.3,20), dd_bin, 1L, 1e4L, 0L, 1)); show("dd N=1", e)
e <- W(.augment_tree_bdi(brts_dd20, c(0.8,0.3), cr_bin, 5L, 0L, 0L, 1)); show("cr N=5 max_missing=0", e)
e <- W(.augment_tree_bdi(brts_dd20, c(0.8,0.0), cr_bin, 5L, 0L, 0L, 1)); show("cr mu=0 N=5 max_missing=0", e)
e <- W(.augment_tree_bdi(brts_dd20, dd_pars(0.8,0.3,20), dd_bin, 5L, 0L, 0L, 1)); show("dd N=5 max_missing=0", e)
set.seed(6); big <- sort(ape::branching.times(ape::rphylo(150, 0.6, 0.2)), decreasing = TRUE); big <- as.numeric(big)
t0 <- proc.time()[3]
e <- W(.augment_tree_bdi(big, c(0.6,0.2), cr_bin, 20L, 1e4L, 0L, 1)); show(sprintf("cr 150-tip N=20 ref=%.3f", DDD::bd_loglik(c(0.6,0.2,0,0), c(0,0,1,0,2), big, 0)), e)
e <- W(.augment_tree_bdi(big, dd_pars(0.6,0.2,200), dd_bin, 50L, 1e4L, 0L, 1)); show(sprintf("dd 150-tip N=50 ref=%.3f", DDD::dd_loglik(c(0.6,0.2,200), c(300,1,0,1,0,2), big, 0)), e)
cat("   time", proc.time()[3] - t0, "\n")
cat("\n== 6. DD low acceptance / budget exhaustion ==\n")
set.seed(7)
for (th in list(c(1.5,1.4,20), c(2,1.9,20), c(1.0,0.95,20))) {
  p <- dd_pars(th[1],th[2],th[3])
  e <- W(.augment_tree_bdi(brts_dd20, p, dd_bin, 20L, 1e4L, 0L, 1)); show(sprintf("dd (%.2f,%.2f,K=%g) N=20", th[1],th[2],th[3]), e)
  cat(sprintf("   returned weights finite: %d/%d ; lengths logf/logg/trees: %d/%d/%d\n", sum(is.finite(e$weights)), length(e$weights), length(e$logf), length(e$logg), length(e$trees)))
}
cat("\n== 7. n_valid=0 cases ==\n")
e <- W(.augment_tree_bdi(brts_cr20, c(0.5,0.4), cr_bin, 3L, 0L, 0L, 1)); show("cr all overflow", e); str(e[c("fhat","acc","trees","weights")])
cat("\n== 8. simulate_tree(method='bdi') ==\n")
set.seed(8)
s <- W(simulate_tree(brts_dd20, pars = c(0.8,0.3), model = cr_bin, n_trees = 3L, method = "bdi")); cat("cr sim class:", class(s), " len:", length(s), "\n")
s <- W(tryCatch(simulate_tree(brts_dd20, pars = dd_pars(0.8,0.3,20), model = dd_bin, n_trees = 3L, method = "bdi"), error = function(e) conditionMessage(e))); cat("dd sim class:", class(s), " len:", length(s), "\n")
s <- W(tryCatch(simulate_tree(brts_dd20, pars = dd_pars(2,1.9,20), model = dd_bin, n_trees = 1L, method = "bdi"), error = function(e) conditionMessage(e))); cat("dd low-acc sim n_trees=1:", class(s), "\n")

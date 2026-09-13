args <- commandArgs(trailingOnly = TRUE); lib <- args[1]
.libPaths(c(lib, .libPaths()))
library(emphasis); ns <- asNamespace("emphasis"); attach(ns, name = "ns", warn.conflicts = FALSE)
cat("BUILD:", basename(lib), "\n")
src <- readLines("/Users/pancho/Code/emphasis/tests/testthat/test-bdi-dd.R"); eval(parse(text = src[1:34]))
dd_pars <- function(l0, m0, K) c(l0, -(l0 - m0) / K, m0, 0)
W <- function(expr) withCallingHandlers(expr, warning = function(w) { cat("   [warning] ", substr(conditionMessage(w),1,110), "\n"); invokeRestart("muffleWarning") })
post <- grepl("wave1", lib)

cat("\n== thinning cross-check (names of augment_trees output) ==\n")
set.seed(1); a <- .augment_tree_internal(brts_dd20, dd_pars(0.8,0.3,20), dd_bin, sample_size = 500L, num_threads = 1L)
cat(names(a), "\n"); lw <- a$logf - a$logg; S <- length(a$trees) + .n0(a$rejected_zero_weights); tf <- log(sum(exp(lw - max(lw)))/S) + max(lw)
cat("thin fhat (S=trees+zero_w):", tf, " zero_w:", .n0(a$rejected_zero_weights), " overruns:", .n0(a$rejected_overruns), " DDD ref:", DDD::dd_loglik(c(0.8,0.3,20), c(300,1,0,1,0,2), brts_dd20, 0), "\n")
if (post) { set.seed(1); e <- .augment_tree_bdi(brts_dd20, dd_pars(0.8,0.3,20), dd_bin, 500L, 1e4L, 0L, 1); cat("bdi fhat:", e$fhat, "\n") }

cat("\n== DD max_missing binding (post only) ==\n")
if (post) { set.seed(2); for (mm in c(1e4, 8, 5)) { e <- W(.augment_tree_bdi(brts_dd20, dd_pars(1.0,0.8,20), dd_bin, 300L, mm, 0L, 1)); cat(sprintf("  mm=%5g fhat=%.3f acc=%.3f nv=%d rej=%d mmrej=%d\n", mm, e$fhat, e$acc, e$n_valid, e$n_rejected, e$n_rejected_max_missing)) } }

cat("\n== rho < 1 passthrough (post only, H2 out of scope: just no crash) ==\n")
if (post) { set.seed(3); e <- W(.augment_tree_bdi(brts_dd20, c(0.8,0.3), cr_bin, 10L, 1e4L, 0L, 0.7)); cat("  rho=0.7 cr fhat", e$fhat, "acc", e$acc, "\n") }

cat("\n== simulate_tree(method='bdi') with phylo ==\n")
set.seed(4); phy <- ape::rphylo(20, 0.5, 0.1)
s <- W(tryCatch(simulate_tree(phy, pars = c(0.8,0.3), model = cr_bin, n_trees = 3L, method = "bdi"), error = function(e) paste("ERR", conditionMessage(e)))); cat("  cr sim:", class(s)[1], length(s), "\n")
s <- W(tryCatch(simulate_tree(phy, pars = dd_pars(0.8,0.3,20), model = dd_bin, n_trees = 3L, method = "bdi"), error = function(e) paste("ERR", conditionMessage(e)))); cat("  dd sim:", class(s)[1], length(s), "\n")
s <- W(tryCatch(simulate_tree(phy, pars = c(0.5,0.4), model = cr_bin, n_trees = 3L, method = "bdi", max_missing = 2), error = function(e) paste("ERR", conditionMessage(e)))); cat("  cr sim max_missing=2:", class(s)[1], length(s), "\n")

cat("\n== estimate_rates dd via BDI: reported quantities ==\n")
set.seed(5); t0 <- proc.time()[3]
fit <- W(tryCatch(estimate_rates(brts_dd20, model = dd_bin, method = "mcem", init_pars = dd_pars(0.8,0.3,25),
                    control = list(sampling = "bdi", num_trees = 200, max_iter = 12, num_threads = 1, verbose = FALSE,
                                   lower_bound = c(0, -0.2, 0, 0), upper_bound = c(3, 0, 3, 0))),
       error = function(e) paste("ERR", conditionMessage(e))))
cat("  time", round(proc.time()[3]-t0,1), "\n")
if (is.character(fit)) cat(fit, "\n") else {
  cat("  pars:", round(fit$pars,4), " loglik:", fit$loglik, " loglik_var:", fit$loglik_var, " iter:", fit$iterations, " stop:", fit$stop_reason, "\n")
  m <- fit$mcem; if (!is.null(m)) { cat("  mcem cols:", names(m), "\n"); print(tail(m[, intersect(names(m), c("fhat","delta_max","rejected","n_nonfinite","num_trees","ESS","m_step"))], 4)) }
  fi <- fit$final_IS; if (!is.null(fi)) cat("  final_IS n_rejected:", fi$n_rejected, " zero_w:", fi$rejected_zero_weights, " ESS:", fi$ESS, " n:", length(fi$logf), " fhat:", fi$fhat, "\n")
  p <- fit$pars; K <- -(p[1]-p[3])/p[2]
  cat("  DDD dd_loglik at estimate (btorph=1):", tryCatch(DDD::dd_loglik(c(p[1],p[3],K), c(300,1,0,1,0,2), brts_dd20, 0), error=function(e) NA), "\n")
  cat("  other fields:", setdiff(names(fit), c("mcem","final_IS")), "\n")
  if (!is.null(fit$diagnostics)) str(fit$diagnostics, max.level = 1)
}
cat("\n== estimate_rates at the H11 B1 configuration (9-tip, non-finite draws) ==\n")
set.seed(6)
fit2 <- W(tryCatch(estimate_rates(brts_dd9, model = dd_bin, method = "mcem", init_pars = c(1.5,-0.12,0.4,0),
                    control = list(sampling = "bdi", num_trees = 200, max_iter = 6, num_threads = 1, verbose = FALSE,
                                   lower_bound = c(0, -0.5, 0, 0), upper_bound = c(3, 0, 3, 0))),
       error = function(e) paste("ERR", conditionMessage(e))))
if (is.character(fit2)) cat(fit2, "\n") else {
  cat("  pars:", round(fit2$pars,4), " loglik:", fit2$loglik, " iter:", fit2$iterations, " stop:", fit2$stop_reason, "\n")
  m <- fit2$mcem; if (!is.null(m)) print(m[, intersect(names(m), c("fhat","delta_max","rejected","n_nonfinite","num_trees","ESS","m_step"))])
  fi <- fit2$final_IS; if (!is.null(fi)) cat("  final_IS n_rejected:", fi$n_rejected, " zero_w:", fi$rejected_zero_weights, " n:", length(fi$logf), "\n")
}

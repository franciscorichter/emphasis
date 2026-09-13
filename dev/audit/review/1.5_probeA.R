args <- commandArgs(trailingOnly = TRUE); lib <- args[1]
.libPaths(c(lib, .libPaths())); suppressPackageStartupMessages({library(emphasis); library(testthat)})
ns <- asNamespace("emphasis"); cat("lib:", find.package("emphasis"), "\n")
brts12 <- c(6.385824, 2.063997, 1.19255, 0.923743, 0.884126, 0.822976, 0.725585, 0.718801, 0.539214, 0.077287, 0.068123)
brts_dd <- c(6, 4.848493, 4.401821, 3.108164, 3.073914, 2.835023, 1.838828, 0.50463)
cr8 <- function(l, m) c(l,0,0,0,m,0,0,0)
run <- function(pars, brts = brts12, model = c(0L,0L,0L), link = 0L, lb = cr8(0,0), ub = cr8(3,3), max_iter = 4L, N = 30L, tol = 1e-2, threads = 1L, max_time = NULL, max_missing = 1e4)
  ns$.mcem_bdi(brts, pars = pars, sample_size = N, max_missing = max_missing, lower_bound = lb, upper_bound = ub, max_iter = max_iter, xtol = 1e-3, tol = tol, patience = 3L, num_threads = threads, model = model, link = link, rho = 1, max_time = max_time)
sm <- function(r) cat(sprintf("  stop=%s iters=%d n_failed=%d loglik=%.4f var=%.3g nrow(mcem)=%s pars=%s\n  n_nonfinite col=%s rejected col=%s ESS=%s\n", r$stop_reason, r$iterations, r$n_failed, r$loglik, r$loglik_var, if (is.null(r$mcem)) "NULL" else nrow(r$mcem), paste(signif(r$pars[c(1,2,5,6)],4), collapse=","), if (is.null(r$mcem)) "-" else paste(r$mcem$n_nonfinite, collapse=","), if (is.null(r$mcem)) "-" else paste(r$mcem$rejected, collapse=","), if (is.null(r$final_IS)) "-" else signif(r$final_IS$ESS,4)))

cat("\n[1] n_nonfinite recording under DD/H11 config (1.4 strips draws before .mcem_bdi sees them)\n")
set.seed(11)
e <- ns$.augment_tree_bdi(brts_dd, c(1.5,-0.12,0.4,0), model_bin = c(1L,0L,0L), sample_size = 200L, max_missing = 30L, link = 0L, rho = 1)
cat(sprintf("  .augment_tree_bdi: n_valid=%d n_nonfinite=%d n_rejected=%d n_rej_mm=%d length(trees)=%d any nonfinite logf=%s\n", e$n_valid, e$n_nonfinite, e$n_rejected, e$n_rejected_max_missing, length(e$trees), any(!is.finite(e$logf))))
set.seed(11)
r <- run(c(1.5,-0.12,0.4,0), brts = brts_dd, model = c(1L,0L,0L), lb = c(0.01,-1,0.001,0), ub = c(5,0,2,0), max_iter = 3L, N = 200L, max_missing = 30L)
sm(r); cat("  num_trees col:", paste(r$mcem$num_trees, collapse=","), " final_IS rejected_zero_weights:", r$final_IS$rejected_zero_weights, "\n")

cat("\n[2] final E-step failure -> does estimate_rates leak the lagged fhat as loglik?\n")
real_aug <- ns$.augment_tree_bdi; ncall <- 0L
with_mocked_bindings(.augment_tree_bdi = function(tree, pars, ...) { ncall <<- ncall + 1L; if (ncall == 3L) stop("final E-step fails"); real_aug(tree, pars, ...) }, .package = "emphasis", {
  set.seed(20)
  fit <- estimate_rates(brts12, model = "cr", method = "mcem", init_pars = c(1.2, 0.9), control = list(lower_bound = c(0,0), upper_bound = c(3,3), sampling = "bdi", sample_size = 30L, num_threads = 1L, max_iter = 2L, tol = 1e-2))
  cat(sprintf("  calls=%d details$loglik=%s fit$loglik=%s details$loglik_var=%s fit$loglik_var=%s AIC=%s trace fhat=%s iterations(top)=%s stop(top)=%s\n", ncall, r1 <- fit$details$loglik, fit$loglik, fit$details$loglik_var, fit$loglik_var, fit$AIC, paste(round(fit$details$mcem$fhat,3), collapse=","), fit$iterations, fit$stop_reason))
})

cat("\n[3] sample_size = 1\n"); set.seed(5); r <- run(cr8(0.5,0.2), N = 1L, max_iter = 3L); sm(r)
cat("\n[3b] sample_size = 1 via estimate_rates\n"); set.seed(5)
fit <- tryCatch(estimate_rates(brts12, model = "cr", method = "mcem", init_pars = c(0.5,0.2), control = list(lower_bound = c(0,0), upper_bound = c(3,3), sampling = "bdi", sample_size = 1L, num_threads = 1L, max_iter = 3L)), error = function(e) e)
if (inherits(fit, "error")) cat("  ERROR:", conditionMessage(fit), "\n") else { cat(sprintf("  loglik=%s var=%s iters=%s\n", fit$loglik, fit$loglik_var, fit$details$iterations)); print(fit) }

cat("\n[4] max_iter = 0\n"); set.seed(6); r <- run(cr8(0.5,0.2), max_iter = 0L); sm(r)
cat("  .run_mcem-style loglik from trace:", tail(r$mcem$fhat, 1), "\n")

cat("\n[5] exponential link (1), cr and dd\n")
set.seed(7); r <- run(c(log(0.5),0,0,0,log(0.2),0,0,0), link = 1L, lb = cr8(-5,-5), ub = cr8(3,3), max_iter = 4L); sm(r)
set.seed(7); r <- run(c(log(1.5),-0.05,log(0.4),0), brts = brts_dd, model = c(1L,0L,0L), link = 1L, lb = c(-5,-1,-5,0), ub = c(3,0,3,0), max_iter = 3L, N = 50L, max_missing = 30L); sm(r)

cat("\n[6] boundary starts: lambda = mu, mu > lambda, rates near 0\n")
set.seed(8); r <- run(cr8(0.4,0.4), max_iter = 4L); sm(r)
set.seed(8); r <- run(cr8(0.3,0.5), max_iter = 4L); sm(r)
set.seed(8); r <- run(cr8(1e-4,1e-6), max_iter = 3L); sm(r)

cat("\n[7] num_threads = 4 vs 1, same seed\n")
set.seed(9); r1 <- run(cr8(0.5,0.2), max_iter = 3L, threads = 1L); set.seed(9); r4 <- run(cr8(0.5,0.2), max_iter = 3L, threads = 4L)
cat("  pars1:", signif(r1$pars[c(1,5)],5), " pars4:", signif(r4$pars[c(1,5)],5), " loglik:", r1$loglik, r4$loglik, "\n")

cat("\n[8] max_time = 0\n"); set.seed(10); r <- run(cr8(0.5,0.2), max_iter = 10L, max_time = 0); sm(r)

cat("\n[9] all bounds fixed (lb == ub): M-step returns start -> never 'converged'?\n")
set.seed(12); r <- run(cr8(0.5,0.2), lb = cr8(0.5,0.2), ub = cr8(0.5,0.2), max_iter = 6L); sm(r); cat("  m_moved:", paste(r$mcem$m_moved, collapse=","), "\n")

cat("\n[10] livelock: E-step fails whenever lambda > 0.55 (alternating success/failure)\n")
ncall <- 0L
with_mocked_bindings(.augment_tree_bdi = function(tree, pars, ...) { ncall <<- ncall + 1L; if (pars[1] > 0.55) stop("region failure"); real_aug(tree, pars, ...) }, m_cpp = function(e_step, init_pars, ...) { est <- init_pars; est[1] <- 0.6; list(estimates = est, nlopt = 4L, time = 0) }, .package = "emphasis", {
  set.seed(13); r <- run(cr8(0.5,0.2), max_iter = 30L, tol = 0); sm(r); cat("  calls:", ncall, "\n")
})

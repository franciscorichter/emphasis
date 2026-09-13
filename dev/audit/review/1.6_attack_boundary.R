source("/Users/pancho/Code/emphasis/dev/audit/review/1.6_common.R")
set.seed(11)
d <- drv(brts22, c(1,1), c(0,0), c(4,4)); invariants(d, ex(c(1,1)), "mu==lambda init")
d <- drv(brts22, c(0.5,1.5), c(0,0), c(4,4)); invariants(d, ex(c(0.5,1.5)), "mu>lambda init")
d <- drv(brts22, c(1,0), c(0,0), c(4,4)); invariants(d, ex(c(1,0)), "mu=0 init (lb 0)")
d <- drv(brts22, c(1,1e-3), c(0,0), c(4,4), max_iter = 60L, sample_size = 200L); invariants(d, ex(c(1,1e-3)), "mu near 0, N=200, max_iter 60")
cat("   delta_max tail:", paste(signif(tail(d$mcem$delta_max[!d$mcem$final_estep], 8),2), collapse=" "), " mu tail:", paste(signif(tail(d$mcem$par5[!d$mcem$final_estep], 5),2), collapse=" "), "\n")
d <- drv(brts22, c(1,0.3), c(0,0), c(4,4), sample_size = 1L); invariants(d, ex(c(1,0.3)), "sample_size 1")
d <- drv(brts22, c(1,0.3), c(0,0), c(4,4), sample_size = 20L, maxN = 20L); invariants(d, ex(c(1,0.3)), "maxN == sample_size")
cat("   maxN col:", paste(d$mcem$maxN, collapse=" "), "\n")
d <- drv(c(1.0, 0.5), c(1,0.3), c(0,0), c(4,4)); invariants(d, ex(c(1,0.3)), "3-tip tree")
d <- tryCatch(drv(c(1.0), c(1,0.3), c(0,0), c(4,4)), error = function(e) {cat("2-tip tree ERROR:", conditionMessage(e), "\n"); NULL}); if (!is.null(d)) invariants(d, ex(c(1,0.3)), "2-tip tree")
set.seed(3); big <- ape::rphylo(150, 1, 0.3); bb <- sort(as.numeric(ape::branching.times(big)), decreasing = TRUE)
t0 <- proc.time()[3]; d <- drv(bb, c(1,0.3), c(0,0), c(4,4), sample_size = 30L, max_iter = 3L); invariants(d, ex(c(1,0.3)), "150-tip tree"); cat(sprintf("   %.1fs\n", proc.time()[3]-t0))
# tol Inf, patience 1 -> converged after one iteration
d <- drv(brts22, c(1,0.3), c(0,0), c(4,4), tol = Inf, patience = 1L); invariants(d, ex(c(1,0.3)), "tol Inf patience 1")
# num_threads 2
d <- drv(brts22, c(1,0.3), c(0,0), c(4,4), num_threads = 2L, max_iter = 4L); invariants(d, ex(c(1,0.3)), "num_threads 2")
# all-fail with max_iter 5 (< 8): what stop_reason?
boom <- function(p) stop("x")
d <- suppressWarnings(drv(brts22, c(1,0.3), c(0,0), c(4,4), conditional = boom, max_iter = 5L))
cat(sprintf("all-fail max_iter=5: stop=%s iter=%d n_failed=%d final=%s maxN=%d pars==init %s mcem NULL %s\n", d$stop_reason, d$iterations, d$n_failed, d$final_estep, d$maxN, isTRUE(all.equal(d$pars, ex(c(1,0.3)))), is.null(d$mcem)))
fit <- suppressWarnings(tryCatch(estimate_rates(brts22, model="cr", method="mcem", init_pars=c(1,0.3), control=list(sampling="dynamic_fresh", num_trees=10L, max_iter=5L, lower_bound=c(0,0), upper_bound=c(4,4), num_threads=1L), conditional = boom), error=function(e) e))
if (inherits(fit,"error")) cat("estimate_rates all-fail max_iter 5 ERROR:", conditionMessage(fit), "\n") else cat(sprintf("  via estimate_rates: loglik=%s stop=%s iterations=%s AIC=%s\n", fit$loglik, fit$stop_reason, fit$iterations, fit$AIC))
# user maxN > cap with one injected failure: not shrunk
flag <- new.env(); flag$left <- 1L
cond1 <- function(p) { if (flag$left > 0L) { flag$left <- flag$left - 1L; stop("inj") }; 0 }
d <- drv(brts22, c(1,0.3), c(0,0), c(4,4), conditional = cond1, maxN = 60000L, max_iter = 3L)
cat(sprintf("maxN 60000 + 1 failure: maxN col %s final maxN %d n_failed %d\n", paste(d$mcem$maxN, collapse=" "), d$maxN, d$n_failed))
# 7 failures then success: pars = new estimate from perturbed point; prev not mixture
flag$left <- 7L
d <- drv(brts22, c(1,0.3), c(0,0), c(4,4), conditional = cond1, max_iter = 12L, tol = 1e6, patience = 1L)
m <- d$mcem[!d$mcem$final_estep, ]
cent <- ex(c(2,2)); pert <- ex(c(1,0.3)); for (k in 1:6) pert <- pmin(pmax(0.8*pert + 0.2*cent, ex(c(0,0))), ex(c(4,4)))
cat(sprintf("7 failures then success: n_failed=%d iter=%d stop=%s maxN=%d; delta row1 measured vs 6x-perturbed point? %s (%.4f vs %.4f)\n", d$n_failed, d$iterations, d$stop_reason, d$maxN,
   isTRUE(all.equal(m$delta_max[1], max(abs(as.numeric(m[1,paste0("par",1:8)]) - pert)/pmax(abs(pert),1e-2)))), m$delta_max[1], max(abs(as.numeric(m[1,paste0("par",1:8)]) - pert)/pmax(abs(pert),1e-2))))
# max_time honoured while failing
flag$left <- 100L; t0 <- proc.time()[3]
d <- suppressWarnings(drv(brts22, c(1,0.3), c(0,0), c(4,4), conditional = cond1, max_iter = 50L, max_time = 0))
cat(sprintf("always-fail with max_time=0: stop=%s n_failed=%d iter=%d %.2fs\n", d$stop_reason, d$n_failed, d$iterations, proc.time()[3]-t0))
# final E-step failure path: conditional fails only on the (K+1)-th call
flag2 <- new.env(); flag2$n <- 0L
condK <- function(p) { flag2$n <- flag2$n + 1L; if (flag2$n == 3L) stop("final fails"); 0 }
d <- drv(brts22, c(1,0.3), c(0,0), c(4,4), conditional = condK, max_iter = 2L)
cat(sprintf("final E-step fails: final_estep=%s nrow=%d iter=%d n_failed=%d stop=%s loglik_var=%s pars==row2? %s\n", d$final_estep, nrow(d$mcem), d$iterations, d$n_failed, d$stop_reason, d$loglik_var,
    isTRUE(all.equal(as.numeric(d$mcem[2,paste0("par",1:8)]), as.numeric(d$pars)))))

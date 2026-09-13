lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, .libPaths())); library(emphasis)
ns <- asNamespace("emphasis"); A <- ns$.augment_tree_bdi
brts20 <- c(6, 4.27428208320388, 4.1258975168475, 3.72207978339226,
  3.64993282835291, 3.17502914508657, 2.81163149722479, 2.48700670286076,
  2.38274529774947, 1.94421117175661, 1.86664932095232, 1.81164887493753,
  1.38503498765192, 1.06976980095552, 0.460758774034712, 0.417936203054209,
  0.363178692045526, 0.253829586012209, 0.0517782331486725)
cr <- c(0L,0L,0L); dd <- c(1L,0L,0L)
id <- function(e) {  # self-consistency of the documented fhat formula
  if (length(e$weights) == 0L) return(NA_real_)
  m <- max(e$weights); e$fhat - (log(sum(exp(e$weights - m))/e$n_valid) + m + log(e$acc))
}
cat("=== A. CR boundary: mu vs lambda (exactness must survive the log(acc) factor) ===\n")
for (p in list(c(0.5,0.1), c(0.5,0.5), c(0.5,0.5-1e-9), c(0.3,0.5), c(0.1,0.9),
               c(1e-6,1e-7), c(5,4.9))) {
  set.seed(9)
  e <- tryCatch(A(brts20, p, cr, sample_size = 30L, max_missing = 1e4L, link = 0L, rho = 1),
                error = function(x) x)
  if (inherits(e, "error")) { cat(sprintf("  lam=%g mu=%g  ERROR %s\n", p[1],p[2],conditionMessage(e))); next }
  ref <- DDD::bd_loglik(pars1 = c(p,0,0), pars2 = c(0,0,1,0,2), brts = brts20, missnumspec = 0)
  cat(sprintf("  lam=%-8g mu=%-8g nv=%3d rej=%3d mm=%3d acc=%.4f sd(w)=%.2e fhat-ref=%+.3e chk=%.1e\n",
      p[1],p[2],e$n_valid,e$n_rejected,e$n_rejected_max_missing,e$acc,
      stats::sd(e$weights), e$fhat-ref, id(e)))
}
cat("\n=== B. dd, exponential link 1 (the other supported link) ===\n")
for (p in list(c(log(0.8),-0.02,log(0.3),0), c(log(0.6),-0.05,log(0.1),0))) {
  set.seed(9)
  e <- tryCatch(A(brts20, p, dd, sample_size = 200L, max_missing = 1e4L, link = 1L, rho = 1),
                error = function(x) x)
  if (inherits(e,"error")) { cat("  ERROR", conditionMessage(e), "\n"); next }
  cat(sprintf("  nv=%d nf=%d rej=%d acc=%.4f fhat=%.4f chk=%.1e\n",
      e$n_valid,e$n_nonfinite,e$n_rejected,e$acc,e$fhat,id(e)))
}
cat("\n=== C. sample_size = 1 (acc from a single completed draw) ===\n")
set.seed(4)
for (r in 1:6) {
  e <- A(brts20, c(0.8,-(0.8-0.3)/20,0.3,0), dd, sample_size = 1L, max_missing = 1e4L, link = 0L, rho = 1)
  cat(sprintf("  rep%d nv=%d rej=%d acc=%.3f fhat=%.3f\n", r,e$n_valid,e$n_rejected,e$acc,e$fhat))
}
cat("\n=== D. n_valid = 0 (every attempt over max_missing): acc / fhat / trees ===\n")
set.seed(4)
e <- withCallingHandlers(
  A(brts20, c(0.5,0.4), cr, sample_size = 5L, max_missing = 0L, link = 0L, rho = 1),
  warning = function(w) { cat("  warning:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
cat(sprintf("  nv=%d att=%d mm=%d acc=%s fhat=%s ntrees=%d nlogf=%d\n",
    e$n_valid, e$n_attempts, e$n_rejected_max_missing, format(e$acc), format(e$fhat),
    length(e$trees), length(e$logf)))
cat("\n=== E. budget exhausted under DD (survivors): is acc still the right factor? ===\n")
set.seed(4)
e <- withCallingHandlers(
  A(brts20, c(0.6,-0.5/15,0.1,0), dd, sample_size = 40L, max_missing = 1e4L, link = 0L, rho = 1),
  warning = function(w) { cat("  warning:", substr(conditionMessage(w),1,110), "\n"); invokeRestart("muffleWarning") })
cat(sprintf("  nv=%d att=%d rej=%d nf=%d acc=%s fhat=%s\n", e$n_valid,e$n_attempts,
    e$n_rejected,e$n_nonfinite,format(e$acc),format(e$fhat)))

source("/Users/pancho/Code/emphasis/dev/audit/review/1.6_common.R")
set.seed(1)
# which model names resolve?
for (nm in c("cr","dd","d","pd","ddd","n","m","md","nd","dm")) {
  r <- tryCatch(emphasis:::.resolve_model(nm), error = function(e) NA)
  if (!all(is.na(r))) cat("model", nm, "->", paste(r, collapse=""), "\n")
}
cases <- list(
  list(model="cr", link="linear",      init=c(1,0.3),          lb=c(0,0),         ub=c(4,4)),
  list(model="cr", link="exponential", init=c(0,-1.2),         lb=c(-3,-5),       ub=c(2,2)),
  list(model="cr", link="gaussian",    init=c(1,0.3),          lb=c(0,0),         ub=c(4,4)),
  list(model="dd", link="linear",      init=c(1.5,-0.02,0.3,0),lb=c(0,-0.2,0,0),  ub=c(4,0,4,0)),
  list(model="dd", link="exponential", init=c(0.3,-0.02,-1,0), lb=c(-3,-0.5,-5,0),ub=c(2,0.5,2,0))
)
for (cs in cases) {
  t0 <- proc.time()[3]
  fit <- tryCatch(estimate_rates(brts22, model = cs$model, link = cs$link, method = "mcem", init_pars = cs$init,
     control = list(sampling = "dynamic_fresh", num_trees = 30L, max_iter = 6L, tol = 1e-2, patience = 3L,
                    lower_bound = cs$lb, upper_bound = cs$ub, num_threads = 1L, max_time = 100)), error = function(e) e)
  if (inherits(fit, "error")) { cat("ERROR", cs$model, cs$link, conditionMessage(fit), "\n"); next }
  d <- fit$details; mb <- emphasis:::.resolve_model(cs$model); attr(d, "mb") <- mb
  invariants(d, ex(cs$init, mb), paste(cs$model, cs$link))
  cat(sprintf("   fit$loglik=%.4f last fhat=%.4f top-level stop=%s iterations=%s  %.1fs\n", fit$loglik, tail(d$mcem$fhat,1), fit$stop_reason, fit$iterations, proc.time()[3]-t0))
}
# D-dependent model (thinning is the only sampler): try model names with D
for (nm in c("d","pd","md")) {
  mb <- tryCatch(emphasis:::.resolve_model(nm), error = function(e) NULL); if (is.null(mb)) next
  np <- 2 + 2*sum(mb); init <- numeric(np); init[1] <- 1; init[np/2+1] <- 0.3
  lb <- rep(-0.5, np); ub <- rep(4, np); lb[1] <- 0; lb[np/2+1] <- 0
  fit <- tryCatch(estimate_rates(brts22, model = nm, link = "linear", method = "mcem", init_pars = init,
     control = list(sampling = "dynamic_fresh", num_trees = 30L, max_iter = 4L, lower_bound = lb, upper_bound = ub, num_threads = 1L, max_time = 100)), error = function(e) e)
  if (inherits(fit, "error")) { cat("ERROR model", nm, conditionMessage(fit), "\n"); next }
  d <- fit$details; attr(d, "mb") <- mb; invariants(d, ex(init, mb), paste("D-model", nm))
  cat(sprintf("   fit$loglik=%.4f stop=%s\n", fit$loglik, fit$stop_reason))
  break
}
# Real H20 check: is the final row an E-step AT theta_K (not the lagged one)?
set.seed(5)
d <- drv(brts22, c(3, 2), c(0,0), c(4,4), sample_size = 50L, max_iter = 2L)
m <- d$mcem; K <- d$iterations
cat(sprintf("\nH20 lag check: iter rows fhat = %s ; final row fhat = %.4f at pars (%s)\n", paste(round(m$fhat[1:K],4), collapse=", "), m$fhat[K+1], paste(round(emphasis:::.contract_pars(d$pars, c(0L,0L,0L)),3), collapse=",")))
reps <- replicate(8, { r <- emphasis:::em_cpp(brts = brts22, init_pars = d$pars, sample_size = 50L, maxN = 2000L, max_missing = 1e4, max_lambda = 1e6,
   lower_bound = ex(c(0,0)), upper_bound = ex(c(4,4)), xtol_rel = 1e-3, num_threads = 1L, copy_trees = FALSE, model = c(0L,0L,0L), link = 0L); r$fhat })
cat(sprintf("independent E-steps at returned pars: mean %.4f sd %.4f  -> final row within noise? %s ; lagged row K within noise? %s\n",
    mean(reps), sd(reps), abs(m$fhat[K+1]-mean(reps)) < 3*sd(reps)+0.05, abs(m$fhat[K]-mean(reps)) < 3*sd(reps)+0.05))

## H21 independent replication on a DIFFERENT tree (20 tips, lambda=.6, mu=.3, seed 7).
## 1. BDI fits from a far start with box [0,3] / [0,30] / [0,300]: does "converged" fire during drift?
## 2. Is the small per-iteration step at the wide-box stop point genuine EM drift (same with N=800)
##    or MC noise / a truncated M-step?
## 3. Thinning (dynamic_fresh) sampler, wide box, far start.
## 4. Side claim: default midpoint init (lambda == mu) => e_step_failure with BDI.
## 5. What box does auto_bounds() give for this tree (is the wide box a pipeline default)?
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(TreeSim); library(DDD); library(ape)})
T0 <- proc.time()[3]
set.seed(7)
t <- TreeSim::sim.bd.taxa(20, 1, 0.6, 0.3, complete = FALSE)[[1]]
brts <- sort(as.numeric(ape::branching.times(t)), decreasing = TRUE)
cat(sprintf("tree: n_tips=%d crown=%.3f\n", length(brts) + 1, brts[1]))
ll <- function(p) DDD::bd_loglik(c(p[1], p[2], 0, 0), c(0, 0, 1, 0, 2), brts, 0)
ml <- DDD::bd_ML(brts, initparsopt = c(0.4, 0.1), idparsopt = 1:2, parsfix = c(0, 0), idparsfix = 3:4,
                 cond = 0, btorph = 1, soc = 2, verbose = FALSE)
mle <- c(ml$lambda0, ml$mu0); llmax <- ml$loglik
H <- optimHess(mle, function(p) -ll(p)); se <- sqrt(diag(solve(H)))
cat(sprintf("bd_ML: lambda=%.4f mu=%.4f loglik=%.3f SE=(%.3f,%.3f)\n", mle[1], mle[2], llmax, se[1], se[2]))
tol <- 1e-3
mb <- c(0L,0L,0L)

fit_once <- function(ub, init, samp = "bdi", max_time = 30) {
  t0 <- proc.time()[3]
  f <- estimate_rates(brts, model = "cr", method = "mcem", init_pars = init,
        control = list(lower_bound = c(0,0), upper_bound = c(ub,ub), sampling = samp,
                       num_threads = 1L, max_time = max_time, verbose = FALSE))
  m <- f$details$mcem
  steps <- if (nrow(m) > 1) pmax(abs(diff(m$par1)), abs(diff(m$par5))) else NA
  data.frame(ub = ub, samp = samp, stop = f$details$stop_reason, iters = f$details$iterations,
             lambda = f$pars[1], mu = f$pars[2],
             d_lam_SE = (f$pars[1] - mle[1]) / se[1], d_mu_SE = (f$pars[2] - mle[2]) / se[2],
             ll_deficit = llmax - ll(f$pars), last_delta = tail(m$delta_max, 1),
             last_abs_step = tail(steps, 1), secs = proc.time()[3] - t0, row.names = NULL)
}

## ---- 1. far start, BDI, three boxes ----
init <- c(2.5 * mle[1], 3.5 * mle[2])
cat(sprintf("\n=== 1. BDI far start (%.3f,%.3f), loglik(start)=%.3f (deficit %.2f) ===\n",
            init[1], init[2], ll(init), llmax - ll(init)))
res1 <- rbind(fit_once(3, init, max_time = 25),
              fit_once(30, init), fit_once(30, init),
              fit_once(300, init), fit_once(300, init))
print(res1, digits = 4, row.names = FALSE)
stop_pt <- c(res1$lambda[2], res1$mu[2])

## ---- 2. is the step at the wide-box stop point EM drift or noise? ----
cat(sprintf("\n=== 2. single E+M step from the ub=30 stop point (%.4f,%.4f): N=200 x4 vs N=800 x2 ===\n",
            stop_pt[1], stop_pt[2]))
lb8 <- emphasis:::.expand_pars(c(0,0), mb); ub8 <- emphasis:::.expand_pars(c(30,30), mb)
one_step_bdi <- function(theta, N) {
  th8 <- emphasis:::.expand_pars(theta, mb)
  e <- emphasis:::.augment_tree_bdi(brts, th8, mb, as.integer(N), 1e4, 0L, 1.0)
  lw <- e$weights; w <- exp(lw - max(lw)); w <- w / sum(w) * length(w)
  es <- list(trees = e$trees, weights = w, rejected = 0L, rejected_overruns = 0L,
             rejected_lambda = 0L, rejected_zero_weights = 0L, time = 0, fhat = e$fhat)
  m <- emphasis:::m_cpp(es, th8, "rpd1", lb8, ub8, 1e-3, 1L, mb, 0L, 1.0, NULL)
  as.numeric(m$estimates)[c(1, 5)]
}
s200 <- t(replicate(4, one_step_bdi(stop_pt, 200)))
s800 <- t(replicate(2, one_step_bdi(stop_pt, 800)))
rep_step <- function(s, lab) {
  d <- sweep(s, 2, stop_pt)
  cat(sprintf("%s: mean step=(%+.4f,%+.4f)  sd=(%.4f,%.4f)  |step|/range30=(%.1e,%.1e)  all steps toward MLE? %s\n",
              lab, mean(d[,1]), mean(d[,2]), sd(d[,1]), sd(d[,2]),
              abs(mean(d[,1]))/30, abs(mean(d[,2]))/30,
              all(sign(d[,1]) == sign(mle[1]-stop_pt[1])) && all(sign(d[,2]) == sign(mle[2]-stop_pt[2]))))
}
rep_step(s200, "N=200"); rep_step(s800, "N=800")

## ---- 3. thinning sampler, wide box, far start ----
cat("\n=== 3. thinning (dynamic_fresh), far start, ub=30 vs ub=3 ===\n")
res3 <- rbind(fit_once(30, init, "dynamic_fresh", max_time = 20),
              fit_once(30, init, "dynamic_fresh", max_time = 20),
              fit_once(3, init, "dynamic_fresh", max_time = 20))
print(res3, digits = 4, row.names = FALSE)

## ---- 4. default init (box midpoint, lambda == mu) ----
cat("\n=== 4. default init with symmetric box [0,2]^2 (BDI) ===\n")
f4 <- estimate_rates(brts, model = "cr", method = "mcem",
        control = list(lower_bound = c(0,0), upper_bound = c(2,2), sampling = "bdi",
                       num_threads = 1L, max_time = 30, verbose = FALSE))
cat(sprintf("stop=%s iters=%d pars=(%.4f,%.4f)\n", f4$details$stop_reason, f4$details$iterations, f4$pars[1], f4$pars[2]))
e4 <- tryCatch(emphasis:::.augment_tree_bdi(brts, emphasis:::.expand_pars(c(1,1), mb), mb, 50L, 1e4, 0L, 1.0), error = function(e) e)
cat("direct E-step at (1,1):", if (inherits(e4, "error")) paste("ERROR:", conditionMessage(e4)) else "ok", "\n")

## ---- 5. auto_bounds box for this tree ----
cat("\n=== 5. auto_bounds(cr) box ===\n")
ab <- tryCatch(auto_bounds(t, model = "cr"), error = function(e) e)
if (inherits(ab, "error")) cat("auto_bounds ERROR:", conditionMessage(ab), "\n") else {
  cat("lb:", round(ab$lower_bound, 4), " ub:", round(ab$upper_bound, 4),
      " range:", round(ab$upper_bound - ab$lower_bound, 4), " tol*range:", signif(tol*(ab$upper_bound - ab$lower_bound), 3), "\n")
}
cat(sprintf("\n(total %.0fs)\n", proc.time()[3] - T0))

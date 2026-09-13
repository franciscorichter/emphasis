## H21: constant N=200, convergence = max_j |dtheta_j|/range_j < 1e-3 for 3 consecutive
## iterations; claim: MC noise of one MCEM step exceeds 1e-3*range, so "converged" fires
## on three quiet draws (or never) and `pars` is one draw from the MCEM stationary law.
## Tests on ONE fixed 25-tip CR tree (lambda=.5, mu=.2), rho=1, cond=NULL, num_threads=1:
##  A. single-step MC noise at theta = MLE (DDD::bd_ML): sd(theta_k)/range and
##     P(delta_max < tol) for the BDI (default) and thinning (em_cpp) E/M steps.
##  B. R full BDI fits started AT the MLE: stop_reason, iterations, spread of `pars`,
##     distance to bd_ML, compared with the MLE's own SE (numerical Hessian).
##  C. same as B with ub = 30 (range x10): the rule's meaning depends on the user box.
##  D. a few thinning (dynamic_fresh) fits from the MLE.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(TreeSim); library(DDD); library(ape)})
lam <- 0.5; mu <- 0.2
set.seed(21)
t <- TreeSim::sim.bd.taxa(25, 1, lam, mu, complete = FALSE)[[1]]
brts <- sort(as.numeric(ape::branching.times(t)), decreasing = TRUE)
cat(sprintf("tree: n_tips=%d crown=%.3f\n", length(brts) + 1, brts[1]))

## reference MLE (unconditioned, crown, soc=2) + numerical SE
ll <- function(p) DDD::bd_loglik(c(p[1], p[2], 0, 0), c(0, 0, 1, 0, 2), brts, 0)
ml <- DDD::bd_ML(brts, initparsopt = c(0.4, 0.1), idparsopt = 1:2, parsfix = c(0, 0), idparsfix = 3:4,
                 cond = 0, btorph = 1, soc = 2, verbose = FALSE)
mle <- c(ml$lambda0, ml$mu0)
H <- optimHess(mle, function(p) -ll(p))
se <- sqrt(diag(solve(H)))
cat(sprintf("DDD::bd_ML: lambda=%.4f mu=%.4f loglik=%.3f ; SE(lambda)=%.4f SE(mu)=%.4f\n",
            mle[1], mle[2], ml$loglik, se[1], se[2]))

lb <- c(0, 0); ub <- c(3, 3); range3 <- ub - lb; tol <- 1e-3
lb8 <- emphasis:::.expand_pars(lb, c(0L,0L,0L)); ub8 <- emphasis:::.expand_pars(ub, c(0L,0L,0L))
mle8 <- emphasis:::.expand_pars(mle, c(0L,0L,0L))

## ---------- A. one-step MC noise at theta = MLE ----------
one_step_bdi <- function(theta8, N = 200L) {
  e <- emphasis:::.augment_tree_bdi(brts, theta8, c(0L,0L,0L), N, 1e4, 0L, 1.0)
  lw <- e$weights; w <- exp(lw - max(lw)); w <- w / sum(w) * length(w)
  es <- list(trees = e$trees, weights = w, rejected = 0L, rejected_overruns = 0L,
             rejected_lambda = 0L, rejected_zero_weights = 0L, time = 0, fhat = e$fhat)
  m <- emphasis:::m_cpp(es, theta8, "rpd1", lb8, ub8, 1e-3, 1L, c(0L,0L,0L), 0L, 1.0, NULL)
  as.numeric(m$estimates)[c(1, 5)]
}
one_step_thin <- function(theta8, N = 200L) {
  r <- emphasis:::em_cpp(brts, theta8, N, 2000L, 1e4, 1e6, lb8, ub8, 1e-3, 1L, FALSE,
                         c(0L,0L,0L), 0L, 1.0, NULL)
  as.numeric(r$estimates)[c(1, 5)]
}
report_noise <- function(th, label) {
  sdv <- apply(th, 2, sd); mn <- colMeans(th)
  d <- abs(diff(th)) / matrix(range3, nrow(th) - 1, 2, byrow = TRUE)   # successive-draw delta
  dmax <- apply(d, 1, max)
  p_quiet <- mean(dmax < tol)
  cat(sprintf("%-10s R=%d  mean=(%.4f,%.4f)  sd=(%.4f,%.4f)  sd/range=(%.1e,%.1e)  sd/SE=(%.2f,%.2f)\n",
              label, nrow(th), mn[1], mn[2], sdv[1], sdv[2], sdv[1]/range3[1], sdv[2]/range3[2], sdv[1]/se[1], sdv[2]/se[2]))
  cat(sprintf("           successive delta_max: median=%.1e  min=%.1e  P(delta_max<tol)=%.2f  -> P(3 in a row)~%.3f ; tol=%.0e\n",
              median(dmax), min(dmax), p_quiet, p_quiet^3, tol))
  invisible(list(sd = sdv, p = p_quiet))
}
cat("\n=== A. single E+M step from theta=MLE, N=200 ===\n")
tA <- proc.time()[3]
thB <- t(replicate(20, one_step_bdi(mle8)))
nB <- report_noise(thB, "BDI")
thT <- t(replicate(12, one_step_thin(mle8)))
nT <- report_noise(thT, "thinning")
cat(sprintf("(A took %.0fs)\n", proc.time()[3] - tA))

## ---------- B. full BDI fits started at the MLE, box [0,3]^2 ----------
fit_once <- function(samp, ub, init, max_time) {
  t0 <- proc.time()[3]
  f <- estimate_rates(brts, model = "cr", method = "mcem", init_pars = init,
        control = list(lower_bound = c(0,0), upper_bound = ub, sampling = samp,
                       num_threads = 1L, max_time = max_time, verbose = FALSE))
  m <- f$details$mcem
  data.frame(stop = f$details$stop_reason, iters = f$details$iterations,
             lambda = f$pars[1], mu = f$pars[2], loglik = f$loglik,
             frac_quiet = mean(m$delta_max < tol), med_delta = median(m$delta_max),
             secs = proc.time()[3] - t0, row.names = NULL)
}
summarise_fits <- function(res, label) {
  cat(sprintf("\n--- %s : %d fits ---\n", label, nrow(res)))
  print(res, digits = 4, row.names = FALSE)
  cat("stop_reason:", paste(names(table(res$stop)), table(res$stop), sep = "=", collapse = " "), "\n")
  cat(sprintf("iterations: median=%d range=[%d,%d]\n", as.integer(median(res$iters)), min(res$iters), max(res$iters)))
  cat(sprintf("pars: mean=(%.4f,%.4f)  sd=(%.4f,%.4f)  sd/SE_MLE=(%.2f,%.2f)\n",
              mean(res$lambda), mean(res$mu), sd(res$lambda), sd(res$mu), sd(res$lambda)/se[1], sd(res$mu)/se[2]))
  cat(sprintf("|mean - bd_ML| = (%.4f,%.4f) ; max |pars - bd_ML| over fits = (%.4f,%.4f)\n",
              abs(mean(res$lambda) - mle[1]), abs(mean(res$mu) - mle[2]),
              max(abs(res$lambda - mle[1])), max(abs(res$mu - mle[2]))))
  cat(sprintf("loglik at returned pars vs bd_ML max: mean deficit = %.4f (max %.4f)\n",
              mean(ml$loglik - sapply(seq_len(nrow(res)), function(i) ll(c(res$lambda[i], res$mu[i])))),
              max(ml$loglik - sapply(seq_len(nrow(res)), function(i) ll(c(res$lambda[i], res$mu[i]))))))
}
cat("\n=== B. full BDI fits, init = MLE, box [0,3]^2, default control ===\n")
tB <- proc.time()[3]
resB <- do.call(rbind, lapply(1:12, function(i) fit_once("bdi", c(3,3), mle, 40)))
summarise_fits(resB, "BDI, ub=3")
cat(sprintf("(B took %.0fs)\n", proc.time()[3] - tB))

cat("\n=== C. same, box [0,30]^2 (range x10 => tol means |dtheta| < 0.03) ===\n")
tC <- proc.time()[3]
resC <- do.call(rbind, lapply(1:8, function(i) fit_once("bdi", c(30,30), mle, 40)))
summarise_fits(resC, "BDI, ub=30")
cat(sprintf("(C took %.0fs)\n", proc.time()[3] - tC))

cat("\n=== D. thinning (dynamic_fresh) fits, init = MLE, box [0,3]^2 ===\n")
tD <- proc.time()[3]
resD <- do.call(rbind, lapply(1:3, function(i) fit_once("dynamic_fresh", c(3,3), mle, 45)))
summarise_fits(resD, "thinning, ub=3")
cat(sprintf("(D took %.0fs)\n", proc.time()[3] - tD))

## H4 replication (independent): does the M-step conditional penalty get diluted by Σw?
## Varied vs the verifier: different tree (seed 11, 30 tips, lambda=0.8, mu=0.3),
## two sample sizes (N = 50, 400) to test that dilution scales with N,
## thinning-path weights exp(lw - max) as E_step.cpp:133-139 hands them,
## and a short full-fit comparison (BDI, 40 iterations, 3 replicates).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(DDD); library(TreeSim)})
options(digits = 5)

set.seed(11)
lam_true <- 0.8; mu_true <- 0.3; ntip <- 30
tr   <- TreeSim::sim.bd.taxa(n = ntip, numbsim = 1, lambda = lam_true, mu = mu_true,
                             complete = FALSE)[[1]]
brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
Tc   <- brts[1]
cat("tree: ", ntip, "tips, crown age", round(Tc, 3), "\n")

logPsurv <- function(lam, mu, T = Tc) {
  if (abs(lam - mu) < 1e-8) mu <- lam - 1e-8
  r  <- lam - mu
  p0 <- mu * (1 - exp(-r * T)) / (lam - mu * exp(-r * T))
  2 * log(max(1 - p0, 1e-300))
}
cond8 <- function(pars8) logPsurv(pars8[1], pars8[5])
model_bin <- c(0L, 0L, 0L)
lb8 <- c(1e-3, 0, 0, 0, 1e-3, 0, 0, 0)
ub8 <- c(5,    0, 0, 0, 5,    0, 0, 0)

dd_ll <- function(p, cond) DDD::bd_loglik(pars1 = c(p[1], p[2], 0, 0),
                                          pars2 = c(0, cond, 0, 0, 2),
                                          brts = brts, missnumspec = 0)
mle_ddd <- function(cond) optim(c(0.5, 0.2), function(p) if (any(p <= 0)) 1e10 else -dd_ll(p, cond),
                                control = list(reltol = 1e-10, maxit = 5000))$par
mle0 <- mle_ddd(0); mle1 <- mle_ddd(1)
cat(sprintf("DDD cond=0 MLE: %.4f %.4f | cond=1 MLE: %.4f %.4f\n", mle0[1], mle0[2], mle1[1], mle1[2]))
p <- c(0.6, 0.2)
cat(sprintf("DDD (cond1-cond0) at (0.6,0.2) = %.5f vs -logPsurv = %.5f\n",
            dd_ll(p, 1) - dd_ll(p, 0), -logPsurv(p[1], p[2])))

theta0 <- c(0.8, 0, 0, 0, 0.3, 0, 0, 0)

mstep <- function(trees, w, cond = NULL) {
  es <- list(trees = trees, weights = w, rejected = 0L, rejected_overruns = 0L,
             rejected_lambda = 0L, rejected_zero_weights = 0L, time = 0, fhat = 0)
  r <- emphasis:::m_cpp(e_step = es, init_pars = theta0, plugin = "rpd1",
                        lower_bound = lb8, upper_bound = ub8, xtol_rel = 1e-6,
                        num_threads = 1L, model = model_bin, link = 0L, rho = 1,
                        rconditional = cond)
  r$estimates[c(1, 5)]
}
ref_argmax <- function(trees, w, pen) {
  Qn <- function(pp) {
    lf <- emphasis:::eval_logf(c(pp[1], 0, 0, 0, pp[2], 0, 0, 0), trees,
                               model = model_bin, link = 0L, rho = 1)$logf
    sum(w * lf) / sum(w)
  }
  optim(c(0.6, 0.2), function(pp) {
    if (any(pp < 1e-3) || any(pp > 5)) return(1e10)
    -(Qn(pp) - pen * logPsurv(pp[1], pp[2])) }, control = list(reltol = 1e-12, maxit = 5000))$par
}

## ------------------------------------------------ BDI sampler, N = 50 and N = 400
for (N in c(50L, 400L)) {
  e <- emphasis:::.augment_tree_bdi(brts, theta0, model_bin = model_bin,
                                    sample_size = N, max_missing = 1e4, link = 0L, rho = 1)
  lw <- e$weights
  w <- exp(lw - max(lw)); w <- w / sum(w) * length(w)
  cat(sprintf("\n== BDI N=%d, Σw=%.1f ==\n", N, sum(w)))
  a  <- mstep(e$trees, w)
  b  <- mstep(e$trees, w, cond8)          # package behaviour
  c_ <- mstep(e$trees, w / N, cond8)      # self-normalised
  d  <- mstep(e$trees, w * 5, cond8)
  r0 <- ref_argmax(e$trees, w, 0); r1 <- ref_argmax(e$trees, w, 1); rN <- ref_argmax(e$trees, w, 1 / N)
  cat(sprintf("no cond           : %.4f %.4f\n", a[1], a[2]))
  cat(sprintf("cond, w (pkg)     : %.4f %.4f\n", b[1], b[2]))
  cat(sprintf("cond, w/N         : %.4f %.4f\n", c_[1], c_[2]))
  cat(sprintf("cond, 5w          : %.4f %.4f\n", d[1], d[2]))
  cat(sprintf("ref argmax Qn     : %.4f %.4f\n", r0[1], r0[2]))
  cat(sprintf("ref argmax Qn-logP: %.4f %.4f\n", r1[1], r1[2]))
  cat(sprintf("ref argmax Qn-logP/N: %.4f %.4f\n", rN[1], rN[2]))
}

## ------------------------------------------------ thinning sampler weights, N = 200
N <- 200L
a <- emphasis:::augment_trees(brts, theta0, sample_size = N, maxN = 5000L, max_missing = 1e4,
                              max_lambda = 1e6, num_threads = 1L, model = model_bin, link = 0L, rho = 1)
lf <- emphasis:::eval_logf(theta0, a$trees, model = model_bin, link = 0L, rho = 1)$logf
lwt <- lf - a$logg
wt  <- exp(lwt - max(lwt))
cat(sprintf("\n== thinning N=%d, Σw=%.2f, ESS=%.1f ==\n", N, sum(wt), sum(wt)^2 / sum(wt^2)))
a0 <- mstep(a$trees, wt); b0 <- mstep(a$trees, wt, cond8); c0 <- mstep(a$trees, wt / sum(wt), cond8)
r0 <- ref_argmax(a$trees, wt, 0); r1 <- ref_argmax(a$trees, wt, 1); rS <- ref_argmax(a$trees, wt, 1 / sum(wt))
cat(sprintf("no cond           : %.4f %.4f\n", a0[1], a0[2]))
cat(sprintf("cond, w (pkg)     : %.4f %.4f\n", b0[1], b0[2]))
cat(sprintf("cond, w/Σw        : %.4f %.4f\n", c0[1], c0[2]))
cat(sprintf("ref argmax Qn     : %.4f %.4f\n", r0[1], r0[2]))
cat(sprintf("ref argmax Qn-logP: %.4f %.4f\n", r1[1], r1[2]))
cat(sprintf("ref argmax Qn-logP/Σw: %.4f %.4f\n", rS[1], rS[2]))

## ------------------------------------------------ short full BDI fits, 3 replicates each
run_fit <- function(cond, iters = 40, N = 200L) {
  r <- emphasis:::.mcem_bdi(brts = brts, pars = theta0, sample_size = N, max_missing = 1e4,
                            lower_bound = lb8, upper_bound = ub8, max_iter = iters,
                            xtol = 1e-4, tol = 1e-4, patience = 3, num_threads = 1L,
                            conditional = cond, model = model_bin, link = 0L,
                            max_time = 120, rho = 1)
  m <- r$mcem
  cols <- grep("^(par|theta|lambda|beta|gamma|mu)", names(m), value = TRUE)
  colMeans(utils::tail(m[, cols, drop = FALSE], 10))[c(1, 5)]
}
cat("\n== full BDI MCEM fits, mean of last 10 iterates ==\n")
for (cc in c("none", "exact")) for (k in 1:3) {
  f <- run_fit(if (cc == "none") NULL else cond8)
  cat(sprintf("cond=%-5s rep %d: %.4f %.4f\n", cc, k, f[1], f[2]))
}
cat(sprintf("targets: DDD cond=0 (%.4f, %.4f)   DDD cond=1 (%.4f, %.4f)\n",
            mle0[1], mle0[2], mle1[1], mle1[2]))

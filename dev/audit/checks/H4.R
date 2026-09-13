## H4 — M-step conditional penalty is not scaled by the weight sum.
## Objective in src/M_step.cpp:61 is  -Σ_i w_i log f_i(θ) + log P_tree(θ)
## with w unnormalised (thinning: Σw ∈ [1,N]; BDI: Σw = N).  The
## self-normalised conditioned EM objective needs  -Σ_i (w_i/Σw) log f_i + log P_tree,
## i.e. the penalty is effectively divided by Σw.
##
## Tests:
##  A. Direct m_cpp calls on one fixed E-step sample: rescale w by 1/N, 1, 10 with
##     the same conditional -> argmin must be invariant if the code were correct.
##  B. R reference argmax of the normalised objective via eval_logf + optim.
##  C. Full CR MCEM fits (BDI and thinning) with an EXACT crown-survival
##     conditional vs the DDD::bd_loglik cond=0 / cond=1 maximisers.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(DDD); library(TreeSim)})
options(digits = 5)

set.seed(5)
lam_true <- 1.0; mu_true <- 0.5; ntip <- 25
tr   <- TreeSim::sim.bd.taxa(n = ntip, numbsim = 1, lambda = lam_true, mu = mu_true,
                             complete = FALSE)[[1]]
brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
Tc   <- brts[1]
cat("tree: ", ntip, "tips, crown age", round(Tc, 3), "\n")

## exact log P(both crown lineages survive to present) for CR
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

## ---------------------------------------------------------------- DDD reference
dd_ll <- function(p, cond) DDD::bd_loglik(pars1 = c(p[1], p[2], 0, 0),
                                          pars2 = c(0, cond, 0, 0, 2),
                                          brts = brts, missnumspec = 0)
mle_ddd <- function(cond) {
  o <- optim(c(0.5, 0.2), function(p) if (any(p <= 0)) 1e10 else -dd_ll(p, cond),
             method = "Nelder-Mead", control = list(reltol = 1e-10, maxit = 5000))
  o$par
}
mle0 <- mle_ddd(0); mle1 <- mle_ddd(1)
cat(sprintf("DDD cond=0 MLE: lambda=%.4f mu=%.4f\n", mle0[1], mle0[2]))
cat(sprintf("DDD cond=1 MLE: lambda=%.4f mu=%.4f\n", mle1[1], mle1[2]))
## sanity: DDD cond=1 - cond=0 loglik difference equals -logPsurv at a point
p <- c(0.7, 0.4)
cat(sprintf("check DDD cond1-cond0 at (0.7,0.4): %.5f  vs -logPsurv: %.5f\n",
            dd_ll(p, 1) - dd_ll(p, 0), -logPsurv(p[1], p[2])))

## ---------------------------------------------------------------- A. m_cpp scale test
theta0 <- c(1.0, 0, 0, 0, 0.5, 0, 0, 0)
N <- 200L
e <- emphasis:::.augment_tree_bdi(brts, theta0, model_bin = model_bin,
                                  sample_size = N, max_missing = 1e4, link = 0L, rho = 1)
lw <- e$weights
w_norm <- exp(lw - max(lw)); w_norm <- w_norm / sum(w_norm) * length(w_norm)  # as bdi.R:796-799
cat(sprintf("BDI weights: sum = %.3f (N = %d), sd(lw) = %.2e\n", sum(w_norm), N, sd(lw)))

mstep <- function(w, cond = NULL) {
  es <- list(trees = e$trees, weights = w, rejected = 0L, rejected_overruns = 0L,
             rejected_lambda = 0L, rejected_zero_weights = 0L, time = 0, fhat = e$fhat)
  r <- emphasis:::m_cpp(e_step = es, init_pars = theta0, plugin = "rpd1",
                        lower_bound = lb8, upper_bound = ub8, xtol_rel = 1e-6,
                        num_threads = 1L, model = model_bin, link = 0L, rho = 1,
                        rconditional = cond)
  r$estimates[c(1, 5)]
}
th_a  <- mstep(w_norm)               # no conditional
th_a2 <- mstep(w_norm / N)           # no conditional, rescaled  (must equal th_a)
th_b  <- mstep(w_norm, cond8)        # conditional, Σw = N   (what the package does)
th_c  <- mstep(w_norm / N, cond8)    # conditional, Σw = 1   (self-normalised form)
th_d  <- mstep(w_norm * 10, cond8)   # conditional, Σw = 10N

## ---------------------------------------------------------------- B. R reference
Qn <- function(p) {
  lf <- emphasis:::eval_logf(c(p[1], 0, 0, 0, p[2], 0, 0, 0), e$trees,
                             model = model_bin, link = 0L, rho = 1)$logf
  sum(w_norm * lf) / sum(w_norm)
}
ref <- function(pen) optim(c(0.7, 0.4), function(p) {
  if (any(p < 1e-3) || any(p > 5)) return(1e10)
  -(Qn(p) - pen * logPsurv(p[1], p[2])) }, control = list(reltol = 1e-12, maxit = 5000))$par
ref_uncond <- ref(0)
ref_cond   <- ref(1)
ref_diluted <- ref(1 / N)

cat("\n== A/B: M-step argmin on one fixed E-step sample (lambda, mu) ==\n")
cat(sprintf("no cond,     w (Σw=N)      : %.4f %.4f\n", th_a[1], th_a[2]))
cat(sprintf("no cond,     w/N (Σw=1)    : %.4f %.4f   [invariance without cond]\n", th_a2[1], th_a2[2]))
cat(sprintf("cond,        w (Σw=N)      : %.4f %.4f   <- package behaviour\n", th_b[1], th_b[2]))
cat(sprintf("cond,        w/N (Σw=1)    : %.4f %.4f   <- self-normalised objective\n", th_c[1], th_c[2]))
cat(sprintf("cond,        10w (Σw=10N)  : %.4f %.4f\n", th_d[1], th_d[2]))
cat(sprintf("R ref argmax Q_norm            : %.4f %.4f\n", ref_uncond[1], ref_uncond[2]))
cat(sprintf("R ref argmax Q_norm - logP     : %.4f %.4f\n", ref_cond[1], ref_cond[2]))
cat(sprintf("R ref argmax Q_norm - logP/N   : %.4f %.4f\n", ref_diluted[1], ref_diluted[2]))
cat(sprintf("Q_norm(th_b) - logP(th_b) = %.5f ; Q_norm(th_c) - logP(th_c) = %.5f\n",
            Qn(th_b) - logPsurv(th_b[1], th_b[2]), Qn(th_c) - logPsurv(th_c[1], th_c[2])))

## ---------------------------------------------------------------- C. full MCEM fits
run_fit <- function(sampling, cond, iters = 60) {
  if (sampling == "bdi") {
    r <- emphasis:::.mcem_bdi(brts = brts, pars = theta0, sample_size = N, max_missing = 1e4,
                              lower_bound = lb8, upper_bound = ub8, max_iter = iters,
                              xtol = 1e-4, tol = 1e-4, patience = 3, num_threads = 1L,
                              conditional = cond, model = model_bin, link = 0L,
                              max_time = 200, rho = 1)
  } else {
    r <- emphasis:::.mcem_dynamic_fresh(brts = brts, pars = theta0, sample_size = N,
                              maxN = 5000L, max_missing = 1e4,
                              lower_bound = lb8, upper_bound = ub8, max_iter = iters,
                              xtol = 1e-4, tol = 1e-4, patience = 3, num_threads = 1L,
                              conditional = cond, model = model_bin, link = 0L,
                              max_time = 200, rho = 1)
  }
  ## average the last 10 iterates to tame MC noise
  m <- r$mcem
  cols <- grep("^(par|theta|lambda|beta|gamma|mu)", names(m), value = TRUE)
  list(last = as.numeric(r$pars)[c(1, 5)],
       avg = if (length(cols) >= 5) colMeans(utils::tail(m[, cols, drop = FALSE], 10))[c(1, 5)] else NA,
       n = nrow(m), stop = r$stop_reason)
}
cat("\n== C: full MCEM fits (lambda, mu), last iterate / mean of last 10 ==\n")
for (s in c("bdi", "dynamic_fresh")) for (cc in c("none", "exact")) {
  f <- run_fit(s, if (cc == "none") NULL else cond8)
  cat(sprintf("%-14s cond=%-5s: last %.4f %.4f | avg10 %.4f %.4f | iters %d (%s)\n",
              s, cc, f$last[1], f$last[2], f$avg[1], f$avg[2], f$n, f$stop))
}
cat(sprintf("targets: DDD cond=0 (%.4f, %.4f)   DDD cond=1 (%.4f, %.4f)\n",
            mle0[1], mle0[2], mle1[1], mle1[2]))

## thinning-path Σw at theta0 (weights as E_step.cpp:133-139 hands them to M-step)
a <- emphasis:::augment_trees(brts, theta0, sample_size = N, maxN = 5000L, max_missing = 1e4,
                              max_lambda = 1e6, num_threads = 1L, model = model_bin, link = 0L, rho = 1)
lf <- emphasis:::eval_logf(theta0, a$trees, model = model_bin, link = 0L, rho = 1)$logf
lwt <- lf - a$logg
wt  <- exp(lwt - max(lwt))
cat(sprintf("thinning Σw at theta0 = %.2f of N = %d (ESS %.1f)\n", sum(wt), N, sum(wt)^2 / sum(wt^2)))

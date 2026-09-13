## H2: BDI ignores rho.  Self-contained check.
## Run: Rscript dev/audit/checks/H2.R
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape) })
options(width = 120)

## ---- tree: fixed 16-tip CR tree (R RNG is seedable; the C++ RNG is not) ----
set.seed(11)
phy  <- ape::rphylo(16, birth = 0.5, death = 0.3)
brts <- sort(ape::branching.times(phy), decreasing = TRUE)
n    <- length(brts) + 1L
cat(sprintf("Tree: n_tips = %d, crown age = %.4f\n", n, brts[1]))

## ---- closed form: Nee et al. crown-conditioned density with Bernoulli sampling rho ----
## p1(t): density that a lineage alive t before the present leaves exactly one
## sampled tip; p0(t): no sampled tip.  Labelled-history convention up to a
## theta-independent constant: log L = 2 log p1(tc) + sum_{i>=2} [log lam + log p1(x_i)].
nee_rho <- function(lam, mu, rho, brts) {
  r <- lam - mu
  E <- exp(-r * brts)
  den <- rho * lam + (lam * (1 - rho) - mu) * E
  logp1 <- log(rho) + 2 * log(r) - r * brts - 2 * log(den)
  2 * logp1[1] + sum(log(lam) + logp1[-1])
}

## ---- part A: BDI's own draws scored at rho = 1 and rho = 0.5 ----
lam <- 0.5; mu <- 0.3; rho <- 0.5
bdi1 <- emphasis:::.augment_tree_bdi(brts, pars = c(lam, mu), model_bin = c(0L,0L,0L),
                                     sample_size = 200L, link = 0L, rho = 1.0)
n_unsamp <- sapply(bdi1$trees, function(tr) sum(tr$t_ext == 5e10))
cat(sprintf("A. BDI draws: unsampled-extant nodes per tree: max = %d (must be 0 if BDI never proposes them)\n",
            max(n_unsamp)))
ev05 <- emphasis:::eval_logf(emphasis:::.expand_pars(c(lam, mu), c(0L,0L,0L)), bdi1$trees,
                  model = c(0L,0L,0L), link = 0L, rho = rho)
ev1  <- emphasis:::eval_logf(emphasis:::.expand_pars(c(lam, mu), c(0L,0L,0L)), bdi1$trees,
                  model = c(0L,0L,0L), link = 0L, rho = 1.0)
d <- ev05$logf - ev1$logf
cat(sprintf("A. logf(rho=0.5) - logf(rho=1) on the same BDI trees: range [%.10f, %.10f]; n*log(rho) = %.10f\n",
            min(d), max(d), n * log(rho)))
cat(sprintf("A. sd of BDI log-weights at rho=1: %.3e (zero-variance IS expected under CR)\n",
            sd(bdi1$weights)))
bdi05 <- emphasis:::.augment_tree_bdi(brts, pars = c(lam, mu), model_bin = c(0L,0L,0L),
                                      sample_size = 200L, link = 0L, rho = rho)
cat(sprintf("A. BDI fhat(rho=0.5) = %.6f ; BDI fhat(rho=1) = %.6f ; diff = %.6f ; n*log(rho) = %.6f\n",
            bdi05$fhat, bdi1$fhat, bdi05$fhat - bdi1$fhat, n * log(rho)))
cat(sprintf("A. .bdi_supported(cr, linear) = %s (takes no rho argument: formals = %s)\n",
            emphasis:::.bdi_supported(c(0L,0L,0L), 0L),
            paste(names(formals(emphasis:::.bdi_supported)), collapse = ",")))

## ---- part B: the theta-independent offset between emphasis f and the Nee form ----
grid <- rbind(c(0.5, 0.3), c(0.4, 0.1), c(0.7, 0.5), c(0.25, 0.05), c(0.6, 0.2))
cat("\nB. rho = 1: BDI fhat - nee (must be a constant across theta)\n")
C1 <- apply(grid, 1, function(p) {
  b <- emphasis:::.augment_tree_bdi(brts, pars = p, model_bin = c(0L,0L,0L),
                                    sample_size = 50L, link = 0L, rho = 1.0)
  b$fhat - nee_rho(p[1], p[2], 1, brts)
})
print(cbind(lam = grid[,1], mu = grid[,2], offset = C1))
C <- mean(C1)
cat(sprintf("   spread of offset = %.2e  ->  C = %.6f\n", diff(range(C1)), C))

## thinning fhat helper (C++ proposal, rho-aware), replicated
thin_fhat <- function(p, rho, N = 1000L, reps = 5L) {
  sapply(seq_len(reps), function(i) {
    a <- emphasis:::.augment_tree_internal(brts, pars = p, model_bin = c(0L,0L,0L),
                                           sample_size = N, maxN = 50L * N,
                                           link = 0L, rho = rho, num_threads = 1L)
    lw <- a$logf - a$logg
    m <- max(lw); log(mean(exp(lw - m))) + m
  })
}

cat("\nB. rho = 0.5: fhat - nee(rho=0.5) - C  (0 = targets the rho-model)\n")
res <- t(apply(grid, 1, function(p) {
  target <- nee_rho(p[1], p[2], rho, brts) + C
  b <- emphasis:::.augment_tree_bdi(brts, pars = p, model_bin = c(0L,0L,0L),
                                    sample_size = 50L, link = 0L, rho = rho)
  th <- thin_fhat(p, rho)
  c(lam = p[1], mu = p[2], target = target,
    bdi_gap = b$fhat - target,
    thin_gap_mean = mean(th) - target, thin_gap_sd = sd(th),
    n_log_rho_minus_true_shift = n * log(rho) - (nee_rho(p[1], p[2], rho, brts) - nee_rho(p[1], p[2], 1, brts)))
}))
print(round(res, 4))

## ---- part B': independent derivation check of the Nee-rho form ----
## Stadler (2009): tree density under (lam, mu, rho) = rho^2 * density under (rho*lam, mu - lam(1-rho), 1)
chk <- nee_rho(lam, mu, rho, brts) - nee_rho(rho * lam, mu - lam * (1 - rho), 1, brts) - 2 * log(rho)
cat(sprintf("\nB'. Stadler transform identity residual: %.3e\n", chk))
## and the thinning estimator at (lam,mu,rho) vs BDI (exact) at the transformed pars
b_tr <- emphasis:::.augment_tree_bdi(brts, pars = c(rho * lam, mu - lam * (1 - rho)),
                                     model_bin = c(0L,0L,0L), sample_size = 50L, link = 0L, rho = 1.0)
th <- thin_fhat(c(lam, mu), rho, N = 2000L, reps = 5L)
cat(sprintf("B'. thinning fhat(0.5,0.3,rho=0.5): mean %.4f (sd %.4f); BDI fhat(0.25,0.05,rho=1) + 2log(0.5) = %.4f; BDI fhat(0.5,0.3,rho=0.5) = %.4f\n",
            mean(th), sd(th), b_tr$fhat + 2 * log(rho), bdi05$fhat))

## ---- part C: effect on estimates ----
mle <- function(rho) {
  o <- optim(c(0.4, 0.2), function(p) if (p[1] <= p[2] || p[2] < 0) 1e10 else -nee_rho(p[1], p[2], rho, brts),
             method = "L-BFGS-B", lower = c(1e-3, 0), upper = c(5, 5))
  o$par
}
m1 <- mle(1); m05 <- mle(rho)
cat(sprintf("\nC. closed-form MLE  rho=1: lam=%.4f mu=%.4f | rho=0.5: lam=%.4f mu=%.4f\n", m1[1], m1[2], m05[1], m05[2]))

fit <- function(sampling, reps = 3L) {
  t(sapply(seq_len(reps), function(i) {
    f <- estimate_rates(brts, method = "mcem", model = "cr",
                        init_pars = c(0.4, 0.2),
                        control = list(rho = rho, sampling = sampling, sample_size = 200L,
                                       max_iter = 40L, max_time = 90, num_threads = 1L,
                                       lower_bound = c(1e-3, 0), upper_bound = c(3, 3)))
    c(lam = f$pars[1], mu = f$pars[2], loglik = f$loglik)
  }))
}
cat("C. estimate_rates(rho = 0.5, sampling = 'bdi'):\n");           print(round(fit("bdi"), 4))
cat("C. estimate_rates(rho = 0.5, sampling = 'dynamic_fresh'):\n"); print(round(fit("dynamic_fresh"), 4))

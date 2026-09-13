## H5 (separation from H7) — score augmentations drawn by an EXACT R implementation of the
## intended thinning process (intensity n(t)·lambda·(1-e^{-mu(T-t)}), TruncExp(mu) lifetimes)
## with the package's own f and q (eval_logf -> Model::loglik / Model::sampling_prob).
## If fhat then matches the closed forms across theta, the -log(2·tips+Ne) term in q is the
## correct parent-marginalisation factor and the bias seen in H5.R Part B/C is the sampler's
## envelope defect (H7), not the weight.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({ library(emphasis); library(DDD) })
set.seed(55)

p0 <- function(s, la, mu) { r <- la - mu
  if (abs(r) < 1e-12) la * s / (1 + la * s) else mu * (exp(r * s) - 1) / (la * exp(r * s) - mu) }
log_p1 <- function(Tt, la, mu) { r <- la - mu
  if (abs(r) < 1e-12) -2 * log(1 + la * Tt) else 2 * log(r) + r * Tt - 2 * log(la * exp(r * Tt) - mu) }
nee_lab <- function(brts, la, mu) {
  Tt <- brts[1]; s <- sort(Tt - brts[-1]); bounds <- c(0, s, Tt); nseg <- length(bounds) - 1
  acc <- (length(brts) - 1) * log(la)
  for (k in seq_len(nseg)) acc <- acc - (1 + k) * integrate(function(t) (la + mu) - 2 * la * sapply(t, p0, la = la, mu = mu),
                                                          bounds[k], bounds[k + 1], rel.tol = 1e-10)$value
  acc
}

## exact Ogata thinning of the intended proposal process (CR, rho = 1)
sample_aug <- function(brts, la, mu) {
  Tt <- brts[1]; s_obs <- sort(Tt - brts[-1])
  ts_v <- numeric(0); te_v <- numeric(0)
  n_at <- function(t) 2 + sum(s_obs <= t) + sum(ts_v <= t & te_v > t)
  cbt <- 0
  while (cbt < Tt) {
    nxt <- min(c(s_obs[s_obs > cbt], Tt))
    L <- n_at(cbt) * la * (1 - exp(-mu * (Tt - cbt)))     # dominates on [cbt, nxt): n can only drop, survival factor decreases
    if (L <= 0) { cbt <- nxt; next }
    cand <- cbt - log(runif(1)) / L
    if (cand >= nxt) { cbt <- nxt; next }
    rate <- n_at(cand) * la * (1 - exp(-mu * (Tt - cand)))
    if (runif(1) < rate / L) {
      R <- Tt - cand
      te <- cand - log(1 - runif(1) * (1 - exp(-mu * R))) / mu   # TruncExp(mu) on (0, R)
      ts_v <- c(ts_v, cand); te_v <- c(te_v, te)
    }
    cbt <- cand
  }
  b <- c(s_obs, Tt, ts_v, te_v)
  tx <- c(rep(1e11, length(s_obs) + 1), te_v, rep(0, length(te_v)))
  o <- order(b); b <- b[o]; tx <- tx[o]
  step <- ifelse(tx == 0, -1, 1)
  n <- 2 + c(0, cumsum(step)[-length(step)])
  data.frame(brts = b, n = n, t_ext = tx, pd = 0, tip_start = 0, focal_tip_start = 0, id = -1L, parent_id = -1L)
}
fhat_exact <- function(brts, la, mu, Ns) {
  trees <- replicate(Ns, sample_aug(brts, la, mu), simplify = FALSE)
  ev <- emphasis:::eval_logf(c(la, 0, 0, 0, mu, 0, 0, 0), trees)
  m <- sapply(trees, function(df) sum(df$t_ext == 0))
  c(fhat = emphasis:::.is_fhat(ev$logf, ev$logg), ess = emphasis:::.ess_from_lw(ev$logf - ev$logg), mean_m = mean(m), P0 = mean(m == 0))
}

cat("==== 2-tip tree: exact-sampler fhat vs 2 log p1(T) ====\n")
grid <- rbind(c(1, 1, 0.1), c(1, 1, 0.6), c(1, 1, 1.0), c(1, 2, 1.8), c(3, 1, 0.8), c(3, 0.6, 0.5))
resB <- t(apply(grid, 1, function(g) {
  Tt <- g[1]; la <- g[2]; mu <- g[3]
  reps <- replicate(3, fhat_exact(Tt, la, mu, 20000))
  ex <- 2 * log_p1(Tt, la, mu)
  Iseg <- function(a, b) (b - a) - (exp(-mu * (Tt - b)) - exp(-mu * (Tt - a))) / mu
  c(T = Tt, lambda = la, mu = mu, exact = ex, fhat = mean(reps["fhat", ]), sd_reps = sd(reps["fhat", ]), gap = mean(reps["fhat", ]) - ex,
    mean_m = mean(reps["mean_m", ]), P0_mc = mean(reps["P0", ]), P0_an = exp(-2 * la * Iseg(0, Tt)))
}))
print(round(resB, 4))

cat("\n==== 10-tip CR tree: exact-sampler fhat - Nee labelled / - DDD::bd_loglik ====\n")
brts10 <- c(10, 7.3, 6.1, 4.8, 3.2, 2.5, 1.7, 0.9, 0.4)
gridC <- rbind(c(0.2, 0.05), c(0.2, 0.15), c(0.3, 0.25), c(0.4, 0.35), c(0.25, 0.28))
resC <- t(apply(gridC, 1, function(g) {
  la <- g[1]; mu <- g[2]
  reps <- replicate(2, fhat_exact(brts10, la, mu, 10000))
  nee <- nee_lab(brts10, la, mu)
  ddd <- DDD::bd_loglik(pars1 = c(la, mu, 0, 0), pars2 = c(0, 0, 0, 0, 2), brts = brts10, missnumspec = 0)
  c(lambda = la, mu = mu, fhat = mean(reps["fhat", ]), sd_reps = sd(reps["fhat", ]), ess = mean(reps["ess", ]), mean_m = mean(reps["mean_m", ]),
    nee_lab = nee, gap_nee = mean(reps["fhat", ]) - nee, ddd = ddd, gap_ddd = mean(reps["fhat", ]) - ddd)
}))
print(round(resC, 4))
cat(sprintf("gap_nee range: [%.4f, %.4f]; gap_ddd spread: %.4f (log 9! = %.4f)\n",
            min(resC[, "gap_nee"]), max(resC[, "gap_nee"]), diff(range(resC[, "gap_ddd"])), lgamma(10)))

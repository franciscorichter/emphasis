## H5 — does the thinning proposal density `q` (Model::sampling_prob, model.hpp:255-328)
## match the sampler, and is the IS weight f/q unbiased for the CR (and DD) likelihood?
##
## Part A  2-tip tree: derive the sampler's true density q_true(z) (n·lambda·mu·e^{-mu L},
##         parent uniform over |alive_ids|, compensator exp(-∫ n·lambda·(1-e^{-mu(T-t)}))),
##         check P(m=0), P(m=1) and the t* marginal against Monte Carlo, and show that
##         the code's logg differs from log q_true by exactly -Σ_j log K_j + Σ_j log max(1, A_j).
## Part B  2-tip tree: fhat(theta) vs the closed form 2·log p_1(T) (Kendall) on a theta grid.
## Part C  10-tip tree: fhat(theta) - Nee labelled-history log-likelihood, and - DDD::bd_loglik,
##         on a theta grid: must be a theta-independent constant.
## Part D  10-tip tree, linear DD (lambda = b0 + b1 N): fhat - DDD::dd_loglik(ddmodel=1) constant.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({ library(emphasis); library(DDD) })
set.seed(5)   # only affects R-side draws; the C++ RNG is clock-seeded

T_EXT_TIP <- 1e11; T_EXT_UNS <- 5e10
is_missing_node <- function(df) !(df$t_ext == T_EXT_TIP | df$t_ext == 0 | df$t_ext == T_EXT_UNS)

draw <- function(brts, pars8, Ns, model = c(0L, 0L, 0L), link = 0L) {
  emphasis:::augment_trees(brts = brts, pars = pars8, sample_size = as.integer(Ns),
                           maxN = as.integer(50 * Ns), max_missing = 10000L,
                           max_lambda = 1e6, num_threads = 1L,
                           model = as.integer(model), link = as.integer(link), rho = 1.0)
}
fhat_of <- function(d) emphasis:::.is_fhat(d$logf, d$logg, n_zero_weight = d$rejected_zero_weights)

## ---- closed forms (CR) --------------------------------------------------------------
p0 <- function(s, la, mu) { r <- la - mu
  if (abs(r) < 1e-12) la * s / (1 + la * s) else mu * (exp(r * s) - 1) / (la * exp(r * s) - mu) }
log_p1 <- function(Tt, la, mu) { r <- la - mu
  if (abs(r) < 1e-12) -2 * log(1 + la * Tt) else 2 * log(r) + r * Tt - 2 * log(la * exp(r * Tt) - mu) }
## Nee labelled-history log-likelihood in the code's f convention (one lambda per observed
## node, no combinatorial factor): (n-2) log la - ∫ n_obs(t) [(la+mu) - 2 la p0(T-t)] dt
nee_lab <- function(brts, la, mu) {
  Tt <- brts[1]; s <- sort(Tt - brts[-1]); bounds <- c(0, s, Tt); nseg <- length(bounds) - 1
  acc <- (length(brts) - 1) * log(la)   # n-2 observed nodes = length(brts)-1
  for (k in seq_len(nseg)) {
    nk <- 1 + k
    acc <- acc - nk * integrate(function(t) (la + mu) - 2 * la * sapply(t, p0, la = la, mu = mu),
                                bounds[k], bounds[k + 1], rel.tol = 1e-10)$value
  }
  acc
}

## ---- sampler's true density for a CR augmented tree -------------------------------
## I(a,b) = ∫_a^b (1 - e^{-mu (T-t)}) dt
Iseg <- function(a, b, Tt, mu) (b - a) - (exp(-mu * (Tt - b)) - exp(-mu * (Tt - a))) / mu
log_q_true <- function(df, la, mu) {
  Tt <- max(df$brts); o <- order(df$brts); df <- df[o, ]
  # compensator over segments (node.n = lineages alive on the segment ending at the node)
  prev <- 0; comp <- 0
  for (i in seq_len(nrow(df))) { comp <- comp + df$n[i] * la * Iseg(prev, df$brts[i], Tt, mu); prev <- df$brts[i] }
  lg <- -comp
  mis <- which(is_missing_node(df))
  for (j in mis) {
    ts <- df$brts[j]; te <- df$t_ext[j]
    A <- sum(df$t_ext != 0 & df$brts < ts & df$t_ext > ts)   # alive_ids, augment_tree.cpp:154-166
    lg <- lg + log(df$n[j] * la * mu) - mu * (te - ts) - log(max(1, A))
  }
  lg
}
## the code's charge, recomputed in R from model.hpp:299-309, to confirm the reading
log_q_code <- function(df, la, mu) {
  Tt <- max(df$brts); o <- order(df$brts); df <- df[o, ]
  prev <- 0; comp <- 0
  for (i in seq_len(nrow(df))) { comp <- comp + df$n[i] * la * Iseg(prev, df$brts[i], Tt, mu); prev <- df$brts[i] }
  lg <- -comp; tips <- df$n[1]; Ne <- 0
  for (i in seq_len(nrow(df))) {
    tips <- tips + (df$t_ext[i] == T_EXT_TIP); Ne <- Ne - (df$t_ext[i] == 0)
    if (is_missing_node(df)[i]) {
      lg <- lg + log(df$n[i] * la * mu) - mu * (df$t_ext[i] - df$brts[i]) - log(2 * tips + Ne); Ne <- Ne + 1
    }
  }
  lg
}

cat("==================== PART A: 2-tip tree, is q the sampler's density? ====================\n")
Tt <- 1; la <- 1; mu <- 0.6; Ns <- 20000
d <- draw(Tt, c(la, 0, 0, 0, mu, 0, 0, 0), Ns)
m <- sapply(d$trees, function(df) sum(is_missing_node(df)))
cat(sprintf("rejected_zero_weights=%d overruns=%d lambda=%d other=%d\n",
            d$rejected_zero_weights, d$rejected_overruns, d$rejected_lambda, d$rejected))
cat("m distribution:\n"); print(table(m))
P0_an <- exp(-2 * la * Iseg(0, Tt, Tt, mu))
q1 <- function(ts, te) 2 * la * mu * exp(-mu * (te - ts)) * exp(-la * (2 * Iseg(0, Tt, Tt, mu) + Iseg(ts, te, Tt, mu)))
q1_ts <- function(ts) sapply(ts, function(a) integrate(function(b) q1(a, b), a, Tt)$value)
P1_an <- integrate(q1_ts, 0, Tt)$value
se <- function(p) sqrt(p * (1 - p) / Ns)
cat(sprintf("P(m=0): analytic %.5f  MC %.5f (se %.5f)   z=%.2f\n", P0_an, mean(m == 0), se(P0_an), (mean(m == 0) - P0_an) / se(P0_an)))
cat(sprintf("P(m=1): analytic %.5f  MC %.5f (se %.5f)   z=%.2f\n", P1_an, mean(m == 1), se(P1_an), (mean(m == 1) - P1_an) / se(P1_an)))
## code's charge: for m=0 trees logg = -compensator = log P0 exactly; for m=1 trees Σ exp(logg) should be P1/4
cat(sprintf("m=0 trees: exp(logg_code) = %.6f (all equal: %s) vs P0 analytic %.6f\n",
            exp(d$logg[m == 0][1]), isTRUE(all.equal(sd(d$logg[m == 0]), 0)), P0_an))
one <- which(m == 1)
lqt <- sapply(one, function(i) log_q_true(d$trees[[i]], la, mu))
lqc <- sapply(one, function(i) log_q_code(d$trees[[i]], la, mu))
cat(sprintf("m=1 trees: max|logg_code - my R re-implementation of the code| = %.2e\n", max(abs(d$logg[one] - lqc))))
cat(sprintf("m=1 trees: logg_code - log q_true: range [%.6f, %.6f]  (-log 4 = %.6f)\n",
            min(d$logg[one] - lqt), max(d$logg[one] - lqt), -log(4)))
## t* marginal of the sampler vs analytic (KS)
tstar <- sapply(one, function(i) { df <- d$trees[[i]]; df$brts[is_missing_node(df)] })
cdf_ts <- function(x) sapply(x, function(u) integrate(q1_ts, 0, u)$value / P1_an)
ks <- suppressWarnings(ks.test(tstar, cdf_ts))
cat(sprintf("m=1 trees: KS test of t* against the q_true marginal: D=%.4f p=%.3f (n=%d)\n", ks$statistic, ks$p.value, length(one)))
## general identity for all trees: logg_code = log q_true - Σ log K_j + Σ log max(1, A_j)
idx <- which(m >= 2)[1:min(200, sum(m >= 2))]
gap <- sapply(idx, function(i) {
  df <- d$trees[[i]]; o <- order(df$brts); df <- df[o, ]
  tips <- df$n[1]; Ne <- 0; s <- 0
  for (k in seq_len(nrow(df))) {
    tips <- tips + (df$t_ext[k] == T_EXT_TIP); Ne <- Ne - (df$t_ext[k] == 0)
    if (is_missing_node(df)[k]) {
      A <- sum(df$t_ext != 0 & df$brts < df$brts[k] & df$t_ext > df$brts[k])
      s <- s - log(2 * tips + Ne) + log(max(1, A)); Ne <- Ne + 1
    }
  }
  (d$logg[i] - log_q_true(df, la, mu)) - s
})
cat(sprintf("m>=2 trees (%d checked): max |(logg_code - log q_true) - (-Σ log K_j + Σ log max(1,A_j))| = %.2e\n", length(idx), max(abs(gap))))

cat("\n==================== PART B: 2-tip tree, fhat vs 2 log p1(T) ====================\n")
grid <- rbind(c(1, 1, 0.1), c(1, 1, 0.5), c(1, 1, 1.0), c(1, 1, 1.4), c(1, 2, 1.8), c(3, 1, 0.8), c(3, 0.6, 0.5))
colnames(grid) <- c("T", "lambda", "mu")
resB <- t(apply(grid, 1, function(g) {
  Tt <- g[1]; la <- g[2]; mu <- g[3]
  reps <- replicate(4, { d <- draw(Tt, c(la, 0, 0, 0, mu, 0, 0, 0), 20000); c(fhat_of(d), mean(sapply(d$trees, function(df) sum(is_missing_node(df)))), d$rejected_zero_weights) })
  ex <- log_p1(Tt, la, mu)
  c(T = Tt, lambda = la, mu = mu, exact = ex, fhat = mean(reps[1, ]), sd_over_reps = sd(reps[1, ]),
    gap = mean(reps[1, ]) - ex, mean_m = mean(reps[2, ]), zero_w = sum(reps[3, ]))
}))
print(round(resB, 4))

cat("\n==================== PART C: 10-tip CR tree, fhat - Nee and fhat - DDD::bd_loglik ====================\n")
brts10 <- c(10, 7.3, 6.1, 4.8, 3.2, 2.5, 1.7, 0.9, 0.4)   # 10 tips, crown age 10
gridC <- rbind(c(0.2, 0.05), c(0.2, 0.15), c(0.3, 0.25), c(0.15, 0.1), c(0.4, 0.35), c(0.25, 0.28))
resC <- t(apply(gridC, 1, function(g) {
  la <- g[1]; mu <- g[2]
  reps <- replicate(3, { d <- draw(brts10, c(la, 0, 0, 0, mu, 0, 0, 0), 20000); c(fhat_of(d), emphasis:::.ess_from_lw(d$logf - d$logg), mean(sapply(d$trees, function(df) sum(is_missing_node(df))))) })
  nee <- nee_lab(brts10, la, mu)
  ddd <- DDD::bd_loglik(pars1 = c(la, mu, 0, 0), pars2 = c(0, 0, 0, 0, 2), brts = brts10, missnumspec = 0)
  c(lambda = la, mu = mu, fhat = mean(reps[1, ]), sd_reps = sd(reps[1, ]), ess = mean(reps[2, ]), mean_m = mean(reps[3, ]),
    nee_lab = nee, gap_nee = mean(reps[1, ]) - nee, ddd = ddd, gap_ddd = mean(reps[1, ]) - ddd)
}))
print(round(resC, 4))
cat(sprintf("gap_ddd spread across theta: %.4f ; log((n-1)!) = %.4f ; nee_lab - ddd spread: %.2e\n",
            diff(range(resC[, "gap_ddd"])), lgamma(10), diff(range(resC[, "nee_lab"] - resC[, "ddd"]))))

cat("\n==================== PART D: 10-tip linear-DD tree, fhat - DDD::dd_loglik(ddmodel=1) ====================\n")
## emphasis linear: lambda(N) = max(0, b0 + b1 N), mu = g0.   DDD ddmodel=1: lambda(N) = la0 - (la0-mu) N/K
gridD <- rbind(c(0.5, -0.02, 0.1), c(0.5, -0.03, 0.2), c(0.6, -0.025, 0.05), c(0.4, -0.015, 0.15))
resD <- t(apply(gridD, 1, function(g) {
  b0 <- g[1]; b1 <- g[2]; g0 <- g[3]; K <- (b0 - g0) / (-b1)
  reps <- replicate(3, { d <- draw(brts10, c(b0, b1, 0, 0, g0, 0, 0, 0), 20000, model = c(1L, 0L, 0L)); c(fhat_of(d), emphasis:::.ess_from_lw(d$logf - d$logg), d$rejected_zero_weights) })
  ddd <- DDD::dd_loglik(pars1 = c(b0, g0, K), pars2 = c(500, 1, 0, 0, 0, 2), brts = brts10, missnumspec = 0)
  c(b0 = b0, b1 = b1, g0 = g0, K = K, fhat = mean(reps[1, ]), sd_reps = sd(reps[1, ]), ess = mean(reps[2, ]), zero_w = sum(reps[3, ]), ddd = ddd, gap_ddd = mean(reps[1, ]) - ddd)
}))
print(round(resD, 4))
cat(sprintf("DD gap_ddd spread across theta: %.4f ; log((n-1)!) = %.4f\n", diff(range(resD[, "gap_ddd"])), lgamma(10)))

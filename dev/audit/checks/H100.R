## H100: "No test passes rho; C++ sampling terms (model.hpp:225-228, 251, 282, 295,
## 307-308, 456) are untested."  Self-contained.  Run: Rscript dev/audit/checks/H100.R
## Parts:
##  0. fact: does any test pass rho?
##  1. reference: Nee/Stadler crown likelihood with Bernoulli sampling rho, validated
##     (a) at rho=1 against DDD::bd_loglik and (b) at rho<1 against a brute-force forward
##     simulation (own recursive birth-death simulator with Bernoulli tip sampling; no formula involved).
##  2. the untested C++ terms: thinning-IS fhat at fixed theta for rho in {1, .7, .4}
##     vs reference + C, with C fixed at rho=1 (BDI, exact).  The rho terms are correct
##     iff gap(rho<1) == gap(rho=1) within replicate SE.
##  3. pin the cited lines: recompute logg (sampling_prob) by hand from returned trees;
##     logf(rho) - logf(1) on the same trees vs n_obs log rho + n_unsamp log(1-rho).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({ library(emphasis); library(ape); library(DDD) })
options(width = 130)
t0 <- Sys.time()

## ---------------------------------------------------------------- 0. fact
tf <- list.files("/Users/pancho/Code/emphasis/tests", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
hits <- sapply(tf, function(f) sum(grepl("rho", readLines(f, warn = FALSE))))
cat(sprintf("0. test files: %d; files mentioning 'rho': %d (%s)\n", length(tf), sum(hits > 0),
            paste(basename(names(hits)[hits > 0]), collapse = ",")))

## ---------------------------------------------------------------- reference
## p1_rho(t): density that a lineage alive t before present leaves exactly one sampled tip.
## Crown (labelled-history, up to a theta/rho-independent constant):
##   log L = 2 log p1(tc) + sum_{i>=2} [log lam + log p1(x_i)]
p1_rho <- function(t, lam, mu, rho) {
  r <- lam - mu; E <- exp(-r * t); den <- rho * lam + (lam * (1 - rho) - mu) * E
  rho * r^2 * E / den^2
}
nee_rho <- function(lam, mu, rho, brts)
  2 * log(p1_rho(brts[1], lam, mu, rho)) + sum(log(lam) + log(p1_rho(brts[-1], lam, mu, rho)))

## ---------------------------------------------------------------- 1a. rho = 1 vs DDD
set.seed(5)
phy  <- ape::rphylo(10, birth = 0.5, death = 0.3)
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
n    <- length(brts) + 1L
cat(sprintf("Tree: n_tips = %d, crown age = %.4f\n", n, brts[1]))
grid <- rbind(c(0.5, 0.3), c(0.4, 0.1), c(0.7, 0.5), c(0.25, 0.05))
off <- apply(grid, 1, function(p) nee_rho(p[1], p[2], 1, brts) -
               DDD::bd_loglik(c(p[1], p[2], 0, 0), c(0, 0, 1, 0, 2), brts, 0))
cat(sprintf("1a. nee_rho(rho=1) - DDD::bd_loglik(cond=0,btorph=1,soc=2) over 4 theta: spread = %.2e (constant => same likelihood)\n",
            diff(range(off))))

## ---------------------------------------------------------------- 1b. brute force at rho < 1
## Own recursive simulator (no TreeSim, no formula): a lineage at age a draws Exp(lam+mu);
## if it reaches the present it is sampled with prob rho; a speciation spawns two children.
## Returns k = number of sampled descendants and, for k >= 2, the age of their MRCA.
## Crown = two independent lineages at age Tc; keep draws with k_A >= 1 and k_B >= 1.
## Reference predictions:  P(n=3)/P(n=2) = 2 lam * int_0^T p1_rho(x) dx ;  x | n=3 ~ p1_rho(x)/int p1_rho
sim_lineage <- function(a, lam, mu, rho) {
  tau <- rexp(1, lam + mu)
  if (tau >= a) return(list(k = as.integer(runif(1) < rho), mrca = NA_real_))
  a2 <- a - tau
  if (runif(1) < mu / (lam + mu)) return(list(k = 0L, mrca = NA_real_))
  c1 <- sim_lineage(a2, lam, mu, rho); c2 <- sim_lineage(a2, lam, mu, rho)
  k <- c1$k + c2$k
  mrca <- if (c1$k > 0 && c2$k > 0) a2 else if (c1$k > 0) c1$mrca else c2$mrca
  list(k = k, mrca = mrca)
}
brute <- function(lam, mu, rho, Tc, nsim) {
  n_obs <- integer(0); x3 <- numeric(0)
  for (i in seq_len(nsim)) {
    A <- sim_lineage(Tc, lam, mu, rho); B <- sim_lineage(Tc, lam, mu, rho)
    if (A$k < 1 || B$k < 1) next
    k <- A$k + B$k; n_obs <- c(n_obs, k)
    if (k == 3) x3 <- c(x3, if (A$k == 2) A$mrca else B$mrca)
  }
  list(n_obs = n_obs, x3 = x3)
}
lamB <- 0.6; muB <- 0.3; TcB <- 2.5
for (rhoB in c(1, 0.5)) {
  set.seed(100 + rhoB * 10)
  b <- brute(lamB, muB, rhoB, TcB, 40000)
  k2 <- sum(b$n_obs == 2); k3 <- sum(b$n_obs == 3)
  I <- integrate(p1_rho, 0, TcB, lam = lamB, mu = muB, rho = rhoB)$value
  ratio_ref <- 2 * lamB * I
  ratio_obs <- k3 / k2; se <- ratio_obs * sqrt(1 / k2 + 1 / k3)
  Fx <- function(x) sapply(x, function(u) integrate(p1_rho, 0, u, lam = lamB, mu = muB, rho = rhoB)$value / I)
  ks <- suppressWarnings(ks.test(b$x3, Fx))
  ex <- integrate(function(x) x * p1_rho(x, lamB, muB, rhoB), 0, TcB)$value / I
  cat(sprintf("1b. rho=%.1f: crown-sampled draws %d/40000; #n=2: %d, #n=3: %d; P(3)/P(2) obs = %.4f (se %.4f) vs ref 2*lam*int p1_rho = %.4f  [z = %.2f]\n",
              rhoB, length(b$n_obs), k2, k3, ratio_obs, se, ratio_ref, (ratio_obs - ratio_ref) / se))
  cat(sprintf("    x|n=3: mean obs %.4f (se %.4f) vs ref %.4f [z = %.2f]; KS p = %.3f\n",
              mean(b$x3), sd(b$x3) / sqrt(length(b$x3)), ex, (mean(b$x3) - ex) / (sd(b$x3) / sqrt(length(b$x3))), ks$p.value))
}

## ---------------------------------------------------------------- 2. thinning IS vs reference
pars8 <- function(l, m) c(l, 0, 0, 0, m, 0, 0, 0)
aug <- function(brts, l, m, N, rho, maxN = 20L * N)
  emphasis:::augment_trees(brts, pars8(l, m), as.integer(N), as.integer(maxN), 100000L, 1e6, 1L,
                           c(0L, 0L, 0L), 0L, rho)
fhat_of <- function(r) { lw <- r$logf - r$logg; m <- max(lw)
  log(sum(exp(lw - m))) + m - log(length(lw) + r$rejected_zero_weights) }
ess_of <- function(r) { w <- exp(r$logf - r$logg - max(r$logf - r$logg)); sum(w)^2 / sum(w^2) }
## C from BDI at rho = 1 (exact, zero-variance under CR; verified in H2/map 4)
C <- mean(apply(grid, 1, function(p)
  emphasis:::.augment_tree_bdi(brts, pars = p, model_bin = c(0L, 0L, 0L), sample_size = 30L,
                               link = 0L, rho = 1.0)$fhat - nee_rho(p[1], p[2], 1, brts)))
cat(sprintf("\n2. C (BDI fhat - nee at rho=1) = %.6f\n", C))
N <- 1000L; reps <- 5L
cat(sprintf("   thinning augment_trees: N = %d, reps = %d; gap = fhat - (nee_rho + C); all gaps should share one value (0 if the thinning weight is exact)\n", N, reps))
res <- list()
for (rho in c(1, 0.7, 0.4)) for (i in seq_len(nrow(grid))) {
  p <- grid[i, ]
  g <- t(sapply(seq_len(reps), function(k) {
    r <- aug(brts, p[1], p[2], N, rho)
    nun <- mean(sapply(r$trees, function(tr) sum(tr$t_ext == 5e10)))
    c(fhat = fhat_of(r), ess = ess_of(r), zero = r$rejected_zero_weights, n_unsamp = nun)
  }))
  res[[length(res) + 1]] <- data.frame(rho = rho, lam = p[1], mu = p[2],
    gap = mean(g[, "fhat"]) - (nee_rho(p[1], p[2], rho, brts) + C),
    se = sd(g[, "fhat"]) / sqrt(reps), ess = mean(g[, "ess"]), zero_w = mean(g[, "zero"]),
    unsamp_per_tree = mean(g[, "n_unsamp"]))
}
res <- do.call(rbind, res)
print(res, digits = 4, row.names = FALSE)
g1 <- res[res$rho == 1, ]; g0 <- res[res$rho < 1, ]
cat(sprintf("   gap at rho=1: mean %.4f (range %.4f..%.4f) | gap at rho<1: mean %.4f (range %.4f..%.4f)\n",
            mean(g1$gap), min(g1$gap), max(g1$gap), mean(g0$gap), min(g0$gap), max(g0$gap)))
d <- merge(g0, g1, by = c("lam", "mu"), suffixes = c("", "_r1"))
d$z <- (d$gap - d$gap_r1) / sqrt(d$se^2 + d$se_r1^2)
cat("   (gap(rho<1) - gap(rho=1)) / se  per cell:\n")
print(d[, c("rho", "lam", "mu", "gap", "gap_r1", "z")], digits = 3, row.names = FALSE)
cat(sprintf("   max |z| = %.2f  (< 3 => rho terms consistent with the reference to IS precision)\n", max(abs(d$z))))
## rho-dependence signature the H2 bug would give: n*log(rho); here for comparison
cat(sprintf("   for scale: n*log(0.7) = %.3f, n*log(0.4) = %.3f (the size of the BDI error in H2)\n", n * log(0.7), n * log(0.4)))

## ---------------------------------------------------------------- 2b. does the residual gap shrink with N?
cat("\n2b. same gap with N = 5000 (finite-N IS bias shrinks; a weight-density mismatch does not)\n")
res2 <- list()
for (rho in c(1, 0.4)) for (p in list(c(0.5, 0.3), c(0.7, 0.5))) {
  g <- sapply(seq_len(5L), function(k) { r <- aug(brts, p[1], p[2], 5000L, rho); c(fhat_of(r), ess_of(r)) })
  res2[[length(res2) + 1]] <- data.frame(rho = rho, lam = p[1], mu = p[2],
    gap_N5000 = mean(g[1, ]) - (nee_rho(p[1], p[2], rho, brts) + C), se = sd(g[1, ]) / sqrt(5), ess = mean(g[2, ]),
    gap_N1000 = res$gap[res$rho == rho & res$lam == p[1] & res$mu == p[2]])
}
print(do.call(rbind, res2), digits = 4, row.names = FALSE)

## ---------------------------------------------------------------- 3. pin the cited lines
rho <- 0.5; lam <- 0.5; mu <- 0.3
r <- aug(brts, lam, mu, 200L, rho)
## 3a. logf sampling term (model.hpp:456-466): same trees scored at rho and at 1
e_r <- emphasis:::eval_logf(pars8(lam, mu), r$trees, model = c(0L, 0L, 0L), link = 0L, rho = rho)
e_1 <- emphasis:::eval_logf(pars8(lam, mu), r$trees, model = c(0L, 0L, 0L), link = 0L, rho = 1.0)
pred <- sapply(r$trees, function(tr) n * log(rho) + sum(tr$t_ext == 5e10) * log(1 - rho))
cat(sprintf("\n3a. logf(rho) - logf(1) vs n_obs*log(rho) + n_unsamp*log(1-rho): max |diff| = %.2e over %d trees; unsampled per tree: %s\n",
            max(abs((e_r$logf - e_1$logf) - pred)), length(pred),
            paste(range(sapply(r$trees, function(tr) sum(tr$t_ext == 5e10))), collapse = "..")))
## 3b. logg recomputed by hand from the map's formula (model.hpp:255-328) under CR:
##  log q = sum_missing [log(n mu lam) - mu(t_ext - s) - log(2 tips + Ne)]
##        + sum_unsamp  [log(n lam (1-rho)) - mu(T - s) - log(2 tips + Ne)]
##        - sum_i n_i lam [dt_i - (rho/mu)(e^{-mu(T-s_i)} - e^{-mu(T-s_{i-1})})]
hand_logg <- function(tr) {
  o <- order(tr$brts); tr <- tr[o, ]
  Tt <- tr$brts[nrow(tr)]; prev <- 0; inte <- 0; logg <- 0; tips <- tr$n[1]; Ne <- 0
  for (i in seq_len(nrow(tr))) {
    s <- tr$brts[i]; dt <- s - prev
    inte <- inte + tr$n[i] * lam * (dt - (rho / mu) * (exp(-mu * (Tt - s)) - exp(-mu * (Tt - prev))))
    is_ext <- tr$t_ext[i] == 0                 # extinction node: t_ext_extinct = 0 (model_helpers.hpp:49)
    is_uns <- tr$t_ext[i] == 5e10              # t_ext_unsampled (line 50)
    is_tip <- tr$t_ext[i] == 1e11              # t_ext_tip (line 48)
    is_mis <- !(is_ext || is_uns || is_tip)    # line 108
    tips <- tips + is_tip; Ne <- Ne - is_ext
    if (is_mis) { logg <- logg + log(tr$n[i] * mu * lam) - mu * (tr$t_ext[i] - s) - log(2 * tips + Ne); Ne <- Ne + 1 }
    else if (is_uns) { logg <- logg + log(tr$n[i] * lam * (1 - rho)) - mu * (Tt - s) - log(2 * tips + Ne); Ne <- Ne + 1 }
    prev <- s
  }
  logg - inte
}
hl <- sapply(r$trees, hand_logg)
cat(sprintf("3b. hand logg vs returned logg: max |diff| = %.2e over %d trees (formula as mapped == code)\n",
            max(abs(hl - r$logg)), length(hl)))
cat(sprintf("    sentinel/columns seen: t_ext values in {%s}; names: %s\n",
            paste(head(sort(unique(unlist(lapply(r$trees, function(tr) tr$t_ext[tr$t_ext >= 1e10])))), 3), collapse = ","),
            paste(names(r$trees[[1]]), collapse = ",")))
cat(sprintf("\nelapsed %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

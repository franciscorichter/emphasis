## H14 independent replication.
## Question: is sqrt(loglik_var) the MC s.e. of the reported loglik (LME with n_zero),
## or of another statistic (K=2 moment-corrected mean(lw) + log(1+m2/2))?
## Vary what the verifier did not: two other trees (20 & 40 tips, other seeds), N in {50, 200},
## R = 30 fresh E-steps per case so ground truth is replicate-averaged, and the BDI sampler for dd.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(TreeSim); library(ape)})
options(width = 160)

lme  <- function(lw, nz = 0L) { m <- max(lw); log(sum(exp(lw - m))) + m - log(length(lw) + nz) }
bcK2 <- function(lw) { e <- mean(lw); e + log(1 + mean((lw - e)^2) / 2) }
ess  <- function(lw) { w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }
bootv <- function(lw, f, B = 200L) { n <- length(lw)
  stats::var(vapply(seq_len(B), function(.) f(lw[sample.int(n, n, TRUE)]), 0)) }

## sanity on one vector: package internal == bcK2, != lme
set.seed(1); lw0 <- rnorm(100, 0, 1.5)
cat(sprintf("sanity: .is_fhat(bias_correct=TRUE)=%.6f  bcK2=%.6f  lme=%.6f  .is_fhat(default)=%.6f\n",
            emphasis:::.is_fhat(lw0, rep(0,100), bias_correct = TRUE), bcK2(lw0), lme(lw0),
            emphasis:::.is_fhat(lw0, rep(0,100))))
## n_zero invisible to bootstrap? (same lw, n_zero irrelevant to .bootstrap_fhat_var by signature)
cat("   .bootstrap_fhat_var formals:", paste(names(formals(emphasis:::.bootstrap_fhat_var)), collapse = ","), "\n")

run_case <- function(lab, brts, p8, mb, N, R = 30L) {
  out <- t(vapply(seq_len(R), function(r) {
    a <- emphasis:::augment_trees(brts, as.numeric(p8), as.integer(N), 20000L, 10000L, 1e6, 1L,
                                  as.integer(mb), 0L, 1.0)
    if (r == 1L && lab == "first") cat("augment_trees fields:", paste(names(a), collapse = ","), "\n")
    lw <- a$logf - a$logg
    nz <- if (!is.null(a$rejected_zero_weights)) as.integer(a$rejected_zero_weights) else 0L
    c(lme = lme(lw, nz), bc = bcK2(lw),
      se_pkg = sqrt(emphasis:::.bootstrap_fhat_var(a$logf, a$logg, K = 2L, B = 200L)),
      se_bootlme = sqrt(bootv(lw, lme)), ess = ess(lw), nz = nz, n = length(lw))
  }, numeric(7)))
  data.frame(case = lab, N = N, n_valid = mean(out[, "n"]), ESS = mean(out[, "ess"]), nz = mean(out[, "nz"]),
             ell_LME = mean(out[, "lme"]), ell_BC = mean(out[, "bc"]), gap = mean(out[, "lme"] - out[, "bc"]),
             se_pkg = mean(out[, "se_pkg"]),
             sd_BC_true = sd(out[, "bc"]), sd_LME_true = sd(out[, "lme"]),
             se_bootLME = mean(out[, "se_bootlme"]),
             ratio_pkg_over_true = mean(out[, "se_pkg"]) / sd(out[, "lme"]),
             ratio_pkg_over_trueBC = mean(out[, "se_pkg"]) / sd(out[, "bc"]))
}

res <- list()
for (tr in list(list(seed = 11, n = 20, lam = 0.6, mu = 0.3), list(seed = 22, n = 40, lam = 0.4, mu = 0.1))) {
  set.seed(tr$seed)
  phy  <- TreeSim::sim.bd.taxa(tr$n, 1, tr$lam, tr$mu, complete = FALSE)[[1]]
  brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
  cat(sprintf("\n=== tree: %d tips, crown age %.2f, sim lam=%.1f mu=%.1f ===\n", tr$n, brts[1], tr$lam, tr$mu))
  for (N in c(50L, 200L)) {
    res[[length(res) + 1]] <- run_case(sprintf("cr(%.1f,%.1f) %dt", tr$lam, tr$mu, tr$n), brts,
                                       c(tr$lam, 0,0,0, tr$mu, 0,0,0), c(0,0,0), N)
    res[[length(res) + 1]] <- run_case(sprintf("dd(b0=%.1f,bN=-.01,mu=%.1f) %dt", tr$lam + 0.2, tr$mu, tr$n), brts,
                                       c(tr$lam + 0.2, -0.01, 0,0, tr$mu, 0,0,0), c(1,0,0), N)
    res[[length(res) + 1]] <- run_case(sprintf("cr far (%.1f,%.1f) %dt", 1.2, 1.0, tr$n), brts,
                                       c(1.2, 0,0,0, 1.0, 0,0,0), c(0,0,0), N)
  }
}
tab <- do.call(rbind, res)
print(tab, digits = 3, row.names = FALSE)

## BDI default sampler, dd, on the 20-tip tree: what estimate_rates actually uses for dd/cr at rho=1
set.seed(11)
phy  <- TreeSim::sim.bd.taxa(20, 1, 0.6, 0.3, complete = FALSE)[[1]]
brts <- sort(as.numeric(ape::branching.times(phy)), decreasing = TRUE)
cat("\n=== BDI sampler (default path), dd b0=.8 bN=-.01 mu=.3, 20-tip tree, N=100, R=15 fresh E-steps ===\n")
t0 <- proc.time()[3]
bd <- t(vapply(seq_len(15L), function(r) {
  e <- emphasis:::.augment_tree_bdi(tree = brts, pars = c(0.8, -0.01, 0, 0, 0.3, 0, 0, 0),
                                    model_bin = c(1L,0L,0L), sample_size = 100L, max_missing = 10000L,
                                    link = 0L, rho = 1.0)
  lw <- e$logf - e$logg
  c(lme = lme(lw), bc = bcK2(lw), fhat = e$fhat,
    se_pkg = sqrt(emphasis:::.bootstrap_fhat_var(e$logf, e$logg, K = 2L, B = 200L)),
    se_bootlme = sqrt(bootv(lw, lme)), ess = ess(lw), n = length(lw))
}, numeric(7)))
cat(sprintf("   %.0f s; n_valid=%.0f ESS=%.1f  ell_LME(=e$fhat? max|diff|=%.2e)=%.3f  ell_BC=%.3f  gap=%.3f\n",
            proc.time()[3] - t0, mean(bd[, "n"]), mean(bd[, "ess"]), max(abs(bd[, "lme"] - bd[, "fhat"])),
            mean(bd[, "lme"]), mean(bd[, "bc"]), mean(bd[, "lme"] - bd[, "bc"])))
cat(sprintf("   mean se_pkg=%.4f  true sd(ell_BC)=%.4f  true sd(ell_LME)=%.4f  mean se_bootLME=%.4f  ratio pkg/trueLME=%.2f\n",
            mean(bd[, "se_pkg"]), sd(bd[, "bc"]), sd(bd[, "lme"]), mean(bd[, "se_bootlme"]),
            mean(bd[, "se_pkg"]) / sd(bd[, "lme"])))

## compare_models arithmetic with replicate-averaged s.e. (dd 20-tip N=200 thinning row and BDI row)
pv <- function(dA, v1, v2) 2 * pnorm(-abs(dA) / sqrt(4 * (v1 + v2)))
r <- tab[tab$case == "dd(b0=0.8,bN=-.01,mu=0.3) 20t" & tab$N == 200, ]
for (dA in c(1, 2)) cat(sprintf("   thinning dd 20t N=200, dAIC=%d: p(pkg se)=%.3f  p(true LME sd)=%.3f  p(true BC sd)=%.3f\n",
                                dA, pv(dA, r$se_pkg^2, r$se_pkg^2), pv(dA, r$sd_LME_true^2, r$sd_LME_true^2),
                                pv(dA, r$sd_BC_true^2, r$sd_BC_true^2)))
for (dA in c(1, 2)) cat(sprintf("   BDI dd 20t N=100,       dAIC=%d: p(pkg se)=%.3f  p(true LME sd)=%.3f\n",
                                dA, pv(dA, mean(bd[, "se_pkg"])^2, mean(bd[, "se_pkg"])^2),
                                pv(dA, sd(bd[, "lme"])^2, sd(bd[, "lme"])^2)))

# H26: CEM maxN=10 / num_trees=5 -> theta-dependent NA particles filter the search.
# Part A: NA fraction of augment_trees(ss=5,maxN=10) over a (lambda,mu) grid on one tree.
# Part B: CEM (cr, rho=1, cond=NULL) with maxN=10 vs maxN=200, replicated; NA-particle
#         parameter distribution in the final population; reference MLE from DDD::bd_ML.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages({library(emphasis); library(ape); library(DDD)})
msg <- function(...) cat(sprintf(...), file = stderr())
set.seed(26)
repeat { phy <- rlineage(0.5, 0.2, Tmax = 8); ex <- drop.fossil(phy)
         if (Ntip(ex) >= 15 && Ntip(ex) <= 30) break }
brts <- emphasis:::.extract_brts(ex)
msg("tree: ntips=%d crown=%.3f\n", Ntip(ex), max(brts))

## ---------- Part A: NA map --------------------------------------------------
one_call <- function(lam, mu, ss, maxN) {
  r <- tryCatch(emphasis:::augment_trees(brts, c(lam,0,0,0,mu,0,0,0), ss, maxN,
                                         1e4L, 1e6, 1L, c(0L,0L,0L), 0L, 1.0),
                error = function(e) e)
  if (inherits(r, "error")) {
    m <- conditionMessage(r)
    zw <- as.integer(sub(".*?(\\d+) zero weights.*", "\\1", m))
    ov <- as.integer(sub(".*?(\\d+) overruns.*", "\\1", m))
    c(na = 1, zw = zw, ov = ov)
  } else c(na = 0, zw = r$rejected_zero_weights, ov = r$rejected_overruns)
}
grid <- expand.grid(lam = c(0.2, 0.5, 1.0), ratio = c(0, 0.3, 0.6, 0.8, 0.95, 1.1))
grid$mu <- grid$lam * grid$ratio
R <- 10
msg("\nPart A: NA fraction over %d calls, ss=5, maxN=10 (then maxN=200)\n", R)
msg("%6s %6s %6s | %8s %8s %8s | %8s\n", "lam","mu","mu/lam","NA10","zw/call","ov/call","NA200")
resA <- t(apply(grid, 1, function(g) {
  o10  <- replicate(R, one_call(g["lam"], g["mu"], 5L, 10L))
  o200 <- replicate(R, one_call(g["lam"], g["mu"], 5L, 200L))
  out <- c(lam = g[["lam"]], mu = g[["mu"]], ratio = g[["ratio"]],
           NA10 = mean(o10["na",]), zw = mean(o10["zw",], na.rm=TRUE),
           ov = mean(o10["ov",], na.rm=TRUE), NA200 = mean(o200["na",]))
  msg("%6.2f %6.3f %6.2f | %8.2f %8.1f %8.1f | %8.2f\n",
      out["lam"], out["mu"], out["ratio"], out["NA10"], out["zw"], out["ov"], out["NA200"])
  out
}))

## ---------- Part B: CEM fits -----------------------------------------------
ref <- tryCatch(DDD::bd_ML(brts = brts, cond = 0, btorph = 0, soc = 2,
                           idparsopt = 1:2, initparsopt = c(0.5, 0.1),
                           verbose = FALSE), error = function(e) NULL)
if (!is.null(ref)) msg("\nDDD::bd_ML (cond=0): lambda=%.4f mu=%.4f loglik=%.4f\n",
                       ref$lambda0, ref$mu0, ref$loglik)

lb <- c(0, 0); ub <- c(1.0, 1.0)
fit_cem <- function(maxN) {
  set.seed(NULL)
  f <- estimate_rates(ex, method = "cem", model = "cr",
        control = list(lower_bound = lb, upper_bound = ub, maxN = maxN,
                       num_trees = 5, num_particles = 30, max_iter = 12,
                       max_time = 45, num_threads = 1))
  fp  <- f$details$final_pop
  na  <- is.na(fp$fhat)
  rat <- fp$pars[,5] / pmax(fp$pars[,1], 1e-9)
  c(maxN = maxN, lambda = unname(f$pars[1]), mu = unname(f$pars[2]), loglik = unname(f$loglik),
    iters = length(f$details$best_loglik),
    mean_valid = mean(f$details$history$n_valid),
    nNA_final = sum(na),
    ratio_NA = if (any(na)) mean(rat[na]) else NA,
    ratio_ok = mean(rat[!na]),
    mu_NA = if (any(na)) mean(fp$pars[na,5]) else NA,
    mu_ok = mean(fp$pars[!na,5]))
}
msg("\nPart B: CEM cr rho=1 cond=NULL, particles=30, num_trees=5, max_iter=12\n")
resB <- NULL
for (rep in 1:3) for (mN in c(10L, 200L)) {
  t0 <- proc.time()[3]
  r <- fit_cem(mN)
  msg("rep %d maxN=%3d: lambda=%.3f mu=%.3f loglik=%.3f iters=%d valid/iter=%.1f | final pop: nNA=%d  mu/lam NA=%.2f ok=%.2f  mu NA=%.3f ok=%.3f  (%.0fs)\n",
      rep, mN, r["lambda"], r["mu"], r["loglik"], r["iters"], r["mean_valid"],
      r["nNA_final"], r["ratio_NA"], r["ratio_ok"], r["mu_NA"], r["mu_ok"], proc.time()[3]-t0)
  resB <- rbind(resB, r)
}
msg("\nSummary by maxN (mean over reps):\n")
print(aggregate(resB[, c("lambda","mu","loglik","mean_valid","nNA_final")],
                by = list(maxN = resB[,"maxN"]), FUN = mean), file = stderr())

## ---------- Part C: dd model, linear link -- the zero-weight boundary ------
## lambda = max(0, b0 + bN * N).  K = -b0/bN.  When K is close to N_obs, an
## augmented missing lineage can push N to K, lambda -> 0, logf = -inf, and the
## tree is a zero-weight rejection (counted in maxN).  Map NA fraction vs K.
set.seed(7)
sdd  <- DDD::dd_sim(c(0.8, 0.1, 30), 8)
brts_dd <- sort(sdd$brts, decreasing = TRUE); Nobs <- length(brts_dd) + 1L
msg("
Part C: dd tree  Nobs=%d crown=%.2f (dd_sim lambda0=0.8 mu=0.1 K=30)\n", Nobs, max(brts_dd))
one_dd <- function(b0, bN, g0, ss, maxN, R) {
  o <- replicate(R, {
    r <- tryCatch(emphasis:::augment_trees(brts_dd, c(b0,bN,0,0,g0,0,0,0), ss, maxN,
                                           1e4L, 1e6, 1L, c(1L,0L,0L), 0L, 1.0),
                  error = function(e) e)
    if (inherits(r, "error")) c(NA_real_, 1, as.integer(sub(".*?(\\d+) zero weights.*", "\\1", conditionMessage(r))))
    else c(emphasis:::.is_fhat(r$logf, r$logg, n_zero_weight = r$rejected_zero_weights), 0, r$rejected_zero_weights)
  })
  c(fhat = mean(o[1,], na.rm = TRUE), na = mean(o[2,]), zw = mean(o[3,], na.rm = TRUE))
}
msg("%8s %8s | %7s %7s %9s | %7s %9s\n", "K","K-Nobs","NA10","zw/call","fhat10","NA200","fhat200")
resC <- NULL
for (K in c(Nobs - 1, Nobs, Nobs + 1, Nobs + 2, Nobs + 3, Nobs + 5, Nobs + 8, Nobs + 15, 2*Nobs, 4*Nobs)) {
  b0 <- 0.8; bN <- -b0 / K
  a <- one_dd(b0, bN, 0.1, 5L, 10L, 10); b <- one_dd(b0, bN, 0.1, 5L, 200L, 10)
  msg("%8.1f %8.1f | %7.2f %7.1f %9.3f | %7.2f %9.3f\n", K, K - Nobs, a["na"], a["zw"], a["fhat"], b["na"], b["fhat"])
  resC <- rbind(resC, c(K = K, NA10 = a["na"], zw10 = a["zw"], fhat10 = a["fhat"], NA200 = b["na"], fhat200 = b["fhat"]))
}

## Part D: CEM dd fits, maxN=10 vs 200 -- does the estimate of K move?
phy_dd <- DDD::L2phylo(sdd$L, dropextinct = TRUE)
lbd <- c(0.05, -0.1, 0, 0); ubd <- c(2, 0, 0.8, 0)   # (b0, bN, g0, gN); gN pinned at 0
fit_dd <- function(maxN) {
  f <- estimate_rates(phy_dd, method = "cem", model = "dd",
        control = list(lower_bound = lbd, upper_bound = ubd, maxN = maxN,
                       num_trees = 5, num_particles = 30, max_iter = 12,
                       max_time = 60, num_threads = 1))
  fp <- f$details$final_pop; na <- is.na(fp$fhat)
  Kp <- -fp$pars[,1] / pmin(fp$pars[,2], -1e-9)
  c(maxN = maxN, b0 = unname(f$pars[1]), bN = unname(f$pars[2]), g0 = unname(f$pars[3]),
    K = unname(-f$pars[1]/f$pars[2]), loglik = unname(f$loglik),
    mean_valid = mean(f$details$history$n_valid), nNA_final = sum(na),
    K_NA = if (any(na)) median(Kp[na]) else NA, K_ok = median(Kp[!na]))
}
msg("\nPart D: CEM dd rho=1 cond=NULL, particles=30, num_trees=5, max_iter=12, bounds b0[0.05,2] bN[-0.1,0] g0[0,0.8]\n")
resD <- NULL
for (rep in 1:3) for (mN in c(10L, 200L)) {
  t0 <- proc.time()[3]; r <- fit_dd(mN)
  msg("rep %d maxN=%3d: b0=%.3f bN=%.4f g0=%.3f K=%.1f loglik=%.3f valid/iter=%.1f | final pop nNA=%d  median K: NA=%.1f ok=%.1f (%.0fs)\n",
      rep, mN, r["b0"], r["bN"], r["g0"], r["K"], r["loglik"], r["mean_valid"], r["nNA_final"], r["K_NA"], r["K_ok"], proc.time()[3]-t0)
  resD <- rbind(resD, r)
}
msg("\nSummary dd by maxN (mean over reps):\n")
print(aggregate(resD[, c("b0","bN","g0","K","loglik","mean_valid","nNA_final")],
                by = list(maxN = resD[,"maxN"]), FUN = mean), file = stderr())
saveRDS(list(A = resA, B = resB, C = resC, D = resD, ref = ref, brts = brts, brts_dd = brts_dd),
        "/Users/pancho/Code/emphasis/dev/audit/checks/H26.rds")

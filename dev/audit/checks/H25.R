## H25: CEM final estimate = softmax(fhat_e)-weighted mean of the final elites,
## where each elite fhat is one fresh draw with max(20, sample_size) trees
## (R/de.R:793-849). Claims to test:
##  (i)  elite fhats re-drawn with 20 trees are noisy by units -> softmax weights
##       are near one-hot on a noise-selected winner;
##  (ii) on a ridge (lambda-mu for CR, intercept-slope for dd) the weighted mean
##       lies off the ridge, so the exact loglik at obtained_estim is below the
##       exact loglik at the best elite;
##  (iii) best_IS$fhat (at obtained_estim) differs from the best elite's fhat
##       re-evaluated with the same tree count.
## Method: run emphasis_cem at small settings, then REPLAY the final block
## from the returned final_pop (same code path: .simulate_particle + .is_fhat
## + softmax) R times, and score every elite and every replay's obtained_estim
## with the exact likelihood (DDD::bd_loglik for CR, DDD::dd_loglik for dd,
## cond = 0, btorph = 1, soc = 2; theta-independent constant offsets only).
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(DDD); library(ape)})
options(width = 140)
t0 <- proc.time()

softmax_mean <- function(fh, P) {
  ok <- is.finite(fh); w <- exp(fh[ok] - max(fh[ok])); w <- w / sum(w)
  list(est = as.numeric(colSums(P[ok, , drop = FALSE] * w)), w = w,
       entropy = -sum(w * log(pmax(w, 1e-300))), wmax = max(w),
       winner = which(ok)[which.max(w)])
}
eval_at <- function(pars8, brts, model, n_trees, maxN) {
  raw <- emphasis:::.simulate_particle(brts, pars8, model, 0L, sample_size = n_trees,
                                       maxN = maxN, max_missing = 1e4, max_lambda = 1e6,
                                       num_threads = 1L, rho = 1.0)
  if (is.null(raw) || length(raw$logf) == 0L) return(c(fhat = NA_real_, ess = NA_real_))
  lw <- raw$logf - raw$logg
  c(fhat = emphasis:::.is_fhat(raw$logf, raw$logg,
                               n_zero_weight = emphasis:::.n0(raw$rejected_zero_weights)),
    ess = emphasis:::.ess_from_lw(lw))
}

replay_final <- function(fit, brts, model, exact, R = 6L, n_trees = 20L) {
  fp <- fit$final_pop
  n_el <- max(1L, ceiling(0.5 * sum(!is.na(fp$fhat))))
  vix <- which(!is.na(fp$fhat)); ord <- vix[order(fp$fhat[vix], decreasing = TRUE)]
  eix <- ord[seq_len(min(n_el, length(ord)))]
  P <- as.matrix(fp$pars[eix, , drop = FALSE])
  maxN <- max(10L, n_trees * 5L)
  cat(sprintf("  final population: %d valid of %d; %d elites; search fhat range of elites = [%.2f, %.2f]\n",
              length(vix), nrow(fp$pars), length(eix), min(fp$fhat[eix]), max(fp$fhat[eix])))
  ex_el <- apply(P, 1, exact)
  cat(sprintf("  exact loglik over elites: min=%.3f max=%.3f spread=%.3f\n",
              min(ex_el), max(ex_el), max(ex_el) - min(ex_el)))
  # replay the elite re-evaluation R times
  FH <- matrix(NA_real_, R, nrow(P)); ES <- FH
  est <- matrix(NA_real_, R, ncol(P)); ent <- wmax <- numeric(R); win <- integer(R)
  ex_est <- ex_win <- numeric(R); fh_est <- ess_est <- numeric(R)
  for (r in seq_len(R)) {
    for (e in seq_len(nrow(P))) { v <- eval_at(P[e, ], brts, model, n_trees, maxN); FH[r, e] <- v["fhat"]; ES[r, e] <- v["ess"] }
    sm <- softmax_mean(FH[r, ], P)
    est[r, ] <- sm$est; ent[r] <- sm$entropy; wmax[r] <- sm$wmax; win[r] <- sm$winner
    ex_est[r] <- exact(est[r, ]); ex_win[r] <- ex_el[sm$winner]
    v <- eval_at(est[r, ], brts, model, n_trees, maxN); fh_est[r] <- v["fhat"]; ess_est[r] <- v["ess"]
  }
  noise_sd <- apply(FH, 2, sd, na.rm = TRUE)
  cat(sprintf("  20-tree fhat noise per elite: SD median=%.3f max=%.3f ; mean ESS/20 = %.1f\n",
              median(noise_sd, na.rm = TRUE), max(noise_sd, na.rm = TRUE), mean(ES, na.rm = TRUE)))
  cat(sprintf("  per-replay within-elite fhat range (max-min over elites): %s\n",
              paste(sprintf("%.2f", apply(FH, 1, function(x) diff(range(x, na.rm = TRUE)))), collapse = " ")))
  cat(sprintf("  softmax: entropy = %s  (log(n_el)=%.2f) ; max weight = %s\n",
              paste(sprintf("%.2f", ent), collapse = " "), log(nrow(P)),
              paste(sprintf("%.2f", wmax), collapse = " ")))
  cat(sprintf("  winner elite index per replay: %s  (exact-best elite = %d)\n",
              paste(win, collapse = " "), which.max(ex_el)))
  cat(sprintf("  exact loglik at obtained_estim per replay: %s\n", paste(sprintf("%.3f", ex_est), collapse = " ")))
  cat(sprintf("  exact loglik at winner elite         : %s\n", paste(sprintf("%.3f", ex_win), collapse = " ")))
  cat(sprintf("  exact(obtained_estim) - exact(best elite overall)   : %s\n",
              paste(sprintf("%+.3f", ex_est - max(ex_el)), collapse = " ")))
  cat(sprintf("  20-tree fhat at obtained_estim : %s\n", paste(sprintf("%.2f", fh_est), collapse = " ")))
  cat(sprintf("  fhat(obtained_estim) - exact(obtained_estim) [offset]: %s\n",
              paste(sprintf("%+.2f", fh_est - ex_est), collapse = " ")))
  cat(sprintf("  fhat(winner elite, 20 trees) - fhat(obtained_estim, 20 trees): %s\n",
              paste(sprintf("%+.2f", FH[cbind(seq_len(R), win)] - fh_est), collapse = " ")))
  cat("  obtained_estim per replay (SD across replays = ", paste(sprintf("%.4f", apply(est, 2, sd)), collapse = " "), "):\n")
  print(round(est, 4))
  invisible(list(P = P, FH = FH, est = est, ex_el = ex_el, ex_est = ex_est))
}

## ---------------------------------------------------------------- A. CR tree
cat("=== A. CR tree (30 tips, lambda=0.5, mu=0.3), emphasis_cem cr/linear ===\n")
# pick a tree whose unconditioned MLE has interior mu (so a lambda-mu ridge exists)
for (seed in 25:60) {
  set.seed(seed)
  repeat { tr <- ape::rlineage(0.5, 0.3, Tmax = 6); tr <- ape::drop.fossil(tr); if (ape::Ntip(tr) >= 25 && ape::Ntip(tr) <= 40) break }
  brts_cr <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
  ml <- suppressMessages(invisible(capture.output(mlv <- DDD::bd_ML(brts_cr, initparsopt = c(0.5, 0.3), idparsopt = 1:2,
                 parsfix = c(0, 0), idparsfix = 3:4, cond = 0, btorph = 1, soc = 2, verbose = FALSE))))
  ml <- mlv
  if (ml$mu0 > 0.1) break
}
cat("seed =", seed, " n_tips =", length(brts_cr) + 1, " crown age =", round(brts_cr[1], 3), "\n")
exact_cr <- function(p8) DDD::bd_loglik(c(p8[1], p8[5], 0, 0), c(0, 0, 1, 0, 2), brts_cr, 0)
lb <- c(0.05, 0, 0, 0, 0, 0, 0, 0); ub <- c(2, 0, 0, 0, 1.5, 0, 0, 0)
cat(sprintf("DDD::bd_ML(cond=0): lambda=%.4f mu=%.4f loglik=%.4f\n", ml$lambda0, ml$mu0, ml$loglik))
fit_cr <- emphasis:::emphasis_cem(brts_cr, max_iter = 25L, num_points = 30L, max_missing = 1e4,
                       sd_vec = (ub - lb) / 4, lower_bound = lb, upper_bound = ub,
                       maxN = 10L, sample_size = 1L, disc_prop = 0.5, num_threads = 1L,
                       model = c(0L, 0L, 0L), link = 0L, verbose = FALSE, max_time = 120)
cat(sprintf("emphasis_cem: converged=%s iters=%d obtained_estim=(%.4f, %.4f) best_IS$fhat=%.3f ESS=%.1f ; exact at estim=%.3f\n",
            fit_cr$converged, length(fit_cr$best_loglik), fit_cr$obtained_estim[1], fit_cr$obtained_estim[5],
            fit_cr$best_IS$fhat, fit_cr$best_IS$ESS, exact_cr(fit_cr$obtained_estim)))
A <- replay_final(fit_cr, brts_cr, c(0L, 0L, 0L), exact_cr, R = 6L)
cat(sprintf("  exact loglik at DDD MLE = %.3f ; exact-best elite is %.3f below it\n", ml$loglik, ml$loglik - max(A$ex_el)))
cat(sprintf("  elite (lambda, mu) cloud: lambda in [%.3f, %.3f], mu in [%.3f, %.3f]; cor(lambda, mu) = %.2f\n",
            min(A$P[, 1]), max(A$P[, 1]), min(A$P[, 5]), max(A$P[, 5]), cor(A$P[, 1], A$P[, 5])))
cat("  elapsed:", round((proc.time() - t0)[3]), "s\n")

## ---------------------------------------------------------------- B. dd tree
cat("\n=== B. DD tree (DDD::dd_sim lambda0=0.8 mu0=0.1 K=40, age 8), emphasis_cem dd/linear, gammaN fixed 0 ===\n")
set.seed(1)
s <- DDD::dd_sim(c(0.8, 0.1, 40), 8)
brts_dd <- sort(as.numeric(ape::branching.times(s$tes)), decreasing = TRUE)
cat("n_tips =", length(brts_dd) + 1, "\n")
exact_dd <- function(p8) {   # beta0=lambda0, betaN=-(lambda0-mu0)/K, gamma0=mu0  (ddmodel=1)
  l0 <- p8[1]; bN <- p8[2]; m0 <- p8[5]
  if (!(l0 > m0) || !(bN < 0)) return(-Inf)
  K <- (l0 - m0) / (-bN)
  DDD::dd_loglik(c(l0, m0, K), c(200, 1, 0, 1, 0, 2), brts_dd, 0)
}
lb <- c(0.1, -0.1, 0, 0, 0, 0, 0, 0); ub <- c(3, -1e-4, 0, 0, 0.8, 0, 0, 0)
fit_dd <- emphasis:::emphasis_cem(brts_dd, max_iter = 25L, num_points = 30L, max_missing = 1e4,
                       sd_vec = (ub - lb) / 4, lower_bound = lb, upper_bound = ub,
                       maxN = 10L, sample_size = 1L, disc_prop = 0.5, num_threads = 1L,
                       model = c(1L, 0L, 0L), link = 0L, verbose = FALSE, max_time = 150)
e <- fit_dd$obtained_estim
cat(sprintf("emphasis_cem: converged=%s iters=%d obtained_estim=(b0=%.4f, bN=%.5f, g0=%.4f) -> K=%.1f ; best_IS$fhat=%.3f ESS=%.1f ; exact at estim=%.3f\n",
            fit_dd$converged, length(fit_dd$best_loglik), e[1], e[2], e[5], (e[1] - e[5]) / (-e[2]),
            fit_dd$best_IS$fhat, fit_dd$best_IS$ESS, exact_dd(e)))
B <- replay_final(fit_dd, brts_dd, c(1L, 0L, 0L), exact_dd, R = 6L)
cat(sprintf("  elite cloud: b0 in [%.3f, %.3f], bN in [%.5f, %.5f], g0 in [%.3f, %.3f]; cor(b0, bN) = %.2f ; implied K in [%.1f, %.1f]\n",
            min(B$P[, 1]), max(B$P[, 1]), min(B$P[, 2]), max(B$P[, 2]), min(B$P[, 5]), max(B$P[, 5]),
            cor(B$P[, 1], B$P[, 2]), min((B$P[, 1] - B$P[, 5]) / -B$P[, 2]), max((B$P[, 1] - B$P[, 5]) / -B$P[, 2])))
cat("  elapsed:", round((proc.time() - t0)[3]), "s\n")

## ---------------------------------------------------------------- C. default control, run to convergence
cat("\n=== C. estimate_rates(method='cem') with default control (50 particles, 1 tree, max_iter 50, patience 5) ===\n")
for (case in c("cr", "dd")) {
  if (case == "cr") { brts <- brts_cr; lbc <- c(0.05, 0); ubc <- c(2, 1.5); model <- c(0L,0L,0L); exact <- exact_cr }
  else              { brts <- brts_dd; lbc <- c(0.1, -0.1, 0, 0); ubc <- c(3, -1e-4, 0.8, 0); model <- c(1L,0L,0L); exact <- exact_dd }
  fit <- estimate_rates(brts, model = case, link = "linear", method = "cem",
                        control = list(lower_bound = lbc, upper_bound = ubc, num_threads = 1L, max_time = 200))
  raw <- fit$details
  cat(sprintf("[%s] converged=%s iters=%d pars=(%s) loglik=%.3f ; exact at pars=%.3f ; best_IS ESS=%.1f\n",
              case, raw$converged, length(raw$best_loglik), paste(sprintf("%.4f", fit$pars), collapse = ", "),
              fit$loglik, exact(raw$obtained_estim), raw$best_IS$ESS))
  X <- replay_final(raw, brts, model, exact, R = 4L)
  if (case == "cr") cat(sprintf("  exact loglik at DDD MLE = %.3f ; obtained_estim is %.3f below it\n", ml$loglik, ml$loglik - exact(raw$obtained_estim)))
  cat("  elapsed:", round((proc.time() - t0)[3]), "s\n")
}

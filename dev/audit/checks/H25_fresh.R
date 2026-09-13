## H25 supplement: remove the winner's-curse confound. On the CR tree of H25.R
## (seed 31), take the final CEM population, then (a) R fresh 20-tree
## evaluations at obtained_estim, (b) R fresh 20-tree evaluations at the
## exact-best elite and at the noisy 20-tree winner. Compare mean fhat, its SD,
## and the exact loglik (DDD::bd_loglik, cond=0, btorph=1, soc=2) at each point.
## Also: E[max over elites of a noisy 20-tree fhat] - exact(winner) = selection bias.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressPackageStartupMessages({library(emphasis); library(DDD); library(ape)})
set.seed(31)
repeat { tr <- ape::rlineage(0.5, 0.3, Tmax = 6); tr <- ape::drop.fossil(tr); if (ape::Ntip(tr) >= 25 && ape::Ntip(tr) <= 40) break }
brts <- sort(as.numeric(ape::branching.times(tr)), decreasing = TRUE)
exact <- function(p8) DDD::bd_loglik(c(p8[1], p8[5], 0, 0), c(0, 0, 1, 0, 2), brts, 0)
eval20 <- function(p8, R = 10L) vapply(seq_len(R), function(.) {
  raw <- emphasis:::.simulate_particle(brts, p8, c(0L,0L,0L), 0L, 20L, 100L, 1e4, 1e6, 1L, 1.0)
  if (is.null(raw) || !length(raw$logf)) return(NA_real_)
  emphasis:::.is_fhat(raw$logf, raw$logg, n_zero_weight = emphasis:::.n0(raw$rejected_zero_weights))
}, numeric(1))
lb <- c(0.05, 0, 0, 0, 0, 0, 0, 0); ub <- c(2, 0, 0, 0, 1.5, 0, 0, 0)
fit <- emphasis:::emphasis_cem(brts, max_iter = 25L, num_points = 30L, max_missing = 1e4,
                       sd_vec = (ub - lb) / 4, lower_bound = lb, upper_bound = ub,
                       maxN = 10L, sample_size = 1L, disc_prop = 0.5, num_threads = 1L,
                       model = c(0L,0L,0L), link = 0L, verbose = FALSE, max_time = 120)
fp <- fit$final_pop; vix <- which(!is.na(fp$fhat)); ord <- vix[order(fp$fhat[vix], decreasing = TRUE)]
eix <- ord[seq_len(ceiling(0.5 * length(vix)))]; P <- as.matrix(fp$pars[eix, , drop = FALSE])
ex_el <- apply(P, 1, exact); ib <- which.max(ex_el)
est <- fit$obtained_estim
cat(sprintf("n_tips=%d  elites=%d  obtained_estim=(%.4f,%.4f)  exact(est)=%.3f  exact(best elite)=%.3f  best_IS$fhat=%.3f\n",
            length(brts)+1, nrow(P), est[1], est[5], exact(est), ex_el[ib], fit$best_IS$fhat))
a <- eval20(est); b <- eval20(P[ib, ])
cat(sprintf("fresh 20-tree fhat at obtained_estim : mean=%.3f sd=%.3f  (exact %.3f, mean-exact=%+.3f)\n", mean(a), sd(a), exact(est), mean(a)-exact(est)))
cat(sprintf("fresh 20-tree fhat at exact-best elite: mean=%.3f sd=%.3f  (exact %.3f, mean-exact=%+.3f)\n", mean(b), sd(b), ex_el[ib], mean(b)-ex_el[ib]))
cat(sprintf("difference of means (best elite - estim): %+.3f ; exact difference: %+.3f\n", mean(b)-mean(a), ex_el[ib]-exact(est)))
## winner's curse: one noisy draw per elite, take the max, compare with exact at the winner
wc <- t(replicate(8, { fh <- vapply(seq_len(nrow(P)), function(e) eval20(P[e, ], 1L), numeric(1)); w <- which.max(fh); c(max = fh[w], exact_w = ex_el[w], p = mean(exp(fh - max(fh)) / sum(exp(fh - max(fh))) > 0.5)) }))
cat(sprintf("winner's curse: mean[max_e fhat_e] - mean[exact(winner)] = %+.3f over 8 replays (per replay: %s)\n",
            mean(wc[, 1] - wc[, 2]), paste(sprintf("%+.2f", wc[, 1] - wc[, 2]), collapse = " ")))
cat(sprintf("replays with any softmax weight > 0.5: %d of 8\n", sum(wc[, 3] > 0)))

## H10 — BDI under DD: rejected draws (survivors at tp / max_missing overflow)
## make the effective proposal g/P(accept|theta), but -log P(accept) is not
## added to logg, so fhat is biased upward by -log P(accept|theta).
##
## Tests
##  A. CR control: same (lambda, mu), same target f.  model=cr (exact BDI,
##     acc = 1) vs model=dd with beta_N = gamma_N = 0 (mean-field Gillespie,
##     acc < 1).  fhat_dd - fhat_cr should equal -log(acc) if H10 holds, and
##     fhat_dd + log(acc) - fhat_cr should be ~ 0.
##  B. theta grid vs DDD::dd_loglik (ddmodel = 1, cond = 0, soc = 2):
##     gap_raw  = fhat_raw - dd_loglik  should vary with theta as -log(acc)
##     gap_corr = fhat_raw + log(acc) - dd_loglik  should be a constant
##  C. max_missing rejection channel + missing counter in .augment_tree_bdi.
##
## Everything is pure R (bdi.R) except eval_logf; C++ RNG is not involved in
## the BDI draws (stats::rexp / runif), so set.seed makes this reproducible.

.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
suppressMessages(library(DDD))
options(width = 140)

## ---- fixed tree (seeded through R's RNG) ---------------------------------
set.seed(11)
sim  <- DDD::dd_sim(pars = c(0.8, 0.3, 25), age = 6)
phy  <- sim$tes
brts <- sort(ape::branching.times(phy), decreasing = TRUE)
tp   <- brts[1]
bt   <- sort(tp - brts[-1])
cat(sprintf("Tree: %d tips, crown age %.3f\n\n", ape::Ntip(phy), tp))

lme <- function(x) { m <- max(x); log(mean(exp(x - m))) + m }
ess <- function(lw) { w <- exp(lw - max(lw)); sum(w)^2 / sum(w^2) }

## Draw n_attempt trees with the package's own sampler functions, record how
## many come back NULL (rejected), score the accepted ones with eval_logf
## exactly as .augment_tree_bdi does (bdi.R:684-699).
bdi_run <- function(pars8, model_bin, n_attempt, max_missing = 1e6L, link = 0L) {
  is_cr <- all(model_bin == 0L)
  p_fun <- Nhat_fun <- Phat_fun <- Ehat_fun <- NULL
  if (!is_cr) {
    sol <- emphasis:::.bdi_iterate(pars8, model_bin, link, bt, tp)
    p_fun <- sol$p_fun; Nhat_fun <- sol$Nhat_fun
    Phat_fun <- sol$Phat_fun; Ehat_fun <- sol$Ehat_fun
  }
  trees <- list(); logg <- numeric(0); n_rej <- 0L
  for (i in seq_len(n_attempt)) {
    a <- emphasis:::.bdi_augment_one(bt, pars8, model_bin, link, tp,
                                     p_fun, Nhat_fun, Phat_fun, Ehat_fun,
                                     as.integer(max_missing))
    if (is.null(a)) { n_rej <- n_rej + 1L; next }
    trees[[length(trees) + 1L]] <- emphasis:::.bdi_to_tree_df(a$species, bt, tp)
    logg <- c(logg, a$logg)
  }
  acc <- 1 - n_rej / n_attempt
  if (length(trees) == 0L) return(list(acc = acc, fhat_raw = -Inf, fhat_corr = -Inf, ess = 0, n_ok = 0L))
  ev   <- emphasis:::eval_logf(pars8, trees, model = as.integer(model_bin),
                               link = as.integer(link), rho = 1)
  lw   <- ev$logf - logg
  ok   <- is.finite(lw)
  list(acc = acc, n_ok = length(trees), n_nonfinite = sum(!ok),
       fhat_raw  = lme(lw[ok]),
       fhat_corr = lme(lw[ok]) + log(acc),
       ess = ess(lw[ok]), lw_sd = sd(lw[ok]), trees = trees, logf = ev$logf, logg = logg)
}

## =========================================================================
## A. CR control — identical target f, two proposals
## =========================================================================
cat("== A. CR control: cr (exact BDI) vs dd with beta_N = gamma_N = 0 ==\n")
lam <- 0.8; mu <- 0.3
p_cr <- c(lam, 0, 0, 0, mu, 0, 0, 0)
r_cr <- bdi_run(p_cr, c(0L, 0L, 0L), 300L)
cat(sprintf("cr : acc = %.3f  fhat = %.4f  sd(lw) = %.2e  (exact: lw constant)\n",
            r_cr$acc, r_cr$fhat_raw, r_cr$lw_sd))
## same trees scored under model_bin = dd give the same logf?
ev_dd_on_cr_trees <- emphasis:::eval_logf(p_cr, r_cr$trees[1:20], model = c(1L, 0L, 0L), link = 0L, rho = 1)
cat(sprintf("max |logf_dd - logf_cr| on the same 20 trees = %.2e  (target f identical)\n",
            max(abs(ev_dd_on_cr_trees$logf - r_cr$logf[1:20]))))
set.seed(101)
A <- t(sapply(1:5, function(r) {
  x <- bdi_run(p_cr, c(1L, 0L, 0L), 1500L)
  c(acc = x$acc, neg_log_acc = -log(x$acc), fhat_raw = x$fhat_raw,
    raw_minus_cr = x$fhat_raw - r_cr$fhat_raw,
    corr_minus_cr = x$fhat_corr - r_cr$fhat_raw, ess = x$ess)
}))
print(round(A, 4))
cat(sprintf("mean(raw - cr) = %.4f   mean(-log acc) = %.4f   mean(corr - cr) = %.4f  (se %.4f)\n\n",
            mean(A[, "raw_minus_cr"]), mean(A[, "neg_log_acc"]),
            mean(A[, "corr_minus_cr"]), sd(A[, "corr_minus_cr"]) / sqrt(nrow(A))))

## =========================================================================
## B. theta grid vs DDD::dd_loglik  (ddmodel 1: lambda(N) = l0 - (l0-mu0) N/K)
## =========================================================================
cat("== B. theta grid: BDI fhat vs DDD::dd_loglik(cond = 0, btorph = 1, soc = 2) ==\n")
grid <- expand.grid(K = c(12, 15, 20, 30, 50, 100, 1e4),
                    lm = c("0.8/0.3", "0.6/0.1", "1.2/0.3"), stringsAsFactors = FALSE)
lmv <- do.call(rbind, strsplit(grid$lm, "/")); grid$l0 <- as.numeric(lmv[, 1]); grid$m0 <- as.numeric(lmv[, 2])
n_rep <- 3L; n_att <- 1500L
set.seed(202)
rows <- lapply(seq_len(nrow(grid)), function(i) {
  l0 <- grid$l0[i]; m0 <- grid$m0[i]; K <- grid$K[i]
  pars8 <- c(l0, -(l0 - m0) / K, 0, 0, m0, 0, 0, 0)
  ref <- DDD::dd_loglik(pars1 = c(l0, m0, K), pars2 = c(300, 1, 0, 1, 0, 2),
                        brts = brts, missnumspec = 0)
  reps <- t(sapply(seq_len(n_rep), function(r) {
    x <- bdi_run(pars8, c(1L, 0L, 0L), n_att)
    c(acc = x$acc, fhat_raw = x$fhat_raw, fhat_corr = x$fhat_corr, ess = x$ess, nonfin = x$n_nonfinite)
  }))
  data.frame(l0 = l0, m0 = m0, K = K, dd_loglik = ref,
             acc = mean(reps[, "acc"]), neg_log_acc = -log(mean(reps[, "acc"])),
             fhat_raw = mean(reps[, "fhat_raw"]), sd_raw = sd(reps[, "fhat_raw"]),
             gap_raw = mean(reps[, "fhat_raw"]) - ref,
             gap_corr = mean(reps[, "fhat_corr"]) - ref,
             ess = mean(reps[, "ess"]), nonfin = sum(reps[, "nonfin"]))
})
B <- do.call(rbind, rows)
print(B, digits = 4, row.names = FALSE)
## Rows where lambda(N) = max(0, .) truncates (DDD -Inf or non-finite lw) are
## H11 territory (NaN / -Inf weights), not H10: exclude them from the summary.
Bc <- B[is.finite(B$dd_loglik) & B$nonfin == 0, ]
cat(sprintf("\nSummary over the %d clean grid points (no truncation, all lw finite):\n", nrow(Bc)))
cat(sprintf("  sd(gap_raw)  = %.4f   range = [%.4f, %.4f]\n", sd(Bc$gap_raw), min(Bc$gap_raw), max(Bc$gap_raw)))
cat(sprintf("  sd(gap_corr) = %.4f   range = [%.4f, %.4f]\n", sd(Bc$gap_corr), min(Bc$gap_corr), max(Bc$gap_corr)))
B <- Bc
fit <- lm(gap_raw ~ neg_log_acc, data = B)
cat(sprintf("  lm(gap_raw ~ -log acc): intercept = %.4f  slope = %.3f  R^2 = %.3f\n",
            coef(fit)[1], coef(fit)[2], summary(fit)$r.squared))
cat(sprintf("  cor(gap_raw, -log acc) = %.3f ; cor(gap_corr, -log acc) = %.3f\n\n",
            cor(B$gap_raw, B$neg_log_acc), cor(B$gap_corr, B$neg_log_acc)))

## =========================================================================
## C. max_missing rejection channel, and the missing counter
## =========================================================================
cat("== C. max_missing rejection (bdi.R:508) is the same uncorrected channel ==\n")
l0 <- 0.8; m0 <- 0.3; K <- 20
pars8 <- c(l0, -(l0 - m0) / K, 0, 0, m0, 0, 0, 0)
set.seed(303)
runsC <- lapply(c(1e6, 6, 4, 3, 2), function(mm) bdi_run(pars8, c(1L, 0L, 0L), 1500L, max_missing = mm))
C <- t(mapply(function(x, mm) c(max_missing = mm, acc = x$acc, neg_log_acc = -log(x$acc),
                                 fhat_raw = x$fhat_raw, fhat_corr = x$fhat_corr),
              runsC, c(1e6, 6, 4, 3, 2)))
print(round(C, 4))
cat("  (fhat_raw barely moves while 99.6% of proposal mass is cut: it estimates E_g[f/g | accept], not the truncated likelihood)\n")

cat("\n.augment_tree_bdi with low acceptance (K = 20, max_missing = 6, sample_size = 100): trees returned vs requested\n")
set.seed(404)
A2 <- emphasis:::.augment_tree_bdi(brts, pars8, c(1L, 0L, 0L), sample_size = 100L,
                                   max_missing = 6L, link = 0L, rho = 1)
cat(sprintf("  .augment_tree_bdi returned %d trees (max_tries = 5*100); fields: %s\n",
            length(A2$trees), paste(names(A2), collapse = ", ")))
cat("  -> no acceptance-rate / n_rejected field; .mcem_bdi hard-codes rejected = 0L (bdi.R:762, 826)\n")

## =========================================================================
## D. M-step invariance: weights are self-normalised (bdi.R:796-799)
## =========================================================================
cat("\n== D. M-step weights are invariant to the missing -log P(accept) ==\n")
x0 <- runsC[[1]]; lw <- x0$logf - x0$logg; lw <- lw[is.finite(lw)]
wn <- function(v) { w <- exp(v - max(v)); w / sum(w) * length(w) }
cat(sprintf("  max |w_norm(lw) - w_norm(lw + log acc)| = %.2e   (P(accept) is common to all draws at theta_k)\n",
            max(abs(wn(lw) - wn(lw + log(x0$acc))))))

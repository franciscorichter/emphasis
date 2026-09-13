## H15 — bias_correct reaches only Mode-1 search-time fhat; property test of the
## moment-based estimator: ell_hat + log(1 + sum_{k=2}^K m_k/k!) -> log mean exp(lw) as K -> Inf.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
set.seed(1)
lw <- rnorm(200, sd = 1.2)
lme <- log(mean(exp(lw)))
f0 <- emphasis:::.is_fhat(lw, rep(0, 200), bias_correct = FALSE)
f2 <- emphasis:::.is_fhat(lw, rep(0, 200), bias_correct = TRUE, K = 2L)
cat(sprintf("log-mean-exp      = %.6f\n", lme))
cat(sprintf("bias_correct=F    = %.6f (diff %.2e)\n", f0, f0 - lme))
cat(sprintf("bias_correct=T K2 = %.6f (diff %.2e)  -> changes fhat: %s\n", f2, f2 - lme, f2 != f0))
for (K in c(2, 4, 8, 16, 32, 64)) {
  fk <- emphasis:::.is_fhat(lw, rep(0, 200), bias_correct = TRUE, K = as.integer(K))
  cat(sprintf("  K=%3d: %.8f  |diff| = %.2e\n", K, fk, abs(fk - lme)))
}
## with zero-weight trees the branch silently falls back
fz <- emphasis:::.is_fhat(lw, rep(0, 200), bias_correct = TRUE, K = 2L, n_zero_weight = 3L)
fz0 <- emphasis:::.is_fhat(lw, rep(0, 200), bias_correct = FALSE, n_zero_weight = 3L)
cat(sprintf("n_zero_weight=3: bias_correct=T gives %.6f, =F gives %.6f, identical: %s\n", fz, fz0, identical(fz, fz0)))
## call-site audit: where is bias_correct passed?
src <- readLines("/Users/pancho/Code/emphasis/R/de.R")
hits <- grep("\\.is_fhat\\(", src)
for (h in hits) cat(sprintf("de.R:%d  %s\n", h, trimws(paste(src[h:(h+2)], collapse = " "))))

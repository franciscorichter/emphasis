.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
# 10^4 augmentations at rho = 1: does any augmented node carry t_ext == 5e10 (unsampled sentinel)?
set.seed(1)
phy <- ape::rlineage(0.5, 0.0, Tmax = 6)  # not used; build brts directly
brts <- sort(c(6, runif(19, 0, 6)), decreasing = TRUE)  # 20 tips, crown age 6
for (pars in list(c(0.6,0,0,0, 0.3,0,0,0), c(0.6,0,0,0, 1e-3,0,0,0), c(0.6,0,0,0, 3,0,0,0))) {
  t0 <- proc.time()[3]
  raw <- emphasis:::augment_trees(brts, pars, sample_size = 10000L, maxN = 100000L,
                                   max_missing = 10000L, max_lambda = 1e4, num_threads = 1L,
                                   model = c(0L,0L,0L), link = 0L, rho = 1.0)
  n_unsamp <- sum(vapply(raw$trees, function(df) sum(df$t_ext == 5e10), 0))
  n_ge_T   <- sum(vapply(raw$trees, function(df) sum(df$t_ext != 0 & df$t_ext < 1e11 & df$t_ext >= brts[1]), 0))
  n_aug    <- sum(vapply(raw$trees, function(df) sum(df$t_ext == 0), 0))
  cat(sprintf("mu=%g: trees=%d aug_extinct=%d  t_ext==5e10: %d  t_ext>=T (non-sentinel): %d  nonfinite lw: %d  time=%.1fs\n",
      pars[5], length(raw$trees), n_aug, n_unsamp, n_ge_T, sum(!is.finite(raw$logf - raw$logg)), proc.time()[3]-t0))
}

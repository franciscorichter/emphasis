.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
# Timing of augment_trees when gamma0 = 0 under linear link (mu clamped to 1e-10 inside C++)
set.seed(2)
brts <- sort(c(6, runif(19, 0, 6)), decreasing = TRUE)
run <- function(pars, rho, n = 2000L) {
  t0 <- proc.time()[3]
  raw <- emphasis:::augment_trees(brts, pars, sample_size = n, maxN = 20L*n,
                                   max_missing = 10000L, max_lambda = 1e4, num_threads = 1L,
                                   model = c(0L,0L,0L), link = 0L, rho = rho)
  el <- proc.time()[3]-t0
  n_aug <- mean(vapply(raw$trees, function(df) sum(df$t_ext == 0), 0))
  n_uns <- mean(vapply(raw$trees, function(df) sum(df$t_ext == 5e10), 0))
  cat(sprintf("gamma0=%g rho=%g: %d trees in %.2fs (%.1f ms/tree); mean extinct/tree=%.3f, mean unsampled/tree=%.2f\n",
      pars[5], rho, length(raw$trees), el, 1000*el/max(1,length(raw$trees)), n_aug, n_uns))
}
run(c(0.6,0,0,0, 0,0,0,0), 1.0)
run(c(0.6,0,0,0, 0.3,0,0,0), 1.0)
run(c(0.6,0,0,0, 1e-6,0,0,0), 1.0)
run(c(0.6,0,0,0, 0,0,0,0), 0.8)
run(c(0.6,0,0,0, 0.3,0,0,0), 0.8)
# also negative eta (max(0, eta) = 0 -> clamped)
run(c(0.6,0,0,0, -2,0,0,0), 1.0)

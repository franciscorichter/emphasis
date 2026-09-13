.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
# Scaling of augmentation time with number of augmented nodes (CR model, no PD/EP covariates)
set.seed(3)
brts <- sort(c(10, runif(39, 0, 10)), decreasing = TRUE)   # 40 tips
out <- data.frame()
for (mu in c(0.05, 0.3, 0.6, 0.9, 1.2, 1.5)) {
  pars <- c(0.5,0,0,0, mu,0,0,0)
  t0 <- proc.time()[3]
  raw <- emphasis:::augment_trees(brts, pars, sample_size = 20L, maxN = 400L,
                                   max_missing = 20000L, max_lambda = 1e6, num_threads = 1L,
                                   model = c(0L,0L,0L), link = 0L, rho = 1.0)
  el <- proc.time()[3]-t0
  nn <- vapply(raw$trees, nrow, 0L)
  out <- rbind(out, data.frame(mu = mu, trees = length(raw$trees), mean_nodes = mean(nn), max_nodes = max(nn), sec = el, ms_per_tree = 1000*el/max(1,length(raw$trees))))
  if (el > 120) break
}
print(out)
fit <- lm(log(ms_per_tree) ~ log(mean_nodes), data = out[out$mean_nodes > 100, ])
cat("empirical exponent of time vs nodes:", coef(fit)[2], "\n")

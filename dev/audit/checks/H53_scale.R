.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
# larger trees: does per-tree cost scale superlinearly in node count (CR model, no covariates)?
set.seed(6)
out <- data.frame()
for (ntip in c(25, 50, 100, 200)) {
  brts <- sort(c(10, runif(ntip - 1, 0, 10)), decreasing = TRUE)
  for (mu in c(0.1, 0.3)) {
    pars <- c(0.6,0,0,0, mu,0,0,0)
    t0 <- proc.time()[3]
    raw <- tryCatch(emphasis:::augment_trees(brts, pars, sample_size = 10L, maxN = 5000L,
                                     max_missing = 50000L, max_lambda = 1e7, num_threads = 1L,
                                     model = c(0L,0L,0L), link = 0L, rho = 1.0), error = function(e) {message(conditionMessage(e)); NULL}); if (is.null(raw)) next
    el <- proc.time()[3]-t0
    nn <- vapply(raw$trees, nrow, 0L)
    out <- rbind(out, data.frame(ntip = ntip, mu = mu, trees = length(raw$trees), mean_nodes = mean(nn), ms_per_tree = 1000*el/max(1,length(raw$trees))))
  }
}
print(out)
fit <- lm(log(ms_per_tree) ~ log(mean_nodes), data = out)
cat("empirical exponent (all points):", coef(fit)[2], "\n")

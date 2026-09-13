.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(3)
for (k in 1:3) {
  sim <- NULL; while (is.null(sim$tes) || ape::Ntip(sim$tes) < 10 || ape::Ntip(sim$tes) > 40)
    sim <- simulate_tree(pars = c(0.6, 0.1), max_t = 5, model = "cr")
  tes <- sim$tes; brts <- emphasis:::.extract_brts(sim)
  tip_edges <- tes$edge[, 2] <= ape::Ntip(tes)
  true_mean_pendant <- mean(tes$edge.length[tip_edges])
  oc <- emphasis:::.observed_covariates(brts, c(0L, 1L, 0L))
  cat(sprintf("tips=%d  X_obs(M)=%.4f  mean(brts)=%.4f  true mean pendant=%.4f  ratio=%.2f\n",
      ape::Ntip(tes), oc[[1]]$X_obs, mean(brts), true_mean_pendant, oc[[1]]$X_obs/true_mean_pendant))
}

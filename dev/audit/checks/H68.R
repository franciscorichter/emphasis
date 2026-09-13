.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(3); tr <- ape::drop.fossil(ape::rlineage(0.5, 0.2, Tmax = 8))
while (ape::Ntip(tr) < 15 || ape::Ntip(tr) > 40) tr <- ape::drop.fossil(ape::rlineage(0.5,0.2,Tmax=8))
brts <- sort(ape::branching.times(tr), decreasing = TRUE); cat("tips", ape::Ntip(tr), "\n")
# DD linear link with slope that drives lambda to 0 -> many zero-weight rejections
model <- c(1L,0L,0L)
ip <- emphasis:::.expand_pars(c(0.8, -0.03, 0.2, 0), model)
N <- 100L
stats <- function(nt, reps = 15) {
  out <- t(sapply(seq_len(reps), function(i) {
    a <- emphasis:::augment_trees(brts, ip, N, 5000L, 1e4, 1e6, nt, model, 0L, 1.0)
    c(ntrees = length(a$trees), zero = a$rejected_zero_weights, over = a$rejected_overruns, lam = a$rejected_lambda, rej = a$rejected)
  }))
  cat(sprintf("threads=%d: ntrees range %s; mean zero=%.1f (sd %.1f), overruns=%.1f\n", nt,
      paste(range(out[,"ntrees"]), collapse="-"), mean(out[,"zero"]), sd(out[,"zero"]), mean(out[,"over"])))
}
stats(1); stats(8)

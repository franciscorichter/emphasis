.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths())); library(emphasis)
set.seed(3); tr <- ape::drop.fossil(ape::rlineage(0.5, 0.2, Tmax = 8))
while (ape::Ntip(tr) < 15 || ape::Ntip(tr) > 40) tr <- ape::drop.fossil(ape::rlineage(0.5,0.2,Tmax=8))
brts <- sort(ape::branching.times(tr), decreasing = TRUE); cat("tips", ape::Ntip(tr), "\n")
model <- c(1L,0L,0L); ip <- emphasis:::.expand_pars(c(0.8, -0.03, 0.2, 0), model)
lb <- emphasis:::.expand_pars(c(0,-0.5,0,-0.5), model); ub <- emphasis:::.expand_pars(c(2,0.5,1,0.5), model)
N <- 100L; maxN <- 3000L
for (nt in c(1L, 4L, 8L)) {
  out <- t(sapply(1:10, function(i) { r <- emphasis:::em_cpp(brts, ip, N, maxN, 1e4, 1e6, lb, ub, 1e-3, nt, FALSE, model, 0L, 1.0, NULL); c(trees = r$trees, fhat = r$fhat, zero = r$rejected_zero_weights, ms = r$time) }))
  cat(sprintf("threads=%d: trees returned: %s | fhat mean %.2f (sd %.2f) | zero-w mean %.1f | time ms mean %.0f\n", nt,
              paste(out[,"trees"], collapse=","), mean(out[,"fhat"]), sd(out[,"fhat"]), mean(out[,"zero"]), mean(out[,"ms"])))
}

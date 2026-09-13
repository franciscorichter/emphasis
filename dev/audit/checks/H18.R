## H18 — S = N + rejected_zero_weights excludes overrun/lambda rejections.
## CR tree, low max_missing to force overruns; compare fhat (as coded) and
## fhat_all = log(sum_w / (N + zero + overrun + lambda)) + max against DDD::bd_loglik
## over a mu grid.  Gap constancy across theta is the criterion.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis)); suppressMessages(library(DDD))
set.seed(7)
phy <- ape::drop.fossil(ape::rlineage(0.5, 0.2, Tmax = 5))
while (ape::Ntip(phy) < 10 || ape::Ntip(phy) > 25) phy <- ape::drop.fossil(ape::rlineage(0.5, 0.2, Tmax = 5))
brts <- sort(ape::branching.times(phy), decreasing = TRUE)
cat(sprintf("tree: %d tips crown %.2f\n", ape::Ntip(phy), brts[1]))
lambda <- 0.5
one <- function(mu, max_missing, N = 400L, reps = 4) {
  out <- t(replicate(reps, {
    raw <- emphasis:::augment_trees(as.numeric(brts), c(lambda,0,0,0,mu,0,0,0), sample_size = N, maxN = 200000L,
                         max_missing = max_missing, max_lambda = 1e4, num_threads = 1L,
                         model = c(0L,0L,0L), link = 0L, rho = 1)
    lw <- raw$logf - raw$logg; m <- max(lw); sw <- sum(exp(lw - m))
    S_code <- N + raw$rejected_zero_weights
    S_all  <- S_code + raw$rejected_overruns + raw$rejected_lambda + raw$rejected
    c(fhat_code = log(sw / S_code) + m, fhat_all = log(sw / S_all) + m,
      over = raw$rejected_overruns, zero = raw$rejected_zero_weights, lam = raw$rejected_lambda,
      ess = emphasis:::.ess_from_lw(lw))
  }))
  colMeans(out)
}
ref <- function(mu) DDD::bd_loglik(pars1 = c(lambda, mu, 0, 0), pars2 = c(0, 0, 1, 0, 2), brts = brts, missnumspec = 0)
mus <- c(0.05, 0.15, 0.25, 0.35, 0.45)
for (mm in c(1000L, 8L, 5L)) {
  cat(sprintf("\nmax_missing = %d\n", mm))
  res <- t(sapply(mus, function(mu) { r <- one(mu, mm); c(mu = mu, r, ref = ref(mu), gap_code = unname(r["fhat_code"]) - ref(mu), gap_all = unname(r["fhat_all"]) - ref(mu)) }))
  print(round(res, 3))
  cat(sprintf("  sd(gap_code) = %.3f   sd(gap_all) = %.3f\n", sd(res[, "gap_code"]), sd(res[, "gap_all"])))
}
## inverse-binomial bias at sample_size = 1 (pure R, analytic form): E[w/(1+G)] vs p*w
p <- 0.3; G <- rgeom(2e5, p)
cat(sprintf("\nN=1 stop-at-first-success: E[1/(1+G)] = %.4f vs p = %.4f (ratio %.3f)\n", mean(1/(1+G)), p, mean(1/(1+G))/p))

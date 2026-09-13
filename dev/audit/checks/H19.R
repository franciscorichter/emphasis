## H19 — rho terms: n_obs log rho + n_unsamp log(1-rho), no binomial coefficient.
## Test: fhat(rho) from the thinning sampler vs Stadler (2009) constant-rate
## reconstructed-tree density with Bernoulli sampling; the gap must be rho-independent.
.libPaths(c("/Users/pancho/.claude/jobs/867af780/tmp/rlib", .libPaths()))
suppressMessages(library(emphasis))
set.seed(9)
phy <- ape::drop.fossil(ape::rlineage(0.5, 0.15, Tmax = 4))
while (ape::Ntip(phy) < 6 || ape::Ntip(phy) > 14) phy <- ape::drop.fossil(ape::rlineage(0.5, 0.15, Tmax = 4))
brts <- sort(ape::branching.times(phy), decreasing = TRUE)   # x_1 = crown age, then internal nodes
n <- ape::Ntip(phy)
cat(sprintf("tree: %d tips crown %.2f\n", n, brts[1]))
la <- 0.5; mu <- 0.15
p1 <- function(t, rho) { r <- la - mu; rho * r^2 * exp(-r*t) / (rho*la + (la*(1-rho) - mu)*exp(-r*t))^2 }
stadler <- function(rho) (n - 2) * log(la) + 2 * log(p1(brts[1], rho)) + sum(log(p1(brts[-1], rho)))
one <- function(rho, N = 1500L, reps = 3) {
  v <- replicate(reps, {
    raw <- emphasis:::augment_trees(as.numeric(brts), c(la,0,0,0,mu,0,0,0), sample_size = N, maxN = 300000L,
                         max_missing = 500L, max_lambda = 1e4, num_threads = 1L,
                         model = c(0L,0L,0L), link = 0L, rho = rho)
    nun <- vapply(raw$trees, function(tr) sum(tr$t_ext == 5e10), integer(1))
    c(fhat = emphasis:::.is_fhat(raw$logf, raw$logg, n_zero_weight = raw$rejected_zero_weights), ess = emphasis:::.ess_from_lw(raw$logf - raw$logg), mean_unsamp = mean(nun), zero = raw$rejected_zero_weights)
  })
  c(rowMeans(v), sd_fhat = sd(v["fhat", ]))
}
rhos <- c(1, 0.8, 0.6, 0.4, 0.25)
res <- t(sapply(rhos, function(r) { o <- one(r); c(rho = r, o, stadler = stadler(r), gap = unname(o["fhat"]) - stadler(r)) }))
print(round(res, 3))
cat(sprintf("gap range = %.3f, sd = %.3f (Monte Carlo sd of fhat ~ %.3f)\n", diff(range(res[, "gap"])), sd(res[, "gap"]), mean(res[, "sd_fhat"])))
## what a binomial coefficient would add at each rho (mean over sampled augmentations)
cat("If choose(n_obs+n_unsamp, n_obs) were added, lchoose at mean n_unsamp:\n")
print(round(cbind(rho = rhos, lchoose = lchoose(n + res[, "mean_unsamp"], n)), 3))
